#Fit the localisation + migration model to both mice, with tree uncertainty, simulation validation and posterior
#predictive checks. Usage: Rscript migration_model/run_migration.R [n_boot_trees] [n_sim_reps]
setwd("~/R_work/mouse_phylo/")
suppressMessages({library(ggplot2);library(tidyr);library(parallel);library(phytools)})
for(f in list.files("~/R_work/my_functions/",pattern=".R",full.names=T)) suppressMessages(source(f))
source("migration_model/R/migration_model.R")
args=commandArgs(TRUE); n_boot=ifelse(length(args)>0,as.integer(args[1]),100); n_rep=ifelse(length(args)>1,as.integer(args[2]),40)
out_dir="migration_model/results"; dir.create(out_dir,showWarnings=F)
n_cores=max(1,detectCores()-2); RNGkind("L'Ecuyer-CMRG"); set.seed(20261030) #reproducible mclapply streams

md=suppressMessages(readxl::read_excel("sample_metadata/CAMPUOY000001_REUPLOAD_20231116_4xplates_York_Lily_CASM_STS_Sample_MANIFEST_Rec20231206_WithCGPIDs20231219.xlsx"))
tissue=setNames(md$TISSUE_PHENOTYPE,md$INTERNAL_CASM_SAMPLE_NAME)
bones=c("Femur","Iliac","Spine","Tibia")
sampling_age=0.5 #years (SPECIMEN_AGE_YEARS)
mu_grid=c(0,10^seq(-3.5,0.5,by=0.05))
prior_mu=c(0.5,rep(0.5/(length(mu_grid)-1),length(mu_grid)-1)) #half on "no migration", half log-uniform on mu > 0

#Parametric bootstrap of the corrected tree: Poisson-resample branch lengths, then make ultrametric again
boot_tree=function(tree){tree$edge.length=rpois(length(tree$edge.length),tree$edge.length); tree}

fits=list(); summ=list()
for(m in c("MD7634","MD7635")){
  bc=readRDS(paste0("burden_correction/results/burden_correction_",m,".Rds"))
  u=make_ultra(bc$tree); pt=prep_tree(u,tissue[u$tip.label],levels=bones)
  T_grid=seq(1,floor(pt$D),by=1)
  #Calibrated localisation window: bone-marrow colonisation ~E15.5 to ~P7 (Kapadia et al. 2025 molecular clock), either
  #taking our burdens at face value or scaling our trees so the 6-month depth matches theirs (48.2 + 45.3/2 = 70.9)
  scale=pt$D/age_to_mol(sampling_age)
  win=c(-4/365,7/365)
  priors_T=list(T_1_to_D=rep(1,length(T_grid)),T_5_to_50=as.numeric(T_grid>=5&T_grid<=50),
                calibrated_unscaled=as.numeric(T_grid>=age_to_mol(win[1])&T_grid<=age_to_mol(win[2])),
                calibrated_scaled=as.numeric(T_grid>=age_to_mol(win[1],scale)&T_grid<=age_to_mol(win[2],scale)))
  years_after=list(calibrated_unscaled=function(T) sampling_age-mol_to_age(T),calibrated_scaled=function(T) sampling_age-mol_to_age(T,scale))
  prior_T=priors_T$T_5_to_50
  cat(m,": depth",round(pt$D,1),"| scale",round(scale,3),"| calibrated window unscaled",round(age_to_mol(win),1),"scaled",round(age_to_mol(win,scale),1),"\n")
  ll=loglik_grid(pt,T_grid,mu_grid)
  #Tree uncertainty: average the likelihood over bootstrap trees (heights rescaled to the point-estimate depth)
  ll_boot=mclapply(1:n_boot,function(b){ub=make_ultra(boot_tree(bc$tree)); ub$edge.length=ub$edge.length*pt$D/max(node.depth.edgelength(ub))
    loglik_grid(prep_tree(ub,tissue[ub$tip.label],levels=bones),T_grid,mu_grid,pi=pt$pi_obs)},mc.cores=n_cores)
  ll_avg=log(Reduce(`+`,lapply(ll_boot,function(x) exp(x-max(ll)))) /n_boot)+max(ll)
  fits[[m]]=list(pt=pt,T_grid=T_grid,ll=ll,ll_boot_avg=ll_avg,prior_T=prior_T,priors_T=priors_T,scale=scale)
  ll_indep=sum(log(pt$pi_obs[as.integer(pt$states)]))  #bone independent of the phylogeny
  cat(m,": max logL",round(max(ll),2),"| no clustering",round(ll_indep,2),"| best with no migration",round(max(ll[,1]),2),
      "at T =",T_grid[which.max(ll[,1])],"\n")
  for(v in c("point_tree","tree_uncertainty")) for(pr in names(priors_T)){
    L=if(v=="point_tree") ll else ll_avg
    post=grid_posterior(L,priors_T[[pr]],prior_mu)
    s=post_summary(post,T_grid,mu_grid,pt$D,pt$pi_obs,years_after=years_after[[pr]])
    summ[[paste(m,v,pr)]]=data.frame(mouse=m,tree=v,prior_T=pr,quantity=rownames(s),s,P_no_migration=sum(post[,mu_grid==0]),row.names=NULL)
  }
}
summ=bind_rows(summ); write.csv(summ,file.path(out_dir,"posterior_summary.csv"),row.names=F)
print(summ%>%mutate(across(where(is.numeric),~signif(.x,3))),row.names=F)

#Joint fit (shared T as a fraction of depth, shared mu): sum log-likelihoods on a common relative-T grid
#Priors: 5-50 molecular time (at the mean depth), and the calibrated window in scaled (Kapadia-equivalent) units
rel=seq(0.01,0.99,by=0.005)
ll_joint=Reduce(`+`,lapply(fits,function(f) loglik_grid(f$pt,rel*f$pt$D,mu_grid)))
D_mean=mean(sapply(fits,function(f) f$pt$D)); pi_mean=colMeans(do.call(rbind,lapply(fits,function(f) f$pt$pi_obs)))
D_ref=age_to_mol(sampling_age)
joint_priors=list(T_5_to_50=as.numeric(rel*D_mean>=5&rel*D_mean<=50),
                  calibrated_scaled=as.numeric(rel*D_ref>=age_to_mol(-4/365)&rel*D_ref<=age_to_mol(7/365)))
joint=bind_rows(lapply(names(joint_priors),function(pr){post_joint=grid_posterior(ll_joint,joint_priors[[pr]],prior_mu)
  #T reported in Kapadia-equivalent molecular time (rel x 70.9)
  s=post_summary(post_joint,rel*D_ref,mu_grid,D_ref,pi_mean,years_after=function(T) sampling_age-mol_to_age(T))
  data.frame(mouse="joint",tree="point_tree",prior_T=pr,quantity=rownames(s),s,P_no_migration=sum(post_joint[,mu_grid==0]),row.names=NULL)}))
print(joint%>%mutate(across(where(is.numeric),~signif(.x,3))),row.names=F)
write.csv(joint,file.path(out_dir,"posterior_summary_joint.csv"),row.names=F)

#Conditional on T: posterior of the expected number of moves per lineage after localisation, for a range of fixed T
cond=bind_rows(lapply(names(fits),function(m){f=fits[[m]]; meff=mu_grid*(1-sum(f$pt$pi_obs^2))
  bind_rows(lapply(c(5,10,15,20,25,30,35,40,45,50,55,60),function(Tf){i=which.min(abs(f$T_grid-Tf)); p=exp(f$ll_boot_avg[i,]-max(f$ll_boot_avg[i,]))*prior_mu; p=p/sum(p)
    n_mig=meff*(f$pt$D-Tf); o=order(n_mig); cw=cumsum(p[o])
    data.frame(mouse=m,T=Tf,age_unscaled_days=round(365*mol_to_age(Tf)),age_scaled_days=round(365*mol_to_age(Tf,f$scale)),P_no_migration=p[1],moves_q2.5=n_mig[o][which(cw>=0.025)[1]],moves_median=n_mig[o][which(cw>=0.5)[1]],
               moves_q97.5=n_mig[o][which(cw>=0.975)[1]],max_loglik=max(f$ll_boot_avg[i,]))}))}))
print(cond%>%mutate(across(where(is.numeric),~signif(.x,3))),row.names=F)
write.csv(cond,file.path(out_dir,"moves_given_T.csv"),row.names=F)

#Likelihood surface plots
surf=bind_rows(lapply(names(fits),function(m){f=fits[[m]]; d=as.data.frame(as.table(f$ll_boot_avg)); names(d)=c("T","mu","ll")
  d%>%mutate(mouse=m,T=as.numeric(as.character(T)),mu=as.numeric(as.character(mu)),moves=mu*(1-sum(f$pt$pi_obs^2))*(f$pt$D-T),dll=ll-max(ll))}))
pdf(file.path(out_dir,"migration_model_plots.pdf"),width=9,height=4)
print(ggplot(surf%>%filter(mu>0)%>%mutate(dll=ifelse(is.finite(dll),pmax(dll,-10),-10)),aes(T,mu,fill=dll))+geom_tile()+scale_y_log10()+
  facet_wrap(~mouse,scales="free_x")+geom_contour(aes(z=dll),breaks=c(-3,-1.92),col="white")+scale_fill_viridis_c()+theme_classic()+
  labs(x="Localisation time T (molecular time)",y="Migration event rate mu (per unit molecular time)",fill="Delta log L",
       title="Likelihood surface (averaged over bootstrap trees); contours at -1.92 and -3"))
print(ggplot(cond,aes(T,moves_median,ymin=pmax(moves_q2.5,1e-3),ymax=moves_q97.5))+geom_pointrange()+facet_wrap(~mouse)+scale_y_log10()+
  theme_classic()+labs(x="Localisation time T (fixed)",y="Moves per lineage after T (median, 95% CrI)"))
dev.off()

#---------------- Validation: simulate tip bones on the real MD7634 tree with known T and moves ----------------
f=fits$MD7634; pt=f$pt; meff_fac=1-sum(pt$pi_obs^2)
truth=expand.grid(T=c(10,20,35),moves=c(0,0.1,0.3,1))
val=bind_rows(mclapply(seq_len(nrow(truth)*n_rep),function(k){i=(k-1)%/%n_rep+1; Tt=truth$T[i]; mv=truth$moves[i]
  mu=mv/meff_fac/(pt$D-Tt)
  s=simulate_states(pt,Tt,mu); ps=prep_tree(pt$tree,s,levels=bones)
  post=grid_posterior(loglik_grid(ps,f$T_grid,mu_grid),f$prior_T,prior_mu)
  q=post_summary(post,f$T_grid,mu_grid,pt$D,ps$pi_obs)
  data.frame(T_true=Tt,moves_true=mv,T_med=q["T","q50"],T_lo=q["T","q2.5"],T_hi=q["T","q97.5"],
             moves_med=q["moves_per_lineage","q50"],moves_lo=q["moves_per_lineage","q2.5"],moves_hi=q["moves_per_lineage","q97.5"],
             P_no_migration=sum(post[,mu_grid==0]))},mc.cores=n_cores))
write.csv(val,file.path(out_dir,"validation_simulations.csv"),row.names=F)
val_summ=val%>%group_by(T_true,moves_true)%>%summarise(T_median=median(T_med),T_cover95=mean(T_lo<=T_true&T_hi>=T_true),
  moves_median=median(moves_med),moves_cover95=mean(moves_lo<=moves_true+1e-9&moves_hi>=moves_true-1e-9),
  mean_P_no_migration=mean(P_no_migration),.groups="drop")
print(as.data.frame(val_summ%>%mutate(across(where(is.numeric),~signif(.x,3)))),row.names=F)
write.csv(val_summ,file.path(out_dir,"validation_summary.csv"),row.names=F)

#---------------- Posterior predictive check: same-bone fraction among pairs coalescing after height h ----------------
ppc=list()
for(m in names(fits)){f=fits[[m]]; pt=f$pt; h_cut=seq(0,floor(pt$D)-5,by=2)
  post=grid_posterior(f$ll_boot_avg,f$prior_T,prior_mu); idx=sample(length(post),200,replace=T,prob=post)
  draws=data.frame(T=f$T_grid[row(post)[idx]],mu=mu_grid[col(post)[idx]])
  #Comparison models: best T with no migration, and early localisation (T=5) with no migration
  T0=f$T_grid[which.max(f$ll_boot_avg[,1])]
  models=list(posterior=draws,no_migration_best_T=data.frame(T=T0,mu=0)[rep(1,200),],no_migration_T5=data.frame(T=5,mu=0)[rep(1,200),],
              no_localisation=data.frame(T=min(pt$h[1:pt$n_tip])-1e-6,mu=0)[rep(1,200),])
  obs=pair_stat(pt,setNames(as.character(pt$states),pt$tree$tip.label),h_cut)
  for(nm in names(models)){d=models[[nm]]
    sims=do.call(rbind,mclapply(1:nrow(d),function(r) pair_stat(pt,simulate_states(pt,d$T[r],d$mu[r]),h_cut),mc.cores=n_cores))
    ppc[[paste(m,nm)]]=data.frame(mouse=m,model=nm,h=h_cut,obs=obs,lo=apply(sims,2,quantile,0.025,na.rm=T),
                                  med=apply(sims,2,median,na.rm=T),hi=apply(sims,2,quantile,0.975,na.rm=T))}
  cat(m,": best T with no migration =",T0,"\n")
}
ppc=bind_rows(ppc); write.csv(ppc,file.path(out_dir,"posterior_predictive.csv"),row.names=F)
pdf(file.path(out_dir,"posterior_predictive.pdf"),width=9,height=6)
print(ggplot(ppc,aes(h))+geom_ribbon(aes(ymin=lo,ymax=hi,fill=model),alpha=0.25)+geom_line(aes(y=med,col=model))+
  geom_point(aes(y=obs),size=0.8)+facet_wrap(~mouse)+theme_classic()+
  labs(x="Molecular time h (pairs whose common ancestor is at height >= h)",y="Fraction of pairs from the same bone",
       title="Observed (points) vs simulated (median, 95% interval)"))
dev.off()
saveRDS(list(fits=fits,mu_grid=mu_grid,prior_mu=prior_mu,summary=summ,joint=joint,conditional=cond,validation=val,ppc=ppc),
        file.path(out_dir,"migration_model_fits.Rds"))
