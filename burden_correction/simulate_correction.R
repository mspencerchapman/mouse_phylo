#Simulation test of the burden correction: colonies with known clonal / sub-clonal counts at the real per-sample
#depths of each mouse, the detection rule and partial detection of 3-read calls; compare estimators with the truth
#Usage: Rscript simulate_correction.R [n_reps]
args=commandArgs(trailingOnly=TRUE); n_reps=ifelse(length(args)>0,as.integer(args[1]),2)
setwd("~/R_work/mouse_phylo/")
source("burden_correction/R/burden_correction.R")
out_dir="burden_correction/results/simulation"; dir.create(out_dir,recursive=T,showWarnings=F)

#Rmd method: pooled 2-component binomial mixture at NR > 20 gives the clonal prior; hard calls at VAF 0.5 vs 0.33
rmd_method=function(priv){
  hi=priv[priv$NR>20,]; p=c(0.5,0.3); pr=c(0.5,0.5)
  for(it in 1:500){a=pr[1]*dbinom(hi$NV,hi$NR,p[1]); b=pr[2]*dbinom(hi$NV,hi$NR,p[2]); g=a/(a+b)
    pr=c(mean(g),1-mean(g)); p=c(sum(g*hi$NV/hi$NR)/sum(g),sum((1-g)*hi$NV/hi$NR)/sum(1-g))}
  prior=pr[which.min(abs(p-0.5))]
  call=prior*dbinom(priv$NV,priv$NR,0.5)>(1-prior)*dbinom(priv$NV,priv$NR,0.33)
  tapply(call,priv$sample,sum)
}

#Simulate the detected private mutations of all samples
#vaf_sub: function(n) giving true sub-clonal VAFs; rho: beta-binomial overdispersion; pc_sd: colony-to-colony clonal VAF sd (0: single-cell derived)
simulate_mouse=function(depths,truth,vaf_sub,delta=c(0.72,1),rho=0,pc_sd=0,n_artefact=0,xy_share=0.033){
  lapply(names(depths),function(s){
    t=truth[truth$sample==s,]; pc_s=min(max(rnorm(1,0.5,pc_sd),0.3),0.6)
    vaf=c(rep(pc_s,t$N_c),vaf_sub(t$N_s),runif(n_artefact,0.08,0.2)); n=length(vaf)
    xy=runif(n)<xy_share & length(depths[[s]]$xy)>0
    NR=ifelse(xy,sample(depths[[s]]$xy,n,replace=T),sample(depths[[s]]$auto,n,replace=T))
    p=pmin(ifelse(xy,2,1)*vaf,0.999)
    if(rho>0){a=p*(1-rho)/rho; b=(1-p)*(1-rho)/rho; p=rbeta(n,a,b)}
    NV=rbinom(n,NR,p)
    keep=detect_rule(NV,NR) & runif(n)<detect_weight(NV,delta)
    data.frame(sample=s,NV=NV,NR=NR,cf=ifelse(xy,2,1))[keep,]
  })%>%bind_rows()
}

scenarios=list(
  matched=list(vaf_sub=function(n) rep(0.27,n)),
  subclonal_spread=list(vaf_sub=function(n) sample(c(0.125,0.25,0.375),n,replace=T,prob=c(0.25,0.6,0.15))),
  overdispersed=list(vaf_sub=function(n) rep(0.27,n),rho=0.02),
  low_vaf_artefacts=list(vaf_sub=function(n) rep(0.27,n),n_artefact=10)
)

res=list()
for(sample_ID in c("MD7634","MD7635")){
  load(paste0("annotated_muts/annotated_mut_set_",sample_ID,"_postMS_reduced_a_j_vaf_post_mix_post_dup"))
  tree=read.tree(paste0("tree_files/tree_",sample_ID,"_postMS_reduced_a_j_vaf_post_mix_post_dup.tree"))
  CM=filtered_muts$COMB_mats.tree.build
  samples=setdiff(tree$tip.label,"Ancestral")
  depths=lapply(setNames(samples,samples),function(s) sample_depths(tree,CM,s))
  for(sc in names(scenarios)) for(rep in 1:n_reps){
    set.seed(1000*rep+match(sc,names(scenarios))+ifelse(sample_ID=="MD7634",0,500))
    truth=data.frame(sample=samples,N_c=round(runif(length(samples),35,80)))
    truth$N_s=round(truth$N_c*runif(length(samples),0.3,1.6))
    a=scenarios[[sc]]
    priv=simulate_mouse(depths,truth,a$vaf_sub,rho=ifelse(is.null(a$rho),0,a$rho),pc_sd=ifelse(is.null(a$pc_sd),0,a$pc_sd),
                        n_artefact=ifelse(is.null(a$n_artefact),0,a$n_artefact))
    #Beta-binomial (default) and binomial versions of the new method
    t0=Sys.time(); fit=fit_mixture(priv); el=as.numeric(difftime(Sys.time(),t0,units="secs"))
    fit_b=fit_mixture(priv,overdispersion=FALSE)
    sens=sapply(depths,clonal_sensitivity,fit=fit); sens_b=sapply(depths,clonal_sensitivity,fit=fit_b)
    est=correct_private(fit,sens); est_b=correct_private(fit_b,sens_b)%>%dplyr::select(sample,N_binomial=N_clonal)
    r=truth%>%left_join(est%>%dplyr::select(sample,n_obs,N_clonal,sensitivity),by="sample")%>%left_join(est_b,by="sample")%>%
      mutate(n_obs=coalesce(n_obs,0L),N_clonal=coalesce(N_clonal,0),N_binomial=coalesce(N_binomial,0),
             rmd=as.numeric(rmd_method(priv)[sample]),rmd=coalesce(rmd,0),
             sens_only=n_obs/sensitivity,median_depth=sapply(depths[sample],function(d) median(d$auto)),
             mouse=sample_ID,scenario=sc,rep=rep)
    cat(sample_ID,sc,rep,": fit",round(el),"s | pc",round(fit$pc,3),"ps",round(fit$ps,3),"rho",signif(fit$rho,2),"delta",round(fit$delta,3),"\n")
    res[[length(res)+1]]=r
  }
}
res=bind_rows(res)
write.csv(res,file.path(out_dir,"simulation_per_sample.csv"),row.names=F)

#Summary: relative bias of the mean, per-sample relative RMSE, correlation with truth, Spearman of error with depth
summ=res%>%tidyr::pivot_longer(c(N_clonal,N_binomial,rmd,n_obs,sens_only),names_to="method",values_to="estimate")%>%
  group_by(scenario,method)%>%
  summarise(bias_pct=round(100*(mean(estimate)/mean(N_c)-1),1),
            rmse_pct=round(100*sqrt(mean((estimate-N_c)^2))/mean(N_c),1),
            cor_truth=round(cor(estimate,N_c),3),
            error_vs_depth=round(cor(estimate-N_c,median_depth,method="s"),2),.groups="drop")%>%
  mutate(method=factor(method,levels=c("N_clonal","N_binomial","rmd","sens_only","n_obs"),
                       labels=c("new, beta-binomial","new, binomial","Rmd hard calls","sensitivity only","observed")))%>%arrange(scenario,method)
write.csv(summ,file.path(out_dir,"simulation_summary.csv"),row.names=F)
print(as.data.frame(summ))
