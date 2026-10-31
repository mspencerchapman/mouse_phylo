#Detection-aware correction of phylogeny branch lengths for sensitivity and sub-clonal private mutations
#
#Private mutations are only on the tree if they pass the tree-building 'vaf' filter (NV >= 3 & NV/NR >= 0.2) and,
#in addition, calls with few variant reads are only partly detected by the caller/filters (delta). The private
#mutations of each colony are a mixture of clonal (VAF ~0.5) and sub-clonal (VAF ~0.27) mutations; the mixture is
#fitted conditional on detection, so the sub-clonal fraction is not confounded with depth.
#Corrected private branch = sum of per-mutation clonal probabilities / clonal sensitivity of the sample.
#Shared branches are divided by the probability that at least one descendant sample detects a clonal mutation.
#Male non-PAR X/Y mutations are hemizygous (VAF = 2 x diploid-equivalent VAF, about half depth).
suppressMessages({library(dplyr); library(ape)})

#Detection rule of the tree-building filter
detect_rule=function(NV,NR,min_nv=3,min_vaf=0.2) NV>=min_nv & NV>=min_vaf*NR

#Detection probability given NV: delta[j] for NV = min_nv+j-1, 1 above
detect_weight=function(NV,delta,min_nv=3){
  d=rep(1,length(NV))
  for(j in seq_along(delta)) d[NV==min_nv+j-1]=delta[j]
  d
}

#Beta-binomial density, vectorised (rho = overdispersion; rho = 0 is binomial)
dbb=function(x,n,p,rho=0){
  if(rho<1e-8) return(dbinom(x,n,p))
  a=p*(1-rho)/rho; b=(1-p)*(1-rho)/rho
  exp(lchoose(n,x)+lbeta(x+a,n-x+b)-lbeta(a,b))
}

#P(mutation passes the rule and is detected | depth NR, VAF p) - vectorised over NR and p
pass_prob=function(NR,p,delta=numeric(0),min_nv=3,min_vaf=0.2,rho=0){
  p=rep_len(p,length(NR)); k=pmax(min_nv,ceiling(min_vaf*NR-1e-9))
  if(rho<1e-8) {
    tot=pbinom(k-1,NR,p,lower.tail=F)
    for(j in seq_along(delta)){v=min_nv+j-1; i=k<=v&NR>=v; tot[i]=tot[i]-(1-delta[j])*dbinom(v,NR[i],p[i])}
    return(ifelse(NR<k,0,tot))
  }
  #Beta-binomial: sum over NV for each distinct (NR, p)
  key=paste(NR,p); u=!duplicated(key)
  val=mapply(function(n,k,p) if(n<k) 0 else {v=k:n; sum(dbb(v,n,p,rho)*detect_weight(v,delta,min_nv))},NR[u],k[u],p[u])
  setNames(val,key[u])[key]
}

#Copy factor for male mouse non-PAR X/Y (GRCm38 PAR: chrX 169,969,759-170,931,299; chrY 90,745,845-91,644,698)
copy_factor_mouse=function(chr,pos,sex="male"){
  if(sex!="male") return(rep(1,length(chr)))
  par=(chr=="X"&pos>=169969759&pos<=170931299)|(chr=="Y"&pos>=90745845&pos<=91644698)
  ifelse(chr%in%c("X","Y")&!par,2,1)
}

#Ancestral nodes (excluding the tip itself) of a tip
anc_nodes=function(tree,tip){
  out=c(); nd=tip; root=length(tree$tip.label)+1
  while(nd!=root){nd=tree$edge[tree$edge[,2]==nd,1]; out=c(out,nd)}
  out
}

#Private mutations of each sample (COMB_mats = filtered_muts$COMB_mats.tree.build), with copy factor
get_private=function(tree,COMB_mats,sex="male",outgroup="Ancestral",mut_types="SNV"){
  mat=COMB_mats$mat
  lapply(setdiff(tree$tip.label,outgroup),function(s){
    i=which(mat$node==which(tree$tip.label==s)&mat$Mut_type%in%mut_types); m=mat$mut_ref[i]
    data.frame(sample=s,mut_ref=m,NV=COMB_mats$NV[m,s],NR=COMB_mats$NR[m,s],cf=copy_factor_mouse(mat$Chrom[i],mat$Pos[i],sex))
  })%>%bind_rows()
}

#Depth of a sample at the mutation sites of other branches (not selected on this sample), split autosome / X-Y
sample_depths=function(tree,COMB_mats,s,sex="male"){
  mat=COMB_mats$mat; cf=copy_factor_mouse(mat$Chrom,mat$Pos,sex); keep=mat$node!=which(tree$tip.label==s)
  list(auto=COMB_mats$NR[mat$mut_ref[keep&cf==1],s],xy=COMB_mats$NR[mat$mut_ref[keep&cf==2],s])
}

#Fit the detection-conditioned mixture for all samples of one individual
#priv: data.frame(sample, NV, NR, cf) of detected private mutations
#Global parameters: clonal VAF pc, sub-clonal VAF ps (diploid-equivalent), beta-binomial overdispersion rho_c / rho_s
#(if overdispersion=TRUE), detection delta for NV = min_nv, min_nv+1
#Per-sample parameter: true-space clonal weight w (proportion of clonal among all private mutations, detected or not),
#profiled out for each set of global parameters
#Colonies are single-cell derived, so pc is shared by all samples (impure colonies should be removed beforehand)
fit_mixture=function(priv,n_delta=2,overdispersion=TRUE,min_nv=3,min_vaf=0.2,start=c(pc=0.5,ps=0.27)){
  priv=priv[detect_rule(priv$NV,priv$NR,min_nv,min_vaf),]
  pcm=function(p) pmin(priv$cf*p,0.999)
  setup=function(t) list(
    fc=dbb(priv$NV,priv$NR,pcm(t$pc),t$rho[1]),fs=dbb(priv$NV,priv$NR,pcm(t$ps),t$rho[2]),
    Sc=pass_prob(priv$NR,pcm(t$pc),t$delta,min_nv,min_vaf,t$rho[1]),Ss=pass_prob(priv$NR,pcm(t$ps),t$delta,min_nv,min_vaf,t$rho[2]),
    ld=log(detect_weight(priv$NV,t$delta,min_nv)))
  nll_w=function(z,i) function(w) -sum(log(w*z$fc[i]+(1-w)*z$fs[i])+z$ld[i]-log(w*z$Sc[i]+(1-w)*z$Ss[i]))
  by_s=split(seq_len(nrow(priv)),priv$sample)
  n_rho=ifelse(overdispersion,2,0)
  tr=function(x) list(pc=plogis(x[1]),ps=plogis(x[2]),rho=if(n_rho>0) exp(x[3:4]) else c(0,0),delta=plogis(x[(3+n_rho):length(x)]))
  nll=function(x){z=setup(tr(x)); sum(sapply(by_s,function(i) optimize(nll_w(z,i),c(1e-6,1-1e-6))$objective))}
  o=optim(c(qlogis(start),rep(log(0.005),n_rho),rep(2,n_delta)),nll,control=list(reltol=1e-8,maxit=3000))
  t=tr(o$par); z=setup(t)
  w=sapply(by_s,function(i) optimize(nll_w(z,i),c(1e-6,1-1e-6))$minimum)
  wi=w[priv$sample]
  priv$p_clonal=wi*z$fc/(wi*z$fc+(1-wi)*z$fs)
  list(pc=t$pc,ps=t$ps,rho=t$rho,delta=t$delta,w=w,priv=priv,nll=o$value,convergence=o$convergence,min_nv=min_nv,min_vaf=min_vaf)
}

#Caller sensitivity (germline-SNP pick-up rate of CaVEMan + LCM filters) recovered from a tree built with the pipeline's
#sensitivity adjustment: private (tip) edge = SNV count / SNV sensitivity
caller_sensitivity_from_tree=function(tree,COMB_mats,outgroup="Ancestral",mut_types="SNV"){
  mat=COMB_mats$mat; samples=setdiff(tree$tip.label,outgroup)
  sapply(setNames(samples,samples),function(s){tip=which(tree$tip.label==s)
    sum(mat$node==tip&mat$Mut_type%in%mut_types)/tree$edge.length[tree$edge[,2]==tip]})
}

#Clonal sensitivity of a sample: probability that a clonal mutation is called and passes the tree-building rule
#xy_share = share of clonal mutations on male X/Y (~171 Mb haploid X vs 2 x 2.5 Gb autosomes -> ~0.033)
#caller_sens = NULL: caller losses are only the fitted NV=3/4 detection (delta), from the depth model alone
#caller_sens given (germline-SNP pick-up rate): caller_sens x P(pass rule | clonal); delta is set to 1 here because
#losses of low-NV calls are part of the caller sensitivity
clonal_sensitivity=function(depths,fit,caller_sens=NULL,xy_share=0.033){
  delta=if(is.null(caller_sens)) fit$delta else rep(1,length(fit$delta))
  pp=function(d,p) mean(pass_prob(d,p,delta,fit$min_nv,fit$min_vaf,fit$rho[1]))
  rule=pp(depths$auto,fit$pc)
  if(length(depths$xy)>0) rule=(1-xy_share)*rule+xy_share*pp(depths$xy,min(2*fit$pc,0.999))
  unname(if(is.null(caller_sens)) rule else caller_sens*rule)
}

#Corrected private burden per sample
correct_private=function(fit,sens){
  fit$priv%>%group_by(sample)%>%summarise(n_obs=n(),n_obs_clonal=sum(p_clonal),.groups="drop")%>%
    mutate(w_clonal=fit$w[sample],sensitivity=sens[sample],N_clonal=n_obs_clonal/sensitivity)
}

#Branch lengths from raw mutation counts; shared branches divided by the probability that at least one descendant
#sample detects a clonal mutation; private branches replaced by the corrected clonal burden
corrected_tree=function(tree,COMB_mats,sens,private_df,outgroup="Ancestral",mut_types="SNV"){
  mat=COMB_mats$mat; nt=length(tree$tip.label); tc=tree
  counts=table(factor(mat$node[mat$Mut_type%in%mut_types],levels=seq_len(max(tree$edge))))
  tc$edge.length=as.numeric(counts[tree$edge[,2]])
  for(i in which(tree$edge[,2]>nt)){
    tips=intersect(tree$tip.label[phangorn::Descendants(tree,tree$edge[i,2],"tips")[[1]]],names(sens))
    tc$edge.length[i]=tc$edge.length[i]/(1-prod(1-sens[tips]))
  }
  for(s in private_df$sample) tc$edge.length[tc$edge[,2]==which(tree$tip.label==s)]=private_df$N_clonal[private_df$sample==s]
  if(outgroup%in%tree$tip.label) tc$edge.length[tc$edge[,2]==which(tree$tip.label==outgroup)]=0
  tc
}

#Whole procedure for one individual. Returns per-sample table, fitted model and the corrected tree
#caller_sensitivity: "from_tree" (pipeline germline-SNP sensitivity recovered from the input tree), a named vector,
#or NULL (depth model only)
correct_branch_lengths=function(tree,COMB_mats,sex="male",outgroup="Ancestral",mut_types="SNV",caller_sensitivity="from_tree",
                               n_delta=2,overdispersion=TRUE,xy_share=0.033){
  priv=get_private(tree,COMB_mats,sex,outgroup,mut_types)
  fit=fit_mixture(priv,n_delta,overdispersion)
  samples=setdiff(tree$tip.label,outgroup)
  if(identical(caller_sensitivity,"from_tree")) caller_sensitivity=caller_sensitivity_from_tree(tree,COMB_mats,outgroup,mut_types)
  depths=lapply(setNames(samples,samples),function(s) sample_depths(tree,COMB_mats,s,sex))
  sens=sapply(samples,function(s) clonal_sensitivity(depths[[s]],fit,if(is.null(caller_sensitivity)) NULL else caller_sensitivity[s],xy_share))
  res=correct_private(fit,sens)
  res=data.frame(sample=samples)%>%left_join(res,by="sample")%>%
    mutate(n_obs=coalesce(n_obs,0L),n_obs_clonal=coalesce(n_obs_clonal,0),N_clonal=coalesce(N_clonal,0),
           caller_sensitivity=if(is.null(caller_sensitivity)) NA else caller_sensitivity[sample],
           sensitivity=sens[sample],median_depth=sapply(depths[sample],function(d) median(d$auto)))
  tc=corrected_tree(tree,COMB_mats,sens,res,outgroup,mut_types)
  res$shared_burden=sapply(res$sample,function(s) sum(tc$edge.length[tc$edge[,2]%in%anc_nodes(tc,which(tc$tip.label==s))]))
  res$total_burden=res$shared_burden+res$N_clonal
  list(samples=res,fit=fit[c("pc","ps","rho","delta","w","nll","convergence","min_nv","min_vaf")],priv=fit$priv,tree=tc)
}
