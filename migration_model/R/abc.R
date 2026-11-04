#Simulation-based (grid ABC / synthetic likelihood) inference for the localisation + migration model, using clustering summary statistics
#(clade x bone chi-squared statistics at molecular-time cut-offs, as in Mouse_phylo_clustering.Rmd, optionally with
#same-bone pair fractions) instead of the full likelihood. Needs migration_model.R.

#Vectorised forward simulation: n_rep replicates of the tip states for one (T, mu); returns n_tip x n_rep integer matrix
simulate_states_multi=function(pt,T,mu,n_rep,pi=pt$pi_obs){
  st=matrix(NA_integer_,pt$n_node,n_rep); hp=pt$h[pt$edge[,1]]; hc=pt$h[pt$edge[,2]]
  draw=function(n) sample.int(pt$K,n,replace=TRUE,prob=pi)
  for(i in rev(seq_len(nrow(pt$edge)))){
    p=pt$edge[i,1]; c=pt$edge[i,2]
    if(hc[i]<T) next
    if(hp[i]<T){st[c,]=draw(n_rep); t=hc[i]-T} else {st[c,]=st[p,]; t=pt$len[i]}
    jump=runif(n_rep)>exp(-mu*t)
    if(any(jump)) st[c,jump]=draw(sum(jump))
  }
  st[1:pt$n_tip,,drop=FALSE]
}

#Clade membership of each tip for each cut-off: clades = all branches crossing the cut-off, including single colonies
#(the same clades as get_expanded_clade_nodes(..., min_clonal_fraction=0, min_samples=1) in Mouse_phylo_clustering.Rmd).
#Every tip is then in exactly one clade, so the clade x bone table partitions all colonies. (Before 2026-11-03 single
#colonies were left out; see archive/abc_R_before_singletons.R.)
clade_index=function(pt,cuts){
  tr=pt$tree; hp=pt$h[tr$edge[,1]]; hc=pt$h[tr$edge[,2]]
  lapply(cuts,function(h){
    nodes=tr$edge[hp<h&hc>=h,2]
    id=rep(NA_integer_,pt$n_tip)
    for(k in seq_along(nodes)) id[match(getTips(tr,nodes[k]),tr$tip.label)]=k
    id})
}

#Chi-squared statistic of the clade x bone table for each cut-off; states = integer vector of tip states
chisq_stats=function(states,cl_idx,K){
  sapply(cl_idx,function(id){
    ok=!is.na(id); n_cl=max(c(0,id[ok])); if(n_cl<2) return(NA_real_)
    M=matrix(tabulate((id[ok]-1)*K+states[ok],n_cl*K),n_cl,K,byrow=TRUE); M=M[,colSums(M)>0,drop=FALSE]
    if(ncol(M)<2) return(0)
    E=outer(rowSums(M),colSums(M))/sum(M); sum((M-E)^2/E)})
}

#Same-bone fraction of tip pairs in MRCA-height bands
pair_index=function(pt,bands){ut=upper.tri(pt$mrca_h)
  lapply(bands,function(b) which(ut&pt$mrca_h>=b[1]&pt$mrca_h<b[2],arr.ind=TRUE))}
pair_stats=function(states,pr_idx) sapply(pr_idx,function(ix) if(nrow(ix)==0) NA_real_ else mean(states[ix[,1]]==states[ix[,2]]))

#Summary statistics of one set of tip states
summaries=function(states,ss) c(chisq_stats(states,ss$cl_idx,ss$K),if(!is.null(ss$pr_idx)) pair_stats(states,ss$pr_idx))
make_ss=function(pt,cuts,bands=NULL) list(cl_idx=clade_index(pt,cuts),pr_idx=if(is.null(bands)) NULL else pair_index(pt,bands),K=pt$K,
  names=c(paste0("chisq_",cuts),if(!is.null(bands)) sapply(bands,function(b) paste0("same_",b[1],"_",b[2]))))

#Reference table: n_sim simulations for every (T, mu) grid cell; returns list(grid, S = array cell x sim x stat)
abc_reference=function(pt,ss,T_grid,mu_grid,n_sim,n_cores=1){
  grid=expand.grid(T=T_grid,mu=mu_grid)
  res=parallel::mclapply(seq_len(nrow(grid)),function(g){
    st=simulate_states_multi(pt,grid$T[g],grid$mu[g],n_sim); t(apply(st,2,summaries,ss=ss))},mc.cores=n_cores)
  S=array(unlist(lapply(res,as.vector)),dim=c(n_sim,length(ss$names),nrow(grid)))
  list(grid=grid,S=aperm(S,c(3,1,2)),names=ss$names)
}

#Kernel ABC likelihood of observed summaries for every grid cell: Gaussian kernel on the Mahalanobis distance, whitened
#with the covariance of the statistics in one reference cell (default: the last cell = latest T, highest mu, i.e. ~no
#clustering, where the discrimination matters). Scaling by the SD across the whole table, or by a covariance pooled over
#all cells, is dominated by the noisy strongly-clustered cells and makes every cell look close to the data.
#Bandwidth h = quantile 'q' of all prior-weighted distances (~ fraction of simulations "accepted").
abc_likelihood=function(ref,s_obs,prior_cell=NULL,q=0.01,use=NULL,whiten_cell=dim(ref$S)[1]){
  use=if(is.null(use)) which(!is.na(s_obs)) else intersect(use,which(!is.na(s_obs)))
  S=ref$S[,,use,drop=FALSE]; S[is.na(S)]=0
  n_cell=dim(S)[1]; n_sim=dim(S)[2]
  e=eigen(cov(matrix(S[whiten_cell,,],n_sim)),symmetric=TRUE); W=e$vectors%*%diag(1/sqrt(pmax(e$values,max(e$values)*1e-6)),length(use))
  X=matrix(S,n_cell*n_sim,length(use))                       #rows: cell fastest, then simulation
  d=matrix(sqrt(rowSums((sweep(X,2,s_obs[use])%*%W)^2)),n_cell,n_sim)
  if(is.null(prior_cell)) prior_cell=rep(1,n_cell)
  h=Hmisc_wtd_quantile(as.vector(d),rep(prior_cell/n_sim,n_sim),q)
  list(lik=rowMeans(exp(-0.5*(d/h)^2)),h=h,d=d)
}

#Synthetic likelihood (Wood 2010): in each cell the summaries are approximated by a multivariate normal with the
#cell's own simulated mean and covariance. Chi-squared statistics are square-root transformed (closer to normal).
#synthetic_prep() is done once per reference table; synthetic_loglik() is then cheap for any observed summaries.
sl_transform=function(X,names) {i=grepl("^chisq_",names); X[,i]=sqrt(pmax(X[,i],0)); X}
synthetic_prep=function(ref,use){
  lapply(seq_len(dim(ref$S)[1]),function(k){X=matrix(ref$S[k,,use],dim(ref$S)[2]); X[is.na(X)]=0
    X=sl_transform(X,ref$names[use]); e=eigen(cov(X),symmetric=TRUE)
    list(m=colMeans(X),V=e$vectors,v=pmax(e$values,1e-6))})  #floor: statistics that never vary in a cell
}
synthetic_loglik=function(prep,s_obs,use,names){x=s_obs[use]; x[is.na(x)]=0; x=as.vector(sl_transform(matrix(x,1),names[use]))
  sapply(prep,function(c){z=crossprod(c$V,x-c$m); -0.5*sum(z^2/c$v)-0.5*sum(log(c$v))})}
Hmisc_wtd_quantile=function(x,w,p){o=order(x); cw=cumsum(w[o])/sum(w); x[o][which(cw>=p)[1]]}
