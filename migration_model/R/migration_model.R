#Localisation + migration model for the anatomical (bone) distribution of HSC colonies on an ultrametric phylogeny
#
#Model (time = molecular time on the ultrametric tree, root = 0, tips = D):
# 1. Before the localisation time T lineages have no fixed location. At T, every lineage crossing T is placed in a
#    bone k with probability pi_k, independently of the other lineages.
# 2. After T, each lineage has migration events at rate mu per unit molecular time; at an event the destination bone
#    is drawn from pi (F81 model, so it can be the same bone). Moves to a different bone happen at rate
#    mu * (1 - pi_current); averaged over bones the effective move rate is mu * (1 - sum(pi^2)).
#Because pi is the stationary distribution, the likelihood of the tip tissues is the product over the lineages
#crossing T of sum_k pi_k L_k(subtree), with L computed by Felsenstein pruning below T.
suppressMessages({library(ape);library(dplyr)})

#Corrected ultrametric tree as in Mouse_phylo_clustering.Rmd (make.ultrametric.tree, rescaled to the mean burden)
make_ultra=function(tree,outgroup="Ancestral"){
  u=make.ultrametric.tree(tree); u$coords=NULL
  u$edge.length[is.infinite(u$edge.length)|is.na(u$edge.length)]=0
  u$edge.length=u$edge.length*mean(get_mut_burden(tree))
  if(outgroup%in%u$tip.label) u=drop.tip(u,outgroup)
  u
}

#Precompute what the likelihood needs: postorder edges, node heights, tip states
prep_tree=function(tree,tip_states,levels=NULL){
  tree=reorder(tree,"postorder")
  if(is.null(levels)) levels=sort(unique(tip_states))
  states=factor(tip_states[tree$tip.label],levels=levels)
  if(any(is.na(states))) stop("tips without a state")
  h=node.depth.edgelength(tree)
  mrca_h=matrix(h[mrca(tree)[tree$tip.label,tree$tip.label]],Ntip(tree))  #height of the MRCA of each tip pair
  list(tree=tree,mrca_h=mrca_h,n_tip=Ntip(tree),n_node=Ntip(tree)+tree$Nnode,edge=tree$edge,len=tree$edge.length,
       h=h,D=max(h[1:Ntip(tree)]),states=states,levels=levels,K=length(levels),
       pi_obs=as.numeric(table(states))/length(states))
}

#Felsenstein pruning: partial likelihoods of every node for one migration rate mu.
#
#For each node n and bone k, L[k,n] = P(observed bones of all tips below n | lineage n is in bone k at node n).
#Tips: L = 1 for the observed bone, 0 for the others. Edges are visited in postorder (children before parents), and
#each child's vector is passed up its branch and multiplied into the parent's:
#  L[,parent] = L[,parent] * (P(t) %*% L[,child]),  t = branch length
#P(t) is the F81 transition matrix: with probability exp(-mu t) there is no migration event on the branch (stay in the
#same bone); otherwise the bone after the last event is drawn from pi. So P(t)[k,j] = e*(k==j) + (1-e)*pi_j with
#e = exp(-mu t), and P(t) %*% L = e*L + (1-e)*sum(pi*L), which needs no matrix.
#Partials are computed for the whole tree, ignoring T; loglik_grid() then uses them only below T.
#Numerical underflow: each child vector is rescaled to max 1 once it is complete, and the log scale factors are
#carried up in lsc, so the true partial of node n is L[,n] * exp(lsc[n]).
#Returns list(L = K x n_node matrix, lsc = log scale per node).
partials=function(pt,mu,pi){
  K=pt$K; L=matrix(1,K,pt$n_node); lsc=numeric(pt$n_node)
  L[,1:pt$n_tip]=0; L[cbind(as.integer(pt$states),1:pt$n_tip)]=1
  for(i in seq_len(nrow(pt$edge))){
    p=pt$edge[i,1]; c=pt$edge[i,2]; e=exp(-mu*pt$len[i])
    if(c>pt$n_tip){s=max(L[,c]); if(s==0) s=1; L[,c]=L[,c]/s; lsc[c]=lsc[c]+log(s)} #child complete (postorder): rescale
    v=e*L[,c]+(1-e)*sum(pi*L[,c])            #P(t) %*% L[,c] for the F81 model
    L[,p]=L[,p]*v; lsc[p]=lsc[p]+lsc[c]
  }
  list(L=L,lsc=lsc)
}

#Log-likelihood surface over a grid of T and mu (mu = 0 allowed).
#
#Lineages crossing T are the edges with parent above T (height < T) and child at or below it (height >= T). At T each
#such lineage starts in bone k with probability pi_k, independently of the others, so
#  logL(T, mu) = sum over crossing edges of log( sum_k pi_k * [P(t') %*% L[,child]]_k ),  t' = height(child) - T.
#Because pi is the stationary distribution of the F81 process, pi %*% P(t') = pi for any t', so this equals
#log(sum_k pi_k L[k,child]): the part of the crossing branch between T and the child node drops out. Hence the
#partials (one pruning pass per mu) can be reused for every T, and only the set of crossing edges changes with T.
#Everything above T (before localisation) carries no information about bones and contributes nothing.
#Checked against brute-force enumeration and forward simulation in migration_model/tests/test_likelihood.R.
loglik_grid=function(pt,T_grid,mu_grid,pi=pt$pi_obs){
  hp=pt$h[pt$edge[,1]]; hc=pt$h[pt$edge[,2]]
  out=matrix(NA,length(T_grid),length(mu_grid),dimnames=list(T=T_grid,mu=mu_grid))
  for(j in seq_along(mu_grid)){
    P=partials(pt,mu_grid[j],pi)
    node_ll=log(colSums(pi*P$L))+P$lsc        #log(sum_k pi_k L[k,n]) for every node n
    for(i in seq_along(T_grid)){
      cr=pt$edge[hp<T_grid[i]&hc>=T_grid[i],2] #child nodes of the edges crossing T
      out[i,j]=sum(node_ll[cr])
    }
  }
  out
}

#Forward simulation of tip states under the model
simulate_states=function(pt,T,mu,pi=pt$pi_obs){
  st=rep(NA_integer_,pt$n_node); hp=pt$h[pt$edge[,1]]; hc=pt$h[pt$edge[,2]]
  for(i in rev(seq_len(nrow(pt$edge)))){     #preorder = reversed postorder
    p=pt$edge[i,1]; c=pt$edge[i,2]
    if(hc[i]<T) next                          #not yet localised
    if(hp[i]<T){st[c]=sample.int(pt$K,1,prob=pi); t=hc[i]-T} else {st[c]=st[p]; t=pt$len[i]}
    if(runif(1)>exp(-mu*t)) st[c]=sample.int(pt$K,1,prob=pi)
  }
  setNames(pt$levels[st[1:pt$n_tip]],pt$tree$tip.label)
}

#Posterior on the grid given priors (vectors over T_grid and mu_grid, need not be normalised)
grid_posterior=function(ll,prior_T,prior_mu){
  lp=ll+outer(log(prior_T),log(prior_mu),"+")
  p=exp(lp-max(lp,na.rm=T)); p[is.na(p)]=0; p/sum(p)
}

#Summaries of a posterior grid
#Calendar calibration from Kapadia et al. 2025 (Nature 641:681; C57BL/6J HSC/MPP colonies): ~48.2 mutations at birth
#(mixed-model intercept), 45.3 per year after birth, ~10 at gastrulation (~E6.5-7), so ~2.9 per day before birth
#(E19.5). 'scale' converts molecular time on our trees to Kapadia-equivalent units (1 = take our burdens at face value).
mol_to_age=function(m,scale=1,birth=48.2,rate_year=45.3,rate_day_prenatal=2.9){
  k=m/scale; ifelse(k>=birth,(k-birth)/rate_year,-(birth-k)/rate_day_prenatal/365)}  #age in years (negative = before birth)
age_to_mol=function(age,scale=1,birth=48.2,rate_year=45.3,rate_day_prenatal=2.9){
  scale*ifelse(age>=0,birth+age*rate_year,birth+age*365*rate_day_prenatal)}

post_summary=function(post,T_grid,mu_grid,D,pi,probs=c(0.025,0.1,0.5,0.9,0.975),years_after=NULL){
  q=function(x,w){o=order(x);cw=cumsum(w[o])/sum(w);sapply(probs,function(p) x[o][which(cw>=p)[1]])}
  meff=mu_grid*(1-sum(pi^2))
  Tm=matrix(T_grid,length(T_grid),length(mu_grid)); Mm=matrix(meff,length(T_grid),length(mu_grid),byrow=T)
  n_mig=Mm*pmax(D-Tm,0)                         #expected moves per lineage from T to sampling
  p_stay=exp(-n_mig)                            #probability a sampled HSC never moved since T
  out=rbind(T=q(Tm,post),move_rate=q(Mm,post),moves_per_lineage=q(n_mig,post),p_never_moved=q(p_stay,post))
  #years_after: function of T giving the calendar time from localisation to sampling; adds moves per year
  if(!is.null(years_after)) out=rbind(out,years_since_T=q(years_after(Tm),post),moves_per_year=q(n_mig/years_after(Tm),post))
  out%>%
    `colnames<-`(paste0("q",probs*100))
}

#Clustering statistic: for tip pairs whose MRCA is at height >= h, the fraction from the same bone
pair_stat=function(pt,states,h_cut){
  st=as.character(states[pt$tree$tip.label]); same=outer(st,st,"==")
  ut=upper.tri(same); mh=pt$mrca_h[ut]; sm=same[ut]
  sapply(h_cut,function(x){i=mh>=x; if(!any(i)) NA else mean(sm[i])})
}
