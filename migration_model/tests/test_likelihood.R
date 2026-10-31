#Checks of the exact likelihood in migration_model/R/migration_model.R (loglik_grid / partials).
#Run from the project directory: Rscript migration_model/tests/test_likelihood.R  (~10 s). Stops with an error if a
#check fails. Also sourced by Migration_model.Rmd.
#
#1. Brute force: on a 5-tip toy tree, sum the probability of every possible history explicitly (bone of each lineage
#   at T and of every internal node below T), without pruning, and compare with loglik_grid().
#2. Normalisation: the probabilities of all possible tip patterns sum to 1.
#3. Forward simulation: the frequency of tip patterns in simulations from simulate_states() matches the exact
#   probabilities.
#4. Large tree (checks the numerical rescaling): with very fast migration, or localisation just before sampling, the
#   bones are independent (rate >> 1/shortest branch), so logL = sum over tips of log(pi of the tip's bone).
suppressMessages(library(ape))
source("migration_model/R/migration_model.R")
set.seed(20261031)
pass=function(ok,msg) {if(!isTRUE(ok)) stop("FAILED: ",msg); cat("ok:",msg,"\n")}

#Toy tree: root at 0; internal nodes at heights 2 (a,b), 0.5 (c + (d,e)) and 2.5 (d,e); tips at 3
tr=read.tree(text="((a:1,b:1):2,(c:2.5,(d:0.5,e:0.5):2):0.5);")
bones=c("F","I","T"); pi=c(0.2,0.3,0.5)
tips=c("a","b","c","d","e")

#F81 transition matrix
P=function(t,mu) exp(-mu*t)*diag(length(pi))+(1-exp(-mu*t))*matrix(pi,length(pi),length(pi),byrow=TRUE)

#Brute-force likelihood: enumerate the bone at T of every crossing lineage and of every internal node below T
brute=function(tree,states,T,mu){
  h=node.depth.edgelength(tree); nt=Ntip(tree); E=tree$edge; hp=h[E[,1]]; hc=h[E[,2]]
  cross=which(hp<T&hc>=T); inner=which(hp>=T)                        #crossing edges; edges entirely below T
  free=unique(c(E[cross,2],E[inner,2]))
  free=free[free>nt]                                                 #internal nodes whose bone is unknown
  grid=expand.grid(rep(list(seq_along(pi)),length(cross)+length(free)))
  tip_state=match(states[tree$tip.label],bones)
  total=0
  for(r in seq_len(nrow(grid))){
    s_T=unlist(grid[r,seq_along(cross)]); s_node=rep(NA,max(E)); s_node[seq_len(nt)]=tip_state
    if(length(free)>0) s_node[free]=unlist(grid[r,length(cross)+seq_along(free)])
    p=prod(pi[s_T])*prod(sapply(seq_along(cross),function(k){i=cross[k]; P(hc[i]-T,mu)[s_T[k],s_node[E[i,2]]]}))
    if(length(inner)>0) p=p*prod(sapply(inner,function(i) P(tree$edge.length[i],mu)[s_node[E[i,1]],s_node[E[i,2]]]))
    total=total+p
  }
  total
}

#1. Brute force vs pruning, several localisation times and migration rates (including mu = 0)
for(st in list(c(a="F",b="F",c="T",d="I",e="F"),c(a="I",b="I",c="T",d="T",e="T"))){
  pt=prep_tree(tr,st[tips],levels=bones)
  for(T in c(0.2,1,2.2,2.9)) for(mu in c(0,0.4,3)){
    b=brute(pt$tree,st,T,mu); g=exp(loglik_grid(pt,T,mu,pi))
    pass(isTRUE(all.equal(b,unname(g[1,1]),tolerance=1e-10)),sprintf("brute force = pruning (%s, T=%.1f, mu=%.1f): %.6g",paste(st,collapse=""),T,mu,b))
  }
}

#2. Probabilities of all 3^5 tip patterns sum to 1
patterns=expand.grid(rep(list(bones),length(tips)),stringsAsFactors=FALSE); names(patterns)=tips
for(T in c(0.2,2.2)) for(mu in c(0,0.4)){
  tot=sum(apply(patterns,1,function(x) exp(loglik_grid(prep_tree(tr,setNames(x,tips),levels=bones),T,mu,pi))))
  pass(abs(tot-1)<1e-10,sprintf("probabilities sum to 1 (T=%.1f, mu=%.1f): %.12f",T,mu,tot))
}

#3. Forward simulation vs exact probabilities of the most common patterns
T=0.2; mu=0.4; n=1e5
pt=prep_tree(tr,c(a="F",b="F",c="T",d="I",e="F"),levels=bones)
sims=replicate(n,paste(simulate_states(pt,T,mu,pi)[tips],collapse=""))
tab=sort(table(sims),decreasing=TRUE)[1:8]
for(pat in names(tab)){
  e=exp(loglik_grid(prep_tree(tr,setNames(strsplit(pat,"")[[1]],tips),levels=bones),T,mu,pi))[1,1]
  z=(tab[[pat]]-n*e)/sqrt(n*e*(1-e))
  pass(abs(z)<4,sprintf("simulated frequency %s: %.4f vs exact %.4f (z = %.1f)",pat,tab[[pat]]/n,e,z))
}

#4. Large tree: limits where bones are independent of the tree (checks the rescaling against underflow)
big=rcoal(300); st=setNames(sample(bones,300,replace=TRUE,prob=pi),big$tip.label)
pt=prep_tree(big,st,levels=bones); indep=sum(log(pi[match(st,bones)]))
pass(abs(loglik_grid(pt,0.01,1e3/min(big$edge.length),pi)[1,1]-indep)<1e-6,"300 tips, very fast migration: logL = independent bones")
pass(abs(loglik_grid(pt,pt$D-1e-9,0,pi)[1,1]-indep)<1e-6,"300 tips, localisation just before sampling: logL = independent bones")
cat("All likelihood checks passed\n")
