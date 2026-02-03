library(dplyr)

source("/lustre/scratch126/casm/team154pc/ms56/my_programs/Prolonged_persistence_functions.R") #Source functions needed for the script

load("filtering_runs/annotated_muts/annotated_mut_set_MD7634_postMS_reduced_a_j_vaf_post_mix_post_dup")
tree=read.tree("filtering_runs/tree_files/tree_MD7634_postMS_reduced_a_j_vaf_post_mix_post_dup.tree")
details<-filtered_muts$COMB_mats.tree.build$mat

gfile="/lustre/scratch124/casm/team78pipelines/canpipe/live/ref/Mus_musculus/GRCm38/genome.fa"

plot_96profile=function(df,colnames=c("sampleID","chr","pos","ref","alt"),genomeFile) {
  require(Rsamtools)
  require(GenomicRanges)
  require(IRanges)
  mutations = data.frame(sampleID=df[[colnames[1]]],
                         chr=df[[colnames[2]]],
                         pos=df[[colnames[3]]],
                         ref=df[[colnames[4]]],
                         trinuc_ref= "-",
                         mut=df[[colnames[5]]])
  mutations_GRange<-GRanges(mutations$chr, IRanges::IRanges(mutations$pos-1, mutations$pos+1))
  mutations$trinuc_ref = as.vector(Rsamtools::scanFa(file=genomeFile,mutations_GRange ))
  
  ntcomp = c(T="A",G="C",C="G",A="T")
  mutations$sub = paste(mutations$ref,mutations$mut,sep=">")
  mutations$trinuc_ref_py = mutations$trinuc_ref
  for (j in 1:nrow(mutations)) {
    if (mutations$ref[j] %in% c("A","G")) { # Purine base
      mutations$sub[j] = paste(ntcomp[mutations$ref[j]],ntcomp[mutations$mut[j]],sep=">")
      mutations$trinuc_ref_py[j] = paste(ntcomp[rev(strsplit(mutations$trinuc_ref[j],split="")[[1]])],collapse="")
    }
  }
  freqs = table(paste(mutations$sub,paste(substr(mutations$trinuc_ref_py,1,1),substr(mutations$trinuc_ref_py,3,3),sep="-"),sep=","))
  sub_vec = c("C>A","C>G","C>T","T>A","T>C","T>G")
  ctx_vec = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
  full_vec = paste(rep(sub_vec,each=16),rep(ctx_vec,times=6),sep=",")
  freqs_full = freqs[full_vec]; freqs_full[is.na(freqs_full)] = 0; names(freqs_full) = full_vec
  
  xstr = paste(substr(full_vec,5,5), substr(full_vec,1,1), substr(full_vec,7,7), sep="")
  
  #pdf("Mut_Sig_fetal_shared.pdf")
  colvec = rep(c("dodgerblue","black","red","grey70","olivedrab3","plum2"),each=16)
  y = freqs_full; maxy = max(y)
  h = barplot(y, las=2, col=colvec, border=NA, ylim=c(0,maxy*1.5), space=1, cex.names=0.6, names.arg=xstr, ylab="Number mutations")
  for (j in 1:length(sub_vec)) {
    xpos = h[c((j-1)*16+1,j*16)]
    rect(xpos[1]-0.5, maxy*1.2, xpos[2]+0.5, maxy*1.3, border=NA, col=colvec[j*16])
    text(x=mean(xpos), y=maxy*1.3, pos=3, label=sub_vec[j])
  }
}

plot_96profile_mutref=function(mutref_vec,genomeFile) {
  mut_mat=stringr::str_split(mutref_vec,pattern="-",simplify=T)
  colnames(mut_mat)=c("chr","pos","ref","alt")
  df=as.data.frame(mut_mat)%>%
    mutate(pos=as.integer(pos))%>%
    mutate(sampleID="this_sample",.before=1)
  plot_96profile(df,genomeFile = genomeFile)
}

vaf_density_plot_final=function(sample,tree,COMB_mats,private_only=F){
  node <- which(tree$tip.label==sample)
  if(private_only) {relevant_nodes<-node} else {relevant_nodes<-get_ancestral_nodes(node,tree$edge)}
  sample_muts <- COMB_mats$mat$mut_ref[COMB_mats$mat$node %in% relevant_nodes]
  COMB_mats$NR[COMB_mats$NR == 0] <- 1
  dens <- density((COMB_mats$NV/COMB_mats$NR)[sample_muts,sample])
  plot(dens,xlim = c(0,1),main=sample)
  abline(v = dens$x[which.max(dens$y)])
  text(0.7, max(dens$y) - 0.2, paste("Peak VAF dens=",round(dens$x[which.max(dens$y)], digits = 2)),col="red",cex = 0.7)
}

binom_mix = function(x,size,nrange=1:3,criterion="BIC",maxit=5000,tol=1e-6){
  ## Perform the EM algorithm for different numbers of components
  ## Select best fit using the Bayesian Information Criterion (BIC)
  ## or the Akaike information criterion (AIC)
  i=1
  results = list()
  BIC_vec = c()
  AIC_vec = c()
  
  for (n in nrange){
    ## Initialise EM algorithm with values from kmeans clustering
    init = kmeans(x/size,n)
    prop_init = init$size/length(x)
    p_init = init$centers
    
    results[[i]] = em.algo(x,size,prop.vector_inits = prop_init,p.vector_inits=p_init,nclust=n,maxit,tol)
    BIC_vec = c(BIC_vec,results[[i]]$BIC)
    AIC_vec = c(AIC_vec,results[[i]]$AIC)
    i=i+1
  }
  if (criterion=="BIC"){
    results[[which.min(BIC_vec)]]$BIC_vec=BIC_vec
    return(results[[which.min(BIC_vec)]])
  }
  if (criterion=="AIC"){
    return(results[[which.min(AIC_vec)]])
  }
}


par(mfrow=c(2,1))
pdf("Private_muts_sig_MD7634.pdf",width=10,height=4)
plot_96profile(df = details%>%mutate(sampleID="MD3764")%>%filter(node%in%1:length(tree$tip.label)),colnames = c("sampleID","Chrom","Pos","Ref","Alt"),genomeFile = gfile)
dev.off()
pdf("Shared_muts_sig_MD7634.pdf",width=10,height=4)
plot_96profile(df = details%>%mutate(sampleID="MD3764")%>%filter(!node%in%1:length(tree$tip.label)),colnames = c("sampleID","Chrom","Pos","Ref","Alt"),genomeFile = gfile)
dev.off()




pdf("Private_muts_VAF.pdf",width=15,height=15)
par(mfrow=c(4,4))
sapply(1:16,function(i) vaf_density_plot_final(tree$tip.label[i],tree,filtered_muts$COMB_mats.tree.build,private_only = T))
dev.off()


COMB_mats=filtered_muts$COMB_mats.tree.build
all_private_vafs<-lapply(tree$tip.label,function(sample) {
  node <- which(tree$tip.label==sample)
  sample_muts <- COMB_mats$mat$mut_ref[COMB_mats$mat$node %in% node]
  COMB_mats$NR[COMB_mats$NR == 0] <- 1
  return((COMB_mats$NV/COMB_mats$NR)[sample_muts,sample])
})%>%unlist()




all_NV<-lapply(head(tree$tip.label,-1),function(sample) {
  node <- which(tree$tip.label==sample)
  sample_muts <- COMB_mats$mat$mut_ref[COMB_mats$mat$node %in% node]
  COMB_mats$NR[COMB_mats$NR == 0] <- 1
  return(COMB_mats$NV[sample_muts,sample])
})%>%unlist()

all_NR<-lapply(head(tree$tip.label,-1),function(Sample) {
  node <- which(tree$tip.label==Sample)
  sample_muts <- COMB_mats$mat$mut_ref[COMB_mats$mat$node %in% node]
  COMB_mats$NR[COMB_mats$NR == 0] <- 1
  sample_NR=COMB_mats$NR[sample_muts,Sample]
  names(sample_NR)<-sample_muts
  return(sample_NR)
})%>%unlist()

dens <- density(all_NV[all_NR>20]/all_NR[all_NR>20])
plot(dens,xlim = c(0,1),main=sample)
abline(v = dens$x[which.max(dens$y)])
text(0.7, max(dens$y) - 0.2, paste("Peak VAF dens=",round(dens$x[which.max(dens$y)], digits = 2)),col="red",cex = 0.7)



cutoff=20
temp=binom_mix(x=all_NV[all_NR>cutoff],size=all_NR[all_NR>cutoff],nrange=1:4)

clonal_mut_refs=names(all_NR)[all_NR>cutoff][which(temp$Which_cluster%in%2:3)]
subclonal_mut_refs=names(all_NR)[all_NR>cutoff][which(temp$Which_cluster==1)]


pdf("Clonal_vs_nonclonal_muts_sig_MD7634.pdf",width=10,height=8)
par(mfrow=c(2,1))
plot_96profile_mutref(clonal_mut_refs,genomeFile = gfile)
plot_96profile_mutref(subclonal_mut_refs,genomeFile = gfile)
dev.off()

