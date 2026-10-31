#Mutational signatures of clonal vs sub-clonal private mutations, using the per-mutation posterior P(clonal) from the
#detection-aware mixture (burden_correction_<mouse>.Rds). Shared mutations (in vivo by definition) are the reference.
#Trinucleotide contexts are fetched once from the GRCm38 reference on the farm (read-only samtools faidx over ssh)
#and cached in results/signatures/contexts.Rds.
setwd("~/R_work/mouse_phylo/")
suppressMessages({library(dplyr);library(tidyr);library(ggplot2);library(MutationalPatterns)})
out_dir="burden_correction/results/signatures"; dir.create(out_dir,recursive=T,showWarnings=F)
ref_fa="/nfs/cancer_ref02/Mus_musculus/GRCm38/genome.fa"

#---------- Mutation sets ----------
muts=list()
for(m in c("MD7634","MD7635")){
  bc=readRDS(paste0("burden_correction/results/burden_correction_",m,".Rds"))
  load(paste0("annotated_muts/annotated_mut_set_",m,"_postMS_reduced_a_j_vaf_post_mix_post_dup"))
  d=filtered_muts$COMB_mats.tree.build$mat%>%filter(Mut_type=="SNV")
  tr=bc$tree; tip_nodes=which(tr$tip.label!="Ancestral")
  shared=d%>%filter(!node%in%tip_nodes)%>%transmute(mouse=m,mut_ref,set="shared",w=1)
  priv=bc$priv%>%transmute(mouse=m,mut_ref,sample,NV,NR,p_clonal)
  muts[[m]]=list(shared=shared,priv=priv)
}
all_refs=unique(c(unlist(lapply(muts,function(x) x$shared$mut_ref)),unlist(lapply(muts,function(x) x$priv$mut_ref))))

#---------- Trinucleotide contexts (cached) ----------
ctx_file=file.path(out_dir,"contexts.Rds")
ctx=if(file.exists(ctx_file)) readRDS(ctx_file) else data.frame(mut_ref=character(),tri=character())
need=setdiff(all_refs,ctx$mut_ref)
if(length(need)>0){
  p=do.call(rbind,strsplit(need,"-")); regions=paste0(p[,1],":",as.integer(p[,2])-1,"-",as.integer(p[,2])+1)
  tmp=tempfile(); writeLines(regions,tmp)
  fa=system(paste0("ssh -o BatchMode=yes farm22 '/software/sciops/pkgg/samtools/1.21/bin/samtools faidx ",ref_fa," -r /dev/stdin' < ",tmp),intern=TRUE)
  seqs=toupper(fa[!grepl("^>",fa)]); stopifnot(length(seqs)==length(need))
  ctx=rbind(ctx,data.frame(mut_ref=need,tri=seqs)); saveRDS(ctx,ctx_file)
}

#---------- 96-channel classification (pyrimidine reference, MutationalPatterns order) ----------
comp=c(A="T",C="G",G="C",T="A")
rc=function(x) sapply(strsplit(x,""),function(b) paste(rev(comp[b]),collapse=""))
classify=function(mut_ref,tri){p=do.call(rbind,strsplit(mut_ref,"-")); ref=p[,3]; alt=p[,4]
  stopifnot(all(substr(tri,2,2)==ref))
  flip=ref%in%c("G","A"); tri2=ifelse(flip,rc(tri),tri); ref2=ifelse(flip,comp[ref],ref); alt2=ifelse(flip,comp[alt],alt)
  paste0(substr(tri2,1,1),"[",ref2,">",alt2,"]",substr(tri2,3,3))}
ctx$channel=classify(ctx$mut_ref,ctx$tri)
channels=MutationalPatterns:::TRIPLETS_96
stopifnot(all(ctx$channel%in%channels))

#Weighted 96-profiles: shared; private weighted by P(clonal) / P(sub-clonal); hard calls (MAP, and high confidence)
sets=bind_rows(lapply(muts,function(x) bind_rows(
  x$shared,
  x$priv%>%transmute(mouse,mut_ref,set="private: clonal (weighted)",w=p_clonal),
  x$priv%>%transmute(mouse,mut_ref,set="private: sub-clonal (weighted)",w=1-p_clonal),
  x$priv%>%filter(p_clonal>0.5)%>%transmute(mouse,mut_ref,set="private: clonal (MAP)",w=1),
  x$priv%>%filter(p_clonal<=0.5)%>%transmute(mouse,mut_ref,set="private: sub-clonal (MAP)",w=1),
  x$priv%>%filter(p_clonal>0.9)%>%transmute(mouse,mut_ref,set="private: clonal (P>0.9)",w=1),
  x$priv%>%filter(p_clonal<0.1)%>%transmute(mouse,mut_ref,set="private: sub-clonal (P<0.1)",w=1))))%>%
  left_join(ctx%>%dplyr::select(mut_ref,channel),by="mut_ref")
set_levels=c("shared","private: clonal (weighted)","private: sub-clonal (weighted)","private: clonal (MAP)","private: sub-clonal (MAP)",
             "private: clonal (P>0.9)","private: sub-clonal (P<0.1)")
sets$set=factor(sets$set,levels=set_levels)
prof=sets%>%mutate(col=paste(mouse,set,sep=" | "))%>%group_by(col,channel)%>%summarise(n=sum(w),.groups="drop")%>%
  pivot_wider(names_from=col,values_from=n,values_fill=0)
M=as.matrix(prof[,-1]); rownames(M)=prof$channel; M=M[channels,,drop=F]
M=M[,order(match(sub(".*\\| ","",colnames(M)),set_levels),sub(" \\|.*","",colnames(M)))]
both=sapply(set_levels,function(s) rowSums(M[,sub(".*\\| ","",colnames(M))==s,drop=F])); colnames(both)=paste("both |",set_levels)
M_all=cbind(M,both)
write.csv(M_all,file.path(out_dir,"profiles_96.csv"))
cat("Mutations per set (weighted counts):\n"); print(round(colSums(M_all),0))

#---------- Signature fitting ----------
sigs=get_known_signatures(muttype="snv",source="COSMIC",genome="mm10"); rownames(sigs)=channels

#Restricted fit: signatures seen in mouse/human HSCs and in culture (SBS1, SBS5, SBS18), plus a check with SBS40
fit_restricted=function(sel) {f=fit_to_signatures(M_all,sigs[,sel]); round(t(t(f$contribution)/colSums(f$contribution)),3)}
r1=fit_restricted(c("SBS1","SBS5","SBS18")); r2=fit_restricted(c("SBS1","SBS5","SBS18","SBS40"))
recon=fit_to_signatures(M_all,sigs[,c("SBS1","SBS5","SBS18")])$reconstructed
cos_fit=sapply(seq_len(ncol(M_all)),function(i) cos_sim(M_all[,i],recon[,i]))

#Strict fit over all COSMIC signatures (backwards elimination)
strict=fit_to_signatures_strict(M_all,sigs,max_delta=0.004)$fit_res$contribution
strict=round(t(t(strict)/colSums(strict)),3); strict=strict[rowSums(strict)>0,,drop=F]

#Key channel fractions: C>T at CpG (SBS1 hallmark) and C>A
cpg=grepl("\\[C>T\\]G",channels); ca=grepl("C>A",channels)
summ=data.frame(set=colnames(M_all),n=round(colSums(M_all)),CpG_CtoT=round(colSums(M_all[cpg,])/colSums(M_all),3),
  CtoA=round(colSums(M_all[ca,])/colSums(M_all),3),t(r1),cosine_SBS1_5_18_fit=round(cos_fit,3),
  cos_to_shared=round(sapply(colnames(M_all),function(c){mm=sub(" \\|.*","",c); cos_sim(M_all[,c],M_all[,paste(mm,"| shared")])}),3),row.names=NULL)
write.csv(summ,file.path(out_dir,"signature_summary.csv"),row.names=F); write.csv(strict,file.path(out_dir,"strict_fit_all_cosmic.csv"))
write.csv(r2,file.path(out_dir,"restricted_fit_SBS1_5_18_40.csv"))
print(summ,row.names=F); cat("\nStrict fit (all COSMIC), proportions:\n"); print(strict)
cat("\nRestricted fit with SBS40:\n"); print(t(r2))

#Bootstrap CIs (multinomial resampling) for the high-confidence and MAP hard-call sets, both mice combined
hard=grep("^both \\|.*(MAP|P>0.9|P<0.1)|^both \\| shared",colnames(M_all),value=TRUE)
set.seed(1)
boot=fit_to_signatures_bootstrapped(round(M_all[,hard]),sigs[,c("SBS1","SBS5","SBS18")],n_boots=200,method="regular")
bt=as.data.frame(boot)%>%mutate(sample=sub("_[0-9]+$","",rownames(boot)))%>%pivot_longer(-sample,names_to="signature")%>%
  group_by(sample)%>%mutate(tot=ave(value,rep(seq_len(n()/3),each=3),FUN=sum))
bt=as.data.frame(boot); bt$sample=sub("_[0-9]+$","",rownames(bt)); bt[,1:3]=bt[,1:3]/rowSums(bt[,1:3])
ci=bt%>%pivot_longer(-sample,names_to="signature")%>%group_by(sample,signature)%>%
  summarise(median=round(median(value),3),lo=round(quantile(value,0.025),3),hi=round(quantile(value,0.975),3),.groups="drop")
write.csv(ci,file.path(out_dir,"bootstrap_SBS1_5_18.csv"),row.names=F); cat("\nBootstrap 95% CIs (proportions):\n"); print(as.data.frame(ci))

#---------- Plots ----------
pdf(file.path(out_dir,"signatures_by_clonality.pdf"),width=12,height=10)
print(plot_96_profile(M_all[,grep("^both",colnames(M_all))],ymax=0.15,condensed=TRUE))
print(plot_96_profile(M_all[,grep("^MD7634",colnames(M_all))],ymax=0.15,condensed=TRUE))
print(plot_96_profile(M_all[,grep("^MD7635",colnames(M_all))],ymax=0.15,condensed=TRUE))
print(as.data.frame(t(r1))%>%mutate(set=rownames(.))%>%pivot_longer(-set)%>%
  ggplot(aes(value,factor(set,levels=rev(colnames(M_all))),fill=name))+geom_col()+theme_classic()+labs(x="Proportion",y=NULL,fill="Signature"))
dev.off()
#Per-mutation assignments, for other analyses
write.csv(bind_rows(lapply(muts,function(x) x$priv))%>%left_join(ctx%>%dplyr::select(mut_ref,channel),by="mut_ref")%>%
  mutate(assignment=ifelse(p_clonal>0.5,"clonal","sub-clonal"),confident=p_clonal>0.9|p_clonal<0.1),
  file.path(out_dir,"private_mutation_assignments.csv"),row.names=F)
