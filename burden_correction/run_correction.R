#Apply the detection-aware branch-length correction to both mice
#Outputs (burden_correction/results/): per-sample burdens, fitted parameters, corrected trees, diagnostic plots
setwd("~/R_work/mouse_phylo/")
source("burden_correction/R/burden_correction.R")
library(ggplot2)
out_dir="burden_correction/results"; dir.create(out_dir,recursive=T,showWarnings=F)
md=suppressMessages(readxl::read_excel("sample_metadata/CAMPUOY000001_REUPLOAD_20231116_4xplates_York_Lily_CASM_STS_Sample_MANIFEST_Rec20231206_WithCGPIDs20231219.xlsx"))%>%
  dplyr::select(sample=INTERNAL_CASM_SAMPLE_NAME,tissue=TISSUE_PHENOTYPE)

all_samples=list(); all_fits=list(); diag=list()
for(sample_ID in c("MD7634","MD7635")){
  load(paste0("annotated_muts/annotated_mut_set_",sample_ID,"_postMS_reduced_a_j_vaf_post_mix_post_dup"))
  tree=read.tree(paste0("tree_files/tree_",sample_ID,"_postMS_reduced_a_j_vaf_post_mix_post_dup.tree"))
  bc=correct_branch_lengths(tree,filtered_muts$COMB_mats.tree.build,sex="male")
  write.tree(bc$tree,file.path(out_dir,paste0("tree_",sample_ID,"_corrected.tree")))
  saveRDS(bc,file.path(out_dir,paste0("burden_correction_",sample_ID,".Rds")))
  #Pipeline tree (sensitivity-adjusted, no sub-clonal removal) for comparison
  pipeline_private=sapply(bc$samples$sample,function(s) tree$edge.length[tree$edge[,2]==which(tree$tip.label==s)])
  pipeline_total=sapply(bc$samples$sample,function(s){tip=which(tree$tip.label==s); sum(tree$edge.length[tree$edge[,2]%in%c(tip,anc_nodes(tree,tip))])})
  all_samples[[sample_ID]]=bc$samples%>%mutate(mouse=sample_ID,raw_private=n_obs,pipeline_private=pipeline_private,pipeline_total=pipeline_total)%>%left_join(md,by="sample")
  all_fits[[sample_ID]]=data.frame(mouse=sample_ID,pc=bc$fit$pc,ps=bc$fit$ps,rho_clonal=bc$fit$rho[1],rho_subclonal=bc$fit$rho[2],delta_NV3=bc$fit$delta[1],delta_NV4=bc$fit$delta[2],
                                   mean_w_clonal=mean(bc$fit$w),convergence=bc$fit$convergence)
  #Diagnostic: true-space clonal weight refitted within depth bins (flat if the detection model is right)
  p=bc$priv; f=bc$fit; pcm=function(x) pmin(p$cf*x,0.999)
  fc=dbb(p$NV,p$NR,pcm(f$pc),f$rho[1]); fs=dbb(p$NV,p$NR,pcm(f$ps),f$rho[2])
  Sc=pass_prob(p$NR,pcm(f$pc),f$delta,rho=f$rho[1]); Ss=pass_prob(p$NR,pcm(f$ps),f$delta,rho=f$rho[2])
  b=cut(p$NR,c(0,10,14,17,20,25,Inf))
  diag[[sample_ID]]=data.frame(mouse=sample_ID,depth_bin=levels(b),n=as.numeric(table(b)),
    w_clonal=sapply(split(seq_len(nrow(p)),b),function(i) optimize(function(w) -sum(log(w*fc[i]+(1-w)*fs[i])-log(w*Sc[i]+(1-w)*Ss[i])),c(1e-6,1-1e-6))$minimum))
}
samples=bind_rows(all_samples); fits=bind_rows(all_fits); diag=bind_rows(diag)
write.csv(samples,file.path(out_dir,"corrected_burdens.csv"),row.names=F)
write.csv(fits,file.path(out_dir,"fitted_parameters.csv"),row.names=F)
write.csv(diag,file.path(out_dir,"clonal_weight_by_depth.csv"),row.names=F)
print(fits,row.names=F); print(diag%>%mutate(w_clonal=round(w_clonal,3)),row.names=F)
print(as.data.frame(samples%>%group_by(mouse,tissue)%>%
  summarise(n=n(),raw_private_SNVs=round(mean(raw_private),1),pipeline_private=round(mean(pipeline_private),1),corrected_private=round(mean(N_clonal),1),
            pipeline_total=round(mean(pipeline_total),1),corrected_total=round(mean(total_burden),1),sd_total=round(sd(total_burden),1),.groups="drop")))
for(m in unique(samples$mouse)){x=samples[samples$mouse==m,]
  cat(m,": caller sensitivity",round(median(x$caller_sensitivity),3),"| clonal sensitivity",round(median(x$sensitivity),3),
      "| private: raw",round(mean(x$raw_private),1),"pipeline",round(mean(x$pipeline_private),1),"corrected",round(mean(x$N_clonal),1),
      "| total: pipeline",round(mean(x$pipeline_total),1),"corrected",round(mean(x$total_burden),1),"(SEM",round(sd(x$total_burden)/sqrt(nrow(x)),2),")",
      "| Spearman corrected private vs depth",round(cor(x$N_clonal,x$median_depth,method="s"),2),"\n")}

#Plots: raw vs corrected private burden, total burden by tissue
pdf(file.path(out_dir,"burden_correction_plots.pdf"),width=8,height=4)
print(ggplot(samples,aes(raw_private,N_clonal,col=tissue))+geom_point(size=1)+geom_abline(linetype=2)+facet_wrap(~mouse)+
  theme_classic()+labs(x="Observed private SNVs",y="Corrected clonal private burden",col="Tissue"))
print(ggplot(samples,aes(median_depth,N_clonal,col=tissue))+geom_point(size=1)+facet_wrap(~mouse)+
  theme_classic()+labs(x="Median sample depth",y="Corrected clonal private burden",col="Tissue"))
print(ggplot(samples,aes(tissue,total_burden,fill=tissue))+geom_boxplot(outlier.size=0.5)+facet_wrap(~mouse)+
  theme_classic()+labs(x=NULL,y="Corrected total burden (shared + private)")+theme(legend.position="none"))
print(ggplot(diag,aes(depth_bin,w_clonal,group=mouse,col=mouse))+geom_point()+geom_line()+ylim(0,1)+
  theme_classic()+labs(x="Depth of private mutation",y="Clonal weight (refitted per bin)"))
dev.off()
