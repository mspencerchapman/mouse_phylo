#Mouse targeted sequencing data analysis

# ---------------------------------------------------#
# Import packages & custom scripts ----
# ---------------------------------------------------#

library(stringr)
library(ape)
library(seqinr)
library(data.table)
library(tidyr)
library(readr)
library(dplyr)
library(ggplot2)
library(readxl)

my_working_directory = ifelse(Sys.info()["sysname"] == "Darwin","~/R_work/mouse_phylo/","/lustre/scratch126/casm/team154pc/ms56/mouse_phylo/")

R_functions_dir=ifelse(Sys.info()["sysname"] == "Darwin","~/R_work/my_functions","/lustre/scratch126/casm/team154pc/ms56/my_functions")
tree_mut_dir=ifelse(Sys.info()["sysname"] == "Darwin","~/R_work/treemut","/lustre/scratch126/casm/team154pc/ms56/fetal_HSC/treemut")
R_function_files=list.files(R_functions_dir,pattern=".R",full.names = T)
sapply(R_function_files[-2],source)
setwd(tree_mut_dir); source("treemut.R");setwd(my_working_directory)

#Define function to import allele counter data similarly to cgpVAF NV and NR matrices
import_allele_counter_data=function(allele_counter_dir,suffix=".allelecounts.txt",bed_file=SNV_mut_refs,verbose=F){
  require(dplyr)
  require(stringr)
  output_files=list.files(path=allele_counter_dir,pattern=suffix,full.names = T)
  Samples=unlist(lapply(stringr::str_split(output_files,pattern="/"),function(x) gsub(suffix,"",tail(x,n=1))))
  if(verbose) {cat(paste("Importing",length(Samples),"allele counter files"),sep="\n")}
  NV_and_NR=lapply(output_files,function(file) {
    if(verbose) {cat(file,sep="\n")}
    dat=read.delim(file,stringsAsFactors = F)%>%
      dplyr::rename("Chrom"=X.CHR,"Pos"=POS)%>%
      left_join(bed_file,by=c("Chrom","Pos"))
    sample_NR=dat$Good_depth
    sample_NV=sapply(1:nrow(dat),function(i) {return(dat[i,paste0("Count_",dat$Alt[i])])})
    return(list(NV=sample_NV,NR=sample_NR))
  })
  
  #Get the standard set of mut_refs (in the same order as allele counter output) by looking at the first allele counter file
  if(verbose) {cat("Getting the 'mut_refs' in the correct order",sep="\n")}
  mut_refs=read.delim(output_files[1],stringsAsFactors = F)%>%
    dplyr::rename("Chrom"=X.CHR,"Pos"=POS)%>%
    left_join(bed_file,by=c("Chrom","Pos"))%>%
    tidyr::unite(col="mut_ref",Chrom,Pos,Ref,Alt,sep="-")%>%
    pull(mut_ref)
  
  #Bind the NV output and NR outputs into separate matrices
  if(verbose) {cat("Binding the NV and NR files into separate matrices",sep="\n")}
  NV=as.matrix(dplyr::bind_cols(lapply(NV_and_NR,function(list) list$NV)))
  NR=as.matrix(dplyr::bind_cols(lapply(NV_and_NR,function(list) list$NR)))
  
  #Name the columns as the samples and rows as the mutation references ('mut_refs')
  rownames(NV)=rownames(NR)<-mut_refs
  colnames(NV)=colnames(NR)<-Samples
  
  return(list(NV=NV,NR=NR))
}

calculate_cell_frac=function(NV,NR,sex="male") {
  #Remove 0's from the depth, to avoid dividing by 0
  NR[NR==0]<-1
  #For autosomal mutations, cell frac = NV/ (NR/2)
  cell_frac=NV/(NR/2)
  
  if(sex=="Male"|sex=="male"|sex=="M") {
    #Get vectors to select out the autosomal and XY chromosomal mutations
    XY_muts=grepl("X",rownames(NR))|grepl("Y",rownames(NR))
    #For XY mutations, cell frac = NV/NR -> replace these accordingly
    cell_frac[XY_muts,]<-(NV[XY_muts,])/(NR[XY_muts,])
  }
  #Cell frac cannot be greater than 1, therefore if comes out as > 1 (which can happen when dividing the NR by 2), coerce to 1
  cell_frac[cell_frac>1]<-1
  return(cell_frac)
}

get_expanded_clade_nodes=function(tree,height_cut_off=100,min_clonal_fraction=0.02,min_samples=1){
  nodeheights=nodeHeights(tree)
  
  #This pulls out nodes that fulfill on the criteria: branches cross the cut-off & contain the minimum proportion of samples
  nodes=tree$edge[,2][nodeheights[,1] < height_cut_off &
                        !nodeheights[,2] < height_cut_off &
                        sapply(tree$edge[,2],function(node) {length(getTips(tree,node))/length(tree$tip.label)})>min_clonal_fraction &
                        sapply(tree$edge[,2],function(node) {length(getTips(tree,node))})>=min_samples]
  df=data.frame(nodes=nodes,n_samples=sapply(nodes,function(node) {length(getTips(tree,node))}),MRCA_time=sapply(nodes,function(node) {nodeheight(tree,node)}),clonal_fraction=sapply(nodes,function(node) {length(getTips(tree,node))/length(tree$tip.label)}))
  return(df)
}

aggregate_cols_by_tissue=function(NV,NR) {
  Samples=colnames(NV)
  Tissue_IDs=unique(stringr::str_split(Samples,pattern="_",simplify=T)[,1])
  NV_agg=lapply(Tissue_IDs,function(ID){
    Tissue_NV<-rowSums(NV[,grepl(ID,colnames(NV)),drop=F])
    return(Tissue_NV)
  })%>%dplyr::bind_cols()%>%
    as.matrix()
  NR_agg=lapply(Tissue_IDs,function(ID){
    Tissue_NR<-rowSums(NR[,grepl(ID,colnames(NR)),drop=F])
    return(Tissue_NR)
  })%>%dplyr::bind_cols()%>%
    as.matrix()
  rownames(NV_agg)=rownames(NR_agg)=rownames(NV)
  colnames(NV_agg)=colnames(NR_agg)=Tissue_IDs
  return(list(NV=NV_agg,NR=NR_agg))
}

# ---------------------------------------------------#
# Set root directory & import metadata ----
# ---------------------------------------------------#

root_dir=ifelse(Sys.info()['sysname']=="Darwin","~/R_work/mouse_phylo/targeted/","/lustre/scratch126/casm/team154pc/ms56/mouse_phylo/targeted/")

my_theme<-theme(text = element_text(family="Helvetica"),
                axis.text = element_text(size = 5),
                axis.title = element_text(size=7),
                legend.text = element_text(size=5),
                legend.title = element_text(size=7),
                strip.text = element_text(size=7),
                legend.spacing = unit(1,"mm"),
                legend.key.size= unit(5,"mm"))

metadata<-readxl::read_excel(paste0("sample_metadata/targeted_metadata.xlsx"))%>%
  mutate(flow_marker=stringr::str_split(SUPPLIER_SAMPLE_ID,pattern="_",simplify=T)[,1],
         age=stringr::str_split(SUPPLIER_SAMPLE_ID,pattern="_",simplify=T)[,3],
         PDID=substr(INTERNAL_CASM_SAMPLE_NAME,1,6))

baitset_SNV_bed_file_path=paste0(root_dir,"baitset_SNVs.bed")
baitset_INDEL_bed_file_path=paste0(root_dir,"baitset_INDELs.bed")
allele_counter_dir=paste0(root_dir,"allele_counter/")
cgpvaf_SNV_file=paste0(root_dir,"merged_SNVs_targeted.tsv")
cgpvaf_INDEL_file=paste0(root_dir,"merged_indels_targeted.tsv")
plots_dir=paste0(root_dir,"plots/")
matrices_files=paste0(root_dir,"full_matrices.RDS")
matrices_files_cgpVAF=paste0(root_dir,"full_cgpvaf_matrices.RDS")
output_dir=paste0(root_dir,"output/")

if(Sys.info()['sysname'] == "Darwin"){
  tree_folder=paste0(my_working_directory,"/tree_files/")
  annotated_muts_folder=paste0(my_working_directory,"/annotated_muts/")
} else {
  tree_folder=paste0(my_working_directory,"filtering_runs/tree_files")
  annotated_muts_folder=paste0(my_working_directory,"filtering_runs/annotated_muts")
}


# ---------------------------------------------------#
# Create named vector to translate WGS/TGS IDs ----
# ---------------------------------------------------#

#There are two sets of MDIDs as the targeted sequencing samples were mistakenly given new IDs
Phylo_MDIDs=c("MD7634","MD7635")
MDIDs=c("MD7816","MD7817")
names(MDIDs)<-Phylo_MDIDs

# ---------------------------------------------------#
# Import the annotated mut files & and the trees ----
# ---------------------------------------------------#

#Set up the all.muts and all.trees objects - if previously run, read in the RDS object, otherwise import the details matrices from the annotated muts objects
cat("Loading the individual annotated mutation files",sep = "\n")
annotated_muts_paths=list.files(annotated_muts_folder,pattern="_vaf_post_mix_post_dup",full.names = T)
all.muts<-lapply(Phylo_MDIDs,function(ID){
  file_paths<-grep(ID,annotated_muts_paths,value = T)
  load(grep("_m40_",file_paths,invert=T,value = T));return(filtered_muts$COMB_mats.tree.build$mat)
})

tree_paths=list.files(tree_folder,pattern="_vaf_post_mix_post_dup.tree",full.names = T)
all.trees<-lapply(Phylo_MDIDs,function(ID){
  tree_paths<-grep(ID,tree_paths,value = T)
  tree<-read.tree(grep("_m40_",tree_paths,invert=T,value = T))
  return(tree)
})
names(all.muts)<-names(all.trees)<-Phylo_MDIDs

# ---------------------------------------------------#
# Import the allele counter targeted sequencing data ----
# ---------------------------------------------------#

if(!file.exists(matrices_files)){
  #Import the SNV bedfile and name the columns appropriately
  SNV_mut_refs=read.delim(baitset_SNV_bed_file_path,stringsAsFactors = F,header=F)
  colnames(SNV_mut_refs)<-c("Chrom","Pos","Ref","Alt")
  
  matrices=import_allele_counter_data(allele_counter_dir = allele_counter_dir, bed_file=SNV_mut_refs)
  saveRDS(matrices,file=matrices_files)
} else {
  matrices<-readRDS(matrices_files)
}

# ---------------------------------------------------#
# Import the cgpVAF targeted sequencing data ----
# ---------------------------------------------------#

if(T|!file.exists(matrices_files_cgpVAF)){
  cat("Importing cgpVAF files",sep="\n")
  matrices_cgpvaf=import_cgpvaf_SNV_and_INDEL(SNV_output_file = cgpvaf_SNV_file,INDEL_output_file = cgpvaf_INDEL_file)
  saveRDS(matrices_cgpvaf,file=matrices_files_cgpVAF)
} else {
  cat("Reading in previously saved cgpVAF files",sep="\n")
  matrices_cgpvaf<-readRDS(matrices_files_cgpVAF)
}

# ---------------------------------------------------#
# Combine data into complete 'targeted seq res' objects ----
# ---------------------------------------------------#

all_targeted_res=vector(mode = "list",length = length(MDIDs))
names(all_targeted_res)<-MDIDs

for(ID in Phylo_MDIDs){
  cat(ID,sep="\n")
  details=all.muts[[ID]]
  tree=all.trees[[ID]]
  
  #In details file, create column for whether it is included in the targeted sequencing baits
  details$targeted_tree = ifelse(details$mut_ref %in% rownames(matrices_cgpvaf$NV),"YES","NO")
  
  #Get the details matrix only including mutations in the targeted panel
  details_targ = details[details$targeted_tree=="YES",]
  
  #Get the tissue names from this mouse
  this_ID_tissues=metadata%>%filter(PDID==MDIDs[ID])%>%pull(INTERNAL_CASM_SAMPLE_NAME)%>%unique()
  
  #Remove tissue names that failed targeted sequencing
  this_ID_tissues<-this_ID_tissues[this_ID_tissues%in%colnames(matrices_cgpvaf[['NV']])]
  
  matrices_ID=lapply(matrices_cgpvaf[c("NV","NR")],function(mat) {
    mat[details_targ$mut_ref,this_ID_tissues]
  })
  
  all_targeted_res[[ID]]<-list(details_targ=details_targ,matrices=matrices_ID,smry=(metadata%>%filter(PDID==MDIDs[ID])))
}

# ---------------------------------------------------#
# Restructure the data ----
# ---------------------------------------------------#

## Nick's version of the Gibbs sampler requires data in a specific format
## Reforat and save for each individual sample

for(ID in Phylo_MDIDs){
  cat(paste("Starting analysis for",ID),sep="\n")
  details_targ<-all_targeted_res[[ID]]$details_targ
  
  this_mouse_IDs<-grep(MDIDs[ID],colnames(matrices_cgpvaf$NV),value=T)
  other_mouse_IDs<-grep(MDIDs[names(MDIDs)!=ID],colnames(matrices_cgpvaf$NV),value=T)
  
  mtr_other<-rowSums(matrices_cgpvaf[["NV"]][details_targ$mut_ref,other_mouse_IDs])
  depth_other<-rowSums(matrices_cgpvaf[["NR"]][details_targ$mut_ref,other_mouse_IDs])
  
  for(tissueID in this_mouse_IDs) {
    cat(tissueID,sep="\n")
    details_tissue<-cbind(details_targ[,c("Chrom","Pos","Ref","Alt","node")],data.frame(mtr=matrices_cgpvaf[["NV"]][details_targ$mut_ref,tissueID],
                                                                                        depth=matrices_cgpvaf[["NR"]][details_targ$mut_ref,tissueID],
                                                                                        mtr_other=mtr_other,
                                                                                        depth_other=depth_other))
    
    saveRDS(object=list(details=details_tissue,tree=all.trees[[ID]]),file = paste0(output_dir,tissueID,"_gibbs_info.RDS"))
  }
}

