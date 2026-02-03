#========================================#
# Load packages (and install if they are not installed yet) ####
#========================================#
cran_packages=c("dplyr","ggplot2","tidyr","readxl","readr")
bioconductor_packages=c()

for(package in cran_packages){
  if(!require(package, character.only=T,quietly = T, warn.conflicts = F)){
    install.packages(as.character(package),repos = "http://cran.us.r-project.org")
    library(package, character.only=T,quietly = T, warn.conflicts = F)
  }
}

#========================================#
# WGS: Import data and combine ####
#========================================#

root_dir="~/R_work/mouse_phylo/"
setwd(root_dir)
shipping_info_dat<-readxl::read_excel(path=paste0(root_dir,"sample_metadata/shipping_info_doc.xlsx"))
canapps_WGS_info<-readxl::read_excel(path=paste0(root_dir,"canapps_info/Cancer_Pipeline_Reports_3357.xls"),skip = 2)

comb_dat<-left_join(shipping_info_dat,canapps_WGS_info,by=c("INTERNAL_CASM_SAMPLE_NAME"="Sample"))

#Summary of number of colonies from each anatomical location
comb_dat%>%
  mutate(MOUSE_ID=substr(INTERNAL_CASM_SAMPLE_NAME,1,6))%>%
  filter(`Seq X`>4)%>%
  group_by(MOUSE_ID,TISSUE_PHENOTYPE)%>%
  dplyr::summarise(n=n())

comb_dat%>%
  mutate(MOUSE_ID=substr(INTERNAL_CASM_SAMPLE_NAME,1,6))%>%
  filter(`Seq X`>4)%>%
  group_by(MOUSE_ID)%>%
  dplyr::summarise(n=n())

comb_dat%>%
  mutate(MOUSE_ID=substr(INTERNAL_CASM_SAMPLE_NAME,1,6))%>%
  filter(`Seq X`>4)%>%
  ggplot(aes(x=TISSUE_PHENOTYPE,fill=TISSUE_PHENOTYPE))+
  geom_bar()+
  theme_classic()+
  facet_grid(~MOUSE_ID)+
  labs(x="Anatomical location")

#========================================#
# TGS: Import data and combine ####
#========================================#

root_dir="~/R_work/mouse_phylo/"
setwd(root_dir)
header=readxl::read_excel(path = "sample_metadata/CASM_STS_Sample_MANIFEST_Lily_Cabrera_York_July_2024_x2_plates.xlsx",sheet = 1,range = "A1:AL1",col_names = T)
dat<-readxl::read_excel(path = "sample_metadata/CASM_STS_Sample_MANIFEST_Lily_Cabrera_York_July_2024_x2_plates.xlsx",sheet = 1,skip=13,col_names = F)
colnames(dat)<-colnames(header)
dat<-dat%>%
  filter(!is.na(VOLUME_UL))

#========================================#
# TGS: Plot summaries ####
#========================================#

dat%>%
  mutate(total=VOLUME_UL*CONCENTRATION_NG_UL)%>%
  ggplot(aes(x=total))+
  geom_histogram()+
  theme_bw()+
  scale_x_log10()

#If we do every 5ng aliquot separately, how many samples will be submitted in total?
max_ng_per_sample=5
n_colonies_resequenced=30
n_reactions=32
total_samps<-dat%>%
  mutate(total=VOLUME_UL*CONCENTRATION_NG_UL)%>%
  mutate(nsamp=ceiling(total/max_ng_per_sample))%>%
  summarise(total_samp=sum(nsamp))
(total_samps$total_samp+n_colonies_resequenced)/n_reactions
View(dat)
#========================================#
#Estimated coverage
#========================================#
L=200 #total read length if doing 100bp PE reads
N=1e10/8 #Number of reads in the NovaSeq 10B flow cell, divided by the number of lanes => the number of reads per lane
G=2.85e6 #Size of the target region of the bait set
On_target=0.3 #Conservative (hopefully) estimate of the on target %
NSample=total_samps$total_samp+n_colonies_resequenced

cov_per_sample=(On_target*L*N/G)/NSample
cov_per_sample


