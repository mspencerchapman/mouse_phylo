#Set the project and the bait set number
#BAITSET=1, TARGET_PROJECT=2737,BEDFILE_ID=S3324404
#BAITSET=2, TARGET_PROJECT=2712,BEDFILE_ID=S3358882
#BAITSET=3, TARGET_PROJECT=3035,BEDFILE_ID=S3401023
TARG_PROJECT=3523
STUDY_DIR=/lustre/scratch126/casm/team154pc/ms56/mouse_phylo/targeted
INSILICO_NORMAL_ID=MDGRCm38is
ORGANISM=mouse
ORGANISM_LATIN=Mus_musculus
GENOME_REF=GRCm38
GENOME_FILE=/nfs/cancer_ref02/${ORGANISM}/${GENOME_REF}/genome.fa
HIGH_DEPTH_BEDFILE=/lustre/scratch124/casm/team78pipelines/canpipe/live/ref/${ORGANISM_LATIN}/${GENOME_REF}/shared/extremedepth.bed.gz

cd ${STUDY_DIR}
mkdir -p allele_counter
#Set of the allele counter analysis using the Fluidigm module version of allele counter (iterates through each sample one at a time)
module purge
module load Fluidigm
bsub -o $PWD/log.%J \
    -e $PWD/err.%J \
    -q yesterday \
    -G team78-grp \
    -R 'select[mem>=16000] span[hosts=1] rusage[mem=16000]' \
    -M16000 \
    -n1 \
    -J "Allele_Ct" \
    ./FluidigmAlleleCounter.sh $TARG_PROJECT

#Get the reference files required for picard tools metrics analysis
COVERED_BED_FILE=${STUDY_DIR}/baitset_info/merged_probe_file_shareable_Sanger_Mouse_phylo_REDUCED_v2_TE-99621817_mm10_low.bed
COVERED_BED_FOR_PICARD=${STUDY_DIR}/baitset_info/Covered_bed_files_combined.bed

cat $COVERED_BED_FILE |cut -f1-3|sed 's/chr//g'>$COVERED_BED_FOR_PICARD

module load picard-tools/3.1.0
bsub -o $PWD/log.%J \
    -e $PWD/err.%J \
    -q normal \
    -G team78-grp \
    -R 'select[mem>=4000] span[hosts=1] rusage[mem=4000]' \
    -M4000 \
    -n1 \
    -J "bedlist" \
    java -jar $PICARD BedToIntervalList \
      I=$COVERED_BED_FOR_PICARD \
      O=${STUDY_DIR}/baitset_info/Covered_bed_files_combined.interval_list \
      SD=${GENOME_FILE}.dict

TARGETS_BED_FILE=${STUDY_DIR}/baitset_info/all_target_segments_covered_by_probes_Sanger_Mouse_phylo_REDUCED_v2_TE-99621817_mm10_low.bed
TARGETS_BED_FOR_PICARD=${STUDY_DIR}/baitset_info/Targets_from_covered.bed

cat $COVERED_BED_FILE |cut -f4|sed 's/chr//g'|sed 's/,/\n/g'|sed 's/:/\t/g'|sed 's/_/\t/g'|awk -F'\t' 'NF==3 {print}'|sed 's/x/X/' >$TARGETS_BED_FOR_PICARD

bsub -o $PWD/log.%J \
    -e $PWD/err.%J \
    -q normal \
    -G team78-grp \
    -R 'select[mem>=4000] span[hosts=1] rusage[mem=4000]' \
    -M4000 \
    -n1 \
    -J "bedlist" \
    java -jar $PICARD BedToIntervalList \
      I=$TARGETS_BED_FOR_PICARD \
      O=${STUDY_DIR}/baitset_info/Targets_bed_files_combined.interval_list \
      SD=${GENOME_FILE}.dict

#Get list of samples for staging bams
cd $STUDY_DIR
ls /nfs/cancer_ref01/nst_links/live/$TARG_PROJECT>samples.txt
module load dataImportExport
mkdir -p bams
mkdir -p logs
lfs setstripe bams -c -1
stageBam.pl -p $TARG_PROJECT -s samples.txt -o bams -lo logs

#
mkdir -p ${STUDY_DIR}/metrics_files
./hs_metrics_loop.sh $TARG_PROJECT
#Combine individual metrics output into table with headers

ls *.hs_metrics.txt|cut -c1-15>samples.txt
sed -n 7p *.hs_metrics.txt|head -n1>headings.txt
for FILE in $(ls *.hs_metrics.txt); do sed -n 8p $FILE; done>metrics_combined.txt
paste samples.txt metrics_combined.txt > sample_metrics.txt
cat headings.txt sample_metrics.txt>metrics_table.txt

##CGPVaf analysis
module load cgpVAFcommand
CREATE_SPLIT_CONFIG_SCRIPT=/lustre/scratch126/casm/team154pc/ms56/my_programs/create_split_config_ini.R
SNV_BEDFILE_NAME=${STUDY_DIR}/baitset_SNVs.bed
INDEL_BEDFILE_NAME=${STUDY_DIR}/baitset_INDELs.bed

#SNV analysis
mkdir -p $STUDY_DIR/cgpVAF/SNV_analysis
cd $STUDY_DIR/cgpVAF/SNV_analysis
mkdir -p output
echo "3"|createVafCmd.pl \
    -pid $TARG_PROJECT \
    -o output \
    -g ${GENOME_FILE} \
    -hdr $HIGH_DEPTH_BEDFILE \
    -mq 30  \
    -bo 1  \
    -b $SNV_BEDFILE_NAME

cd $STUDY_DIR/cgpVAF/SNV_analysis/output
Rscript $CREATE_SPLIT_CONFIG_SCRIPT -p $TARG_PROJECT -b 10 -n $INSILICO_NORMAL_ID -s ${STUDY_DIR}/samples.txt
cd $STUDY_DIR/cgpVAF/SNV_analysis

echo "3"|createVafCmd.pl \
    -pid $TARG_PROJECT  \
    -o output \
    -i output/${TARG_PROJECT}_cgpVafConfig_split.ini \
    -g ${GENOME_FILE} \
    -hdr $HIGH_DEPTH_BEDFILE \
    -mq 30  \
    -bo 1  \
    -b $SNV_BEDFILE_NAME


sed -e 's/\%5/\%50/g' run_bsub.sh >run_bsub_updated.sh
bash run_bsub_updated.sh

cd $STUDY_DIR/cgpVAF/SNV_analysis/output/output/$INSILICO_NORMAL_ID/snp
ls *_vaf.tsv > files

#for first file
cut -f 3,4,5,6,24,26,39,41,54,56,69,71,84,86,99,101,114,116,129,131,144,146,159,161,174,176,189,191,204,206,219,221,234,236,249,251,264,266,279,281,294,296,309,311,324,326,339,341,354,356,369,371,384,386,399,401,414,416,429,431,444,446,459,461,474,476,489,491,504,506 $(sed -n '1p' files) > temp.1.cut   #Chr, Pos, Ref, Alt, + MTR, DEP for all samples (max 11 samples in one file)

#for subsequent files (exclude Chr, Pos, Ref, Alt and PDv37is)
for FILE in $(tail -n+2 files); do    
    if [ -s temp.$FILE ]
    then
        echo "temp file temp.$FILE already exists. Moving onto next file..."
    else
        echo "temp file temp.$FILE does not yet exist, will be created"
        cut -f 39,41,54,56,69,71,84,86,99,101,114,116,129,131,144,146,159,161,174,176,189,191,204,206,219,221,234,236,249,251,264,266,279,281,294,296,309,311,324,326,339,341,354,356,369,371,384,386,399,401,414,416,429,431,444,446,459,461,474,476,489,491,504,506 $FILE >  temp.$FILE.cut
    fi       
done

#Remove empty rows where header was with awk
for FILE in $(ls temp.*); do
    echo $FILE
    awk 'NF' $FILE > output.$FILE
    rm $FILE
done

#Concatenate output files to one merged file & move to the root directory
paste output.* > merged_SNVs_targeted.tsv && rm output.*

mv $STUDY_DIR/cgpVAF/SNV_analysis/output/output/${INSILICO_NORMAL_ID}/snp/merged_SNVs_targeted.tsv $STUDY_DIR
cd $STUDY_DIR

#----------------INDEL analysis----------------
mkdir -p $STUDY_DIR/cgpVAF/INDEL_analysis
cd $STUDY_DIR/cgpVAF/INDEL_analysis
mkdir -p output
echo "1"|createVafCmd.pl \
    -pid $TARG_PROJECT \
    -o output \
    -g ${GENOME_FILE} \
    -hdr /lustre/scratch124/casm/team78pipelines/canpipe/live/ref/Mus_musculus/${GENOME_REF}/shared/extremedepth.bed.gz \
    -mq 30  \
    -bo 1  \
    -b $INDEL_BEDFILE_NAME

cd $STUDY_DIR/cgpVAF/INDEL_analysis/output
Rscript $CREATE_SPLIT_CONFIG_SCRIPT -p $TARG_PROJECT -b 10 -n $INSILICO_NORMAL_ID -s ${STUDY_DIR}/samples.txt
cd $STUDY_DIR/cgpVAF/INDEL_analysis

echo "1"|createVafCmd.pl \
    -pid $TARG_PROJECT  \
    -o output \
    -i output/${TARG_PROJECT}_cgpVafConfig_split.ini \
    -g ${GENOME_FILE} \
    -hdr $HIGH_DEPTH_BEDFILE \
    -mq 30  \
    -bo 1  \
    -b $INDEL_BEDFILE_NAME

sed -e 's/\%5/\%50/g' run_bsub.sh >run_bsub_updated.sh
bash run_bsub_updated.sh

cd $STUDY_DIR/cgpVAF/INDEL_analysis/output/output/$INSILICO_NORMAL_ID/indel
ls *_vaf.tsv > files

#for first file
cut -f 3,4,5,6,16,18,25,27,34,36,43,45,52,54,61,63,70,72,79,81,88,90,97,99,106,108,115,117,124,126,133,135,142,144,151,153,160,162,169,171,178,180,187,189,196,198,205,207,214,216,223,225,232,234,241,243,250,252,259,261,268,270,277,279,286,288,295,297,304,306,313,315 $(sed -n '1p' files) > temp.1   #Chr, Pos, Ref, Alt, + MTR, DEP for all samples (max 11 samples in one file)

#for subsequent files (exclude Chr, Pos, Ref, Alt and PDv37is)
for FILE in $(tail -n+2 files); do
    cut -f 25,27,34,36,43,45,52,54,61,63,70,72,79,81,88,90,97,99,106,108,115,117,124,126,133,135,142,144,151,153,160,162,169,171,178,180,187,189,196,198,205,207,214,216,223,225,232,234,241,243,250,252,259,261,268,270,277,279,286,288,295,297,304,306,313,315 $FILE >  temp.$FILE
done

#Remove empty rows where header was with awk
for FILE in $(ls temp.*); do
    awk 'NF' $FILE > output.$FILE
    rm $FILE
done

#Concatenate output files to one merged file & move to the root directory
paste output.* > merged_indels_targeted.tsv && rm output.*
mv $STUDY_DIR/cgpVAF/INDEL_analysis/output/output/${INSILICO_NORMAL_ID}/indel/merged_indels_targeted.tsv $STUDY_DIR/
cd $STUDY_DIR


#Run the Gibbs Sampler to get posterior distributions
cd /lustre/scratch126/casm/team154pc/ms56/mouse_phylo/targeted/GibbsSampler

prepend_path LD_LIBRARY_PATH /software/isg/languages/R/4.4.0/exec/lib/R/lib
prepend_path LD_LIBRARY_PATH /nfs/users/nfs_m/ms56/R/x86_64-pc-linux-gnu-library/4.4

module load julia/1.9.2
julia wrap_gibbs.jl test

#RUNNING GIBBS SAMPLER LOCALLY
cd /Users/ms56/R_work/mouse_phylo/targeted/GibbsSampler
SAMPLES=$(cat ../samples.txt)
for ID in $SAMPLES; do julia wrap_gibbs.jl $ID; done

samples<-readLines("../samples.txt")
for (ID in samples[4:123]) {system(paste("julia wrap_gibbs.jl",ID))}