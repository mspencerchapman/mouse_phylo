#Load the required modules
module load cgpVAFcommand
module load mpboot
#module load perl-5.38.0
module load vagrent
module load pcap-core
module load picard-tools
module load canned-queries-client
module load bwa-0.7.17

#Set variables
STUDY_ID="mouse_phylo"
EXP_ID="MD7634" #others are MD7635/4
ALL_PROJECT_NUMBERS=3357
INSILICO_NORMAL_ID=MDGRCm38is
INSILICO_NORMAL_PROJECT=1670
DONOR_AGE=2

#Reference file paths
GENOME_FILE=/lustre/scratch124/casm/team78pipelines/canpipe/live/ref/Mus_musculus/GRCm38/genome.fa
HIGH_DEPTH_REGIONS=/lustre/scratch124/casm/team78pipelines/canpipe/live/ref/Mus_musculus/GRCm38/shared/extremedepth.bed.gz
IMPORT_NEW_SAMPLES_SCRIPT=/lustre/scratch126/casm/team154pc/ms56/my_programs/import_new_samples_only.R
CREATE_SPLIT_CONFIG_SCRIPT=/lustre/scratch126/casm/team154pc/ms56/my_programs/create_split_config_ini.R
MUTATION_PARAMETER_SCRIPT=/lustre/scratch126/casm/team154pc/ms56/my_programs/Mutation_filtering_get_parameters_MOUSE.R
SENSITIVITY_ANALYSIS_SCRIPT=/lustre/scratch119/casm/team154pc/ms56/my_programs/Sensitivity_analysis_from_SNPs.R
GENES_TO_ANNOTATE=/lustre/scratch126/casm/team154pc/ms56/reference_files/chip_drivers.txt
TREE_BUILDING_SCRIPT=

#Secondary variable
PD_NUMBERS=$EXP_ID
STUDY_DIR=/lustre/scratch126/casm/team154pc/ms56/${STUDY_ID}
MATS_AND_PARAMS_DIR=${STUDY_DIR}/filtering_runs/mats_and_params
ROOT_DIR=/lustre/scratch126/casm/team154pc/ms56/${STUDY_ID}/${EXP_ID}
SNV_BEDFILE_NAME=${STUDY_ID}_${EXP_ID}_caveman.bed
INDEL_BEDFILE_NAME=${STUDY_ID}_${EXP_ID}_pindel.bed
MS_FILTERED_BEDFILE_NAME=${STUDY_ID}_${EXP_ID}_postMS_SNVs.bed

RUN_ID=$EXP_ID
RUN_ID_M=${EXP_ID}_m40
RUN_ID_TB=${RUN_ID_M}_postMS_reduced
IFS=',' read -r -a PROJECT_ARRAY <<< "$ALL_PROJECT_NUMBERS"

#--------------------Submit LCM filtering jobs-----------------------------
module load openjdk-16.0.2
module load samtools-1.19/python-3.12.0
mkdir -p $ROOT_DIR/MS_filters/output_files
cd $ROOT_DIR/MS_filters
/lustre/scratch126/casm/team154pc/ms56/my_programs/Submitting_Mathijs_filters_jobs.R -p $ALL_PROJECT_NUMBERS -s $PD_NUMBERS -o $ROOT_DIR/MS_filters/output_files -q long

#--------------------SNV analysis-----------------------------
mkdir -p $ROOT_DIR/caveman_raw/caveman_pileup/output; cd $ROOT_DIR/caveman_raw
Rscript $IMPORT_NEW_SAMPLES_SCRIPT -p $ALL_PROJECT_NUMBERS -s $PD_NUMBERS -t SNVs
cut -f 1,2,4,5 *.caveman_c.annot.vcf_pass_flags|sort|uniq>caveman_pileup/$SNV_BEDFILE_NAME
cd $ROOT_DIR/caveman_raw/caveman_pileup

#Run createVafCmd - with input value '3' (this is the option for selecting the caveman files as input)
for PROJECT_NUMBER in "${PROJECT_ARRAY[@]}"; do
    echo "3"|createVafCmd.pl -pid $PROJECT_NUMBER  -o output -g $GENOME_FILE -hdr $HIGH_DEPTH_REGIONS -mq 30  -bo 1  -b $SNV_BEDFILE_NAME
done

#Then run the "create_split_config_ini.R" script in the output folder
PROJECT_NUMBER=${PROJECT_ARRAY[0]}
cd $ROOT_DIR/caveman_raw/caveman_pileup/output
Rscript $CREATE_SPLIT_CONFIG_SCRIPT -p $PROJECT_NUMBER
cd $ROOT_DIR/caveman_raw/caveman_pileup

#Then re-run the createVafCmd.pl script with the new config file  - with input value '3' (this is the option for selecting the caveman files as input)
echo "3"|createVafCmd.pl -pid $PROJECT_NUMBER  -o output -i output/${PROJECT_NUMBER}_cgpVafConfig_split.ini -g $GENOME_FILE -hdr $HIGH_DEPTH_REGIONS -mq 30  -bo 1  -b $SNV_BEDFILE_NAME

#Update the run_bsub.sh command to allow more jobs in the array to run together and get more memory
sed -e 's/\%5/\%50/g;s/2000/4000/g;s/500/1000/g;' run_bsub.sh >run_bsub_updated.sh

#Run this if need to switch to the long queue (not normally necessary)
#sed -i -e 's/normal/long/g' run_bsub_updated.sh

bash run_bsub_updated.sh

#--------------------Indel analysis-----------------------------
mkdir -p $ROOT_DIR/pindel_raw/pindel_pileup/output; cd $ROOT_DIR/pindel_raw
Rscript $IMPORT_NEW_SAMPLES_SCRIPT -p $ALL_PROJECT_NUMBERS -s $PD_NUMBERS -t indels
cut -f 1,2,4,5 *.pindel.annot.vcf_pass_flags|sort|uniq>pindel_pileup/$INDEL_BEDFILE_NAME
cd $ROOT_DIR/pindel_raw/pindel_pileup

#Run createVafCmd for each project id containing samples - with input value '1' (this is the option for selecting the pindel files as input)
for PROJECT_NUMBER in "${PROJECT_ARRAY[@]}"; do
    echo "1"|createVafCmd.pl -pid $PROJECT_NUMBER  -o output -g $GENOME_FILE -hdr $HIGH_DEPTH_REGIONS -mq 30  -bo 1  -b $INDEL_BEDFILE_NAME
done

#Then run the "create_split_config_ini.R" script in the output folder
PROJECT_NUMBER=${PROJECT_ARRAY[0]}
cd $ROOT_DIR/pindel_raw/pindel_pileup/output
Rscript $CREATE_SPLIT_CONFIG_SCRIPT -p $PROJECT_NUMBER
cd $ROOT_DIR/pindel_raw/pindel_pileup

#Then re-run the createVafCmd.pl script with the new config file - with input value '1' (this is the option for selecting the pindel files as input)
echo "1"|createVafCmd.pl -pid $PROJECT_NUMBER  -o output -i output/${PROJECT_NUMBER}_cgpVafConfig_split.ini -g $GENOME_FILE -hdr $HIGH_DEPTH_REGIONS -mq 30  -bo 1  -b $INDEL_BEDFILE_NAME

#Update the run_bsub.sh command to allow more jobs in the array to run together & get more memory
sed -e 's/\%5/\%50/g;s/2000/4000/g' run_bsub.sh >run_bsub_updated.sh

#Run this if need to switch to the long queue (not normally necessary)
#sed -i -e 's/normal/long/g' run_bsub_updated.sh

bash run_bsub_updated.sh

#--------------------ONCE cgpVAF HAS COMPLETED-----------------------------
#-----------------------------SNV merge-----------------------------
# bsub -o $PWD/log.%J \
#     -e $PWD/err.%J \
#     -q normal \
#     -G team78-grp \
#     -R 'select[mem>=8000] span[hosts=1] rusage[mem=8000]' \
#     -M8000 \
#     -n1 \
#     -J "SNV_merge" \
#     /lustre/scratch119/casm/team154pc/ms56/my_programs/INDEL_merge.sh $EXP_ID $ROOT_DIR

cd $ROOT_DIR/caveman_raw/caveman_pileup/output/output/$INSILICO_NORMAL_ID/snp
ls *_vaf.tsv > files

#for first file
cut -f 3,4,5,6,24,26,39,41,54,56,69,71,84,86,99,101,114,116,129,131,144,146,159,161,174,176 $(sed -n '1p' files) > temp.1.cut   #Chr, Pos, Ref, Alt, + MTR, DEP for all samples (max 11 samples in one file)

#for subsequent files (exclude Chr, Pos, Ref, Alt and PDv37is)
for FILE in $(tail -n+2 files); do    
    if [ -s temp.$FILE ]
    then
        echo "temp file temp.$FILE already exists. Moving onto next file..."
    else
        echo "temp file temp.$FILE does not yet exist, will be created"
        cut -f 39,41,54,56,69,71,84,86,99,101,114,116,129,131,144,146,159,161,174,176 $FILE >  temp.$FILE.cut
    fi       
done

#Remove empty rows where header was with awk
for FILE in $(ls temp.*); do
    echo $FILE
    awk 'NF' $FILE > output.$FILE
    rm $FILE
done

#Concatenate output files to one merged file & move to the root directory
paste output.* > merged_SNVs_${EXP_ID}.tsv && rm output.*

mv $ROOT_DIR/caveman_raw/caveman_pileup/output/output/$INSILICO_NORMAL_ID/snp/merged_SNVs_${EXP_ID}.tsv $ROOT_DIR/

#-----------------------------INDEL merge-----------------------------
cd $ROOT_DIR/pindel_raw/pindel_pileup/output/output/$INSILICO_NORMAL_ID/indel
ls *_vaf.tsv > files

#for first file
cut -f 3,4,5,6,16,18,25,27,34,36,43,45,52,54,61,63,70,72,79,81,88,90,97,99,106,108,115,117 $(sed -n '1p' files) > temp.1   #Chr, Pos, Ref, Alt, + MTR, DEP for all samples (max 11 samples in one file)

#for subsequent files (exclude Chr, Pos, Ref, Alt and PDv37is)
for FILE in $(tail -n+2 files); do
    cut -f 25,27,34,36,43,45,52,54,61,63,70,72,79,81,88,90,97,99,106,108,115,117 $FILE >  temp.$FILE
done

#Remove empty rows where header was with awk
for FILE in $(ls temp.*); do
    awk 'NF' $FILE > output.$FILE
    rm $FILE
done

#Concatenate output files to one merged file & move to the root directory
paste output.* > merged_indels_${EXP_ID}.tsv && rm output.*
mv $ROOT_DIR/pindel_raw/pindel_pileup/output/output/$INSILICO_NORMAL_ID/indel/merged_indels_${EXP_ID}.tsv $ROOT_DIR/

#--------------------AFTER ALL CGPVAF MATRICES ARE GENERATED-----------------------------
cd $ROOT_DIR
mkdir -p log_files
mkdir -p err_files

#Run the "Mutation_filtering_get_params" script. Run twice: (1) not excluding samples, (2) excluding colonies with peak VAF <0.4
bsub -o $ROOT_DIR/log_files/mats_and_params.log.%J -e $ROOT_DIR/err_files/mats_and_params.err.%J \
    -q normal -G team78-grp -R 'select[mem>=16000] span[hosts=1] rusage[mem=16000]' \
    -M16000 -n6 -J GET_PARAMS \
    Rscript $MUTATION_PARAMETER_SCRIPT \
    -r $RUN_ID \
    -u MDGRCm38is \
    -s $ROOT_DIR/merged_SNVs_${EXP_ID}.tsv \
    -i $ROOT_DIR/merged_indels_${EXP_ID}.tsv \
    -o $MATS_AND_PARAMS_DIR
bsub -o $ROOT_DIR/log_files/mats_and_params.log.%J -e $ROOT_DIR/err_files/mats_and_params.err.%J \
    -q normal -G team78-grp -R 'select[mem>=16000] span[hosts=1] rusage[mem=16000]' \
    -M16000 -n6 -J GET_PARAMS \
    Rscript $MUTATION_PARAMETER_SCRIPT \
    -r $RUN_ID_M \
    -u MDGRCm38is \
    -s $ROOT_DIR/merged_SNVs_${EXP_ID}.tsv \
    -i $ROOT_DIR/merged_indels_${EXP_ID}.tsv \
    -o $MATS_AND_PARAMS_DIR \
    -m \
    -v 0.4

#ONCE MS FILTERS JOBS HAVE FINISHED
cd $ROOT_DIR/MS_filters/output_files
perl /lustre/scratch126/casm/team154pc/ms56/my_programs/filter_for_bed_file.pl
cut -f 1,2,4,5 *_passed.vcf_for_bed_file|sort|uniq>$ROOT_DIR/MS_filters/$MS_FILTERED_BEDFILE_NAME

#ONCE Mutation_filtering_get_paramaters.R SCRIPT HAS COMPLETED
bsub -o $ROOT_DIR/log_files/MSReduce.log.%J -e $ROOT_DIR/err_files/MSReduce.err.%J -q yesterday -R 'select[mem>=16000] span[hosts=1] rusage[mem=16000]' -M16000 -n1 -J MSreduce Rscript /lustre/scratch126/casm/team154pc/ms56/my_programs/Reducing_mutset_from_MSfilters.R -r $RUN_ID -b $ROOT_DIR/MS_filters/$MS_FILTERED_BEDFILE_NAME -d $MATS_AND_PARAMS_DIR
bsub -o $ROOT_DIR/log_files/MSReduce.log.%J -e $ROOT_DIR/err_files/MSReduce.err.%J -q yesterday -R 'select[mem>=16000] span[hosts=1] rusage[mem=16000]' -M16000 -n1 -J MSreduce Rscript /lustre/scratch126/casm/team154pc/ms56/my_programs/Reducing_mutset_from_MSfilters.R -r $RUN_ID_M -b $ROOT_DIR/MS_filters/$MS_FILTERED_BEDFILE_NAME -d $MATS_AND_PARAMS_DIR

#Run the sensitivity analysis
bsub -o $ROOT_DIR/log_files/sensitivity.log.%J -e $ROOT_DIR/err_files/sensitivity.err.%J \
    -q normal -R 'select[mem>=4000] span[hosts=1] rusage[mem=4000]' \
    -M4000 -n1 -J sensitivity \
    Rscript /lustre/scratch126/casm/team154pc/ms56/my_programs/Sensitivity_analysis_from_SNPs.R \
    -m $MATS_AND_PARAMS_DIR/mats_and_params_${RUN_ID}_postMS \
    -o $ROOT_DIR -n sensitivity_analysis_${EXP_ID} \
    -i $ROOT_DIR/pindel_raw \
    -s $ROOT_DIR/MS_filters/output_files \
    -x '_passed.vcf_for_bed_file'

#-----------------------Run the tree-building script - this has lots of options

# -i The id for the tree-building run - will be included in the output file names.
# -m Path to mutation filtering output file
# -f Option to do p-value based filtering, or vaf based filtering (pval or vaf)
# -c Cut off to exclude samples with low coverage
# -o Folder for storing script output files. Folder & subfolders will be created if don't already exist.
# -s Path to the sensitivity dataframe
# -p Save trees as polytomous trees (not bifurcating trees)
# -t The age ('time') of the individual - for building the age-adjusted ultrametric tree
# -j Option to do initial tree-building with just the SNVs (i.e. don't' include indels)
# -a Option to keep an ancestral branch
# -q an optional text file of genes of interest to annotate on the tree

bsub -o $ROOT_DIR/log_files/treebuild.log.%J -e $ROOT_DIR/err_files/treebuild.err.%J \
    -q normal -R 'select[mem>=4000] span[hosts=1] rusage[mem=4000]' \
    -M4000 -n1 -J tree_build \
    Rscript $TREE_BUILDING_SCRIPT \
    -i ${RUN_ID_TB}_a_j \
    -m ${STUDY_DIR}/filtering_runs/mats_and_params/mats_and_params_${RUN_ID_TB} \
    -f pval \
    -d $EXP_ID \
    -c 4 \
    -o ${STUDY_DIR}/filtering_runs \
    -s $ROOT_DIR/sensitivity_analysis_${EXP_ID} \
    -p \
    -t $DONOR_AGE \
    -j \
    -q $GENES_TO_ANNOTATE \
    -a 

bsub -o $ROOT_DIR/log_files/treebuild.log.%J -e $ROOT_DIR/err_files/treebuild.err.%J \
    -q normal -R 'select[mem>=4000] span[hosts=1] rusage[mem=4000]' \
    -M4000 -n1 -J tree_build \
    Rscript $TREE_BUILDING_SCRIPT \
    -i ${RUN_ID_TB}_a_j \
    -m ${STUDY_DIR}/filtering_runs/mats_and_params/mats_and_params_${RUN_ID_TB} \
    -f vaf \
    -d $EXP_ID \
    -c 4 \
    -o ${STUDY_DIR}/filtering_runs \
    -s $ROOT_DIR/sensitivity_analysis_${EXP_ID} \
    -p \
    -t $DONOR_AGE \
    -j \
    -q $GENES_TO_ANNOTATE \
    -a


##############



SAMPLES=$(ls|cut -d '.' -f1|grep 'MD'|uniq)
for SAMPLE_ID in $SAMPLES; do
    bsub -o $PWD/log_files/${SAMPLE_ID}.log.%J \
    -e $PWD/log_files/${SAMPLE_ID}.log.%J \
    -q normal \
    -G team78-grp \
    -R 'select[mem>=16000] span[hosts=1] rusage[mem=16000]' \
    -M16000 \
    -n1 \
    -J "MS_filt" \
    /lustre/scratch126/casm/team154pc/ms56/my_programs/LCMfiltering_2024/runlcmfilter.sh /nfs/cancer_ref01/nst_links/live/3357/${SAMPLE_ID}/${SAMPLE_ID}.caveman_c.annot.vcf.gz /nfs/cancer_ref01/nst_links/live/3357/${SAMPLE_ID}/${SAMPLE_ID}.sample.dupmarked.bam output_files/ ${SAMPLE_ID} mm10 3
done


