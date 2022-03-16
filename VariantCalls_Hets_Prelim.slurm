#!/bin/sh
#SBATCH --partition=kamiak,cahnrs,cahnrs_bigmem	# Partition/Queue to use
#SBATCH --job-name=VC	# Job name
#SBATCH --output=VC_%j.out	# Output file (stdout)
#SBATCH --error=VC_%j.err	# Error file (stderr)
#SBATCH --time=24:00:00	# Wall clock time limit Days-HH:MM:SS
#SBATCH --mail-type=ALL	# Email notification: BEGIN,END,FAIL,ALL
#SBATCH --mail-user=lance.merrick@wsu.edu	# Email address for notifications
#SBATCH --nodes=1		# Number of nodes (min-max)
#SBATCH --ntasks-per-node=8	# Number of tasks per node (max)
#SBATCH --ntasks=1		# Number of tasks (processes)
#SBATCH --cpus-per-task=1	# Number of cores per task (threads)
#SBATCH --mem-per-cpu=32G	# Memory per core (gigabytes)

module load fastqc
module load java
module load vcftools
module load perl
module load sqlite3
module load python
module load bowtie2
module load trimmomatic
module load bwa
module load samtools
module load picard
module load gatk
module load beagle
module load plink

#Name for the project to store all the information in a single folder
STUDY='VariantCall_Prelim_Het'
#Path to work directory - VariantCalls
#WD=~/Scripts/$STUDY
#cd $WD

#Name of output
project=Prelim_Het
#Path to Key file and filename
KEY='06_Key/2021Prelim_KeyFile.txt'

#Restriction enzyme used for GBS
E='PstI-MspI'

#Folder with the FASTQ files - Note: DO NOT INLCLUDE FILENAME
FASTQ_Folder='04_fastq'

#REFERENCE GENOME
RG='02_ReferenceGenome/iwgsc_refseqv1.0_all_chromosomes/161010_Chinese_Spring_v1.0_pseudomolecules.fasta'



#Create directory within STUDY folder
mkdir ./database ./logs ./tagsForAlign ./tagsAligned ./SNPQualitySummary ./ProductionSNPCallerPluginV2


#Load Java required for Tassel and get the path and add program name (run_pipeline.pl)
#Note the memory provided to TASSEL here should match with the information in the SLURM header
TASSEL='05_Software/tassel-5-standalone/run_pipeline.pl -Xms32g -Xmx32g'


##Step1: GBSSeqToTagDBPlugin
$TASSEL -fork1 -GBSSeqToTagDBPlugin -deleteOldData true -e $E -i $FASTQ_Folder -db database/$STUDY.db -k $KEY -kmerLength 64 -minKmerL 20 -mnQS 20 -mxKmerNum 100000000 -endPlugin -runfork1 > ./logs/GBSSeqToTagDBPlugin_2014.log
#Default quality score changed from 0 to 20

##Step2: TagExportToFastqPlugin
$TASSEL -fork1 -TagExportToFastqPlugin -db database/$STUDY.db -o tagsForAlign/$STUDY\_tagsForAlign.fa.gz -c 5 -endPlugin -runfork1 > ./logs/TagExportToFastqPlugin.log
#Minimum count of reads across taxa for a tag to be output (Default: 1) changed to 5

#Path to BWA plus the name of the program bwa
BWA='05_Software/bwa/bwa'

##Step3: Alignment of tags to Reference Genome
##Note we have mentioned 8 threads are available in the SLURM header (ntasks-per-node=8)
$BWA aln -t 8 $RG tagsForAlign/$STUDY\_tagsForAlign.fa.gz > tagsAligned/tagsAligned.sai
#Using bwa aln invokes BWA-backtrack algorithm. Note for shorter reads this is the preferred algorithm.

##Step4: Convert the *.sai output file from BWA to *.sam format
$BWA samse $RG tagsAligned/tagsAligned.sai tagsForAlign/$STUDY\_tagsForAlign.fa.gz > tagsAligned/tagsAligned.sam

#3.5 only get unique snps
grep -v "XS:i" tagsAligned/tagsAligned.sam | awk '$4!=0' > tagsAligned/tagsAligned_unique.sam

##Step5: SAMToGBSdbPlugin
$TASSEL -fork1 -SAMToGBSdbPlugin -i tagsAligned/tagsAligned_unique.sam -db database/$STUDY.db -endPlugin -runfork1 > ./logs/SAMToGBSdbPlugin.log

##Step6: DiscoverySNPCallerPluginV2
$TASSEL -fork1 -DiscoverySNPCallerPluginV2 -db database/$STUDY.db -sC "chr1A" -eC "chrUn" -mnLCov 0.01 -mnMAF 0.01 -deleteOldData true -endPlugin -runfork1 > ./logs/DiscoverySNPCallerPluginV2.log

##Step7: SNPQualityProfilerPlugin
$TASSEL -fork1 -SNPQualityProfilerPlugin -db database/$STUDY.db -statFile SNPQualitySummary/outputStats.txt -deleteOldData true -endPlugin -runfork1 > ./logs/SNPQualityProfilerPlugin.log

#Step 7.5
## UpdateSNPPositionQualityPlugin - UPDATE DATABASE WITH QUALITY SCORE fast < 30 minutes 15GB
$TASSEL -fork1 -UpdateSNPPositionQualityPlugin -db database/$STUDY.db -qsFile SNPQualitySummary/SNPqual_stats.txt -endPlugin -runfork1 > ./logs/UpdateSNPPositionQualityPlugin.log

##Step8: ProductionSNPCallerPluginV2
$TASSEL -fork1 -ProductionSNPCallerPluginV2 -ko false -db database/$STUDY.db -e $E -i $FASTQ_Folder -k $KEY -kmerLength 64 -o ProductionSNPCallerPluginV2/${project}.h5 -endPlugin -runfork1 > ./logs/ProductionSNPCallerPluginV2_2014_run1.h5.log

##Step9: Convert the h5 file into vcf and hapmap formats
#Will be used as base comparisons but kept in a different folder
$TASSEL -h5 ProductionSNPCallerPluginV2/${project}.h5 -export ProductionSNPCallerPluginV2/${project} -exportType Hapmap
$TASSEL -h5 ProductionSNPCallerPluginV2/${project}.h5 -export ProductionSNPCallerPluginV2/${project} -exportType VCF

##Step9.5: ProductionSNPCallerPluginV2
$TASSEL -fork1 -ProductionSNPCallerPluginV2 -ko false -db database/$STUDY.db -e $E -i $FASTQ_Folder -k $KEY -kmerLength 64 -o ${project}.vcf -endPlugin -runfork1 > ./logs/ProductionSNPCallerPluginV2_2014_run1_V2.vcf.log

#Step 10 Extract tag sequences using code from GetTagSequenceFromDBPlugin
$TASSEL -fork1 -GetTagSequenceFromDBPlugin \
    -db database/$STUDY.db \
    -o ${project}.tags.txt \
    -endPlugin -runfork1 > ./logs/GetTagSequenceFromDBPlugin.log

##Step 10 and 11 make output that is probably not used, but included anyway
##Step 10 Use code from Liang Gao to get Tag position and merge with tag sequence
## extract tag sequences in the db (sqlite3)
## extract tagid, snpid, table
## Note the location of the script GBSv2_table_join_to_get_tagid.sql might be different for you if you do not use KSU computer cluster

#Scripts for step 10 and 11
JOIN='GBSv2_table_join_to_get_tagid.sql'
ID='snp_to_tagid.py'


sqlite3 -separator $'\t' database/$STUDY.db < \
$JOIN \
        > ${project}.joined.table.txt

##Step 11 Code from Liang Gao
## using the python script to link ID, sequences together.
## Note the location of the script snp_to_tagid.py might be different for you if you donot use KSU computer cluster
python $ID \
        -t ${project}.tags.txt \
        -s  ${project}.joined.table.txt \
        > ${project}.snpid.linked.seq.fa
#    /hom
#Filter out for Cycle 6 data only
#vcftools --vcf ${name}.vcf --keep ../../data/Original_Data/Cycle6_VCF_Filter.txt --out ${name} --recode

#Rename file
#mv ${name}.recode.vcf ${name}.vcf #obliterates original vcf file, but for this if fine

#Get allele depth counts
vcftools --vcf ${project}.vcf --out ${project} --geno-depth

#Python script for filter based on a depth of 4
GDEPTH='Genotype_by_Depth4.pl'
#Filter based on perl
perl $GDEPTH ${project}.vcf > ${project}_Diff_Call.vcf

#Everything below is just if you want to filter out the vcf
#Filter on minimum of two reads and 70% max missing, MAF > 0.01, no indels or multiallelic calls.
vcftools --vcf ${project}_Diff_Call.vcf  --out ${project}_Diff_Filtered  --remove-indels --min-alleles 2 --max-alleles 2 --maf 0.01 --minDP 2 --max-missing 0.3 --recode --recode-INFO-all
#vcftools --vcf ${name}.vcf  --out ${name}_Diff_Filtered  --remove-indels --min-alleles 2 --max-alleles 2 --maf 0.01 --minDP 2 --max-missing 0.3 --recode --recode-INFO-all
#rename file
mv ${project}_Diff_Filtered.recode.vcf ${project}_Diff_Filtered.vcf

#get final depth
vcftools --vcf ${project}_Diff_Filtered.vcf --out ${project}_Diff_Filtered --geno-depth

## Convert to Hapmap format
$TASSEL -fork1 -vcf ${project}_Diff_Filtered.vcf -export ${project}_Diff_Filtered -exportType Hapmap -runfork1 >>  DiscoveryBuild_2014_run1_V2_pipeline.out


###SNP imputation, proceed to GWAS from this step unless calculating LD
beagle gt=${project}_Diff_Filtered.vcf out=${project}_imputation.vcf #map=
#Decompress output
gzip -d ${project}_imputation.vcf.gz

#Below code is not necessary unless you want to look and filter for LD and output for STRUCTURE
###convert to PLINK
vcftools --vcf ${project}_imputation.vcf --plink --out ${project}

###calculate LD for figure, set --ld-window-kb to length of longest chr
plink --file ${project} --r2 --ld-window-r2 0 --ld-window 9999 --ld-window-kb 5292 --out ${project} ###is this actually any good, stolen from forum
plink --file ${project} --r2 --ld-window-r2 0 --ld-window 9999 --ld-window-kb 529 --out ${project}2
plink --file ${project} --r2 --ld-window-r2 0 --out ${project}2

###remove markers in LD
plink --file ${project} --allow-extra-chr --make-founders --indep-pairphase 50 5 0.9 --out LD_check

###thin markers for STRUCTURE
vcftools --vcf ${project}_imputation.vcf --thin 50000 --recode --recode-INFO-all --out ${project}.STRUCTURESNPs

#Clean up compress objects
#zip original files
gzip tagsAligned/tagsAligned.sam
#zip unique sam file
gzip tagsAligned/tagsAligned_unique.sam
#zip full vcf file
gzip ${project}.vcf
#zip depth call
gzip ${project}_Diff_Call.vcf
#zipfull depth call
gzip ${project}.gdepth
#zip filtered final vcf
gzip  ${project}_Diff_Filtered.vcf
#zip imputed filtered final vcf
gzip  ${project}_imputation.vcf.gz
