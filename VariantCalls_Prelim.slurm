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
STUDY='VariantCall_Prelim'
project=Prelim_NF
#Path to Key file and filename
KEY='06_Key/2021Prelim_KeyFile.txt'

#Restriction enzyme used for GBS
E='PstI-MspI'

#Folder with the FASTQ files - Note: DO NOT INLCLUDE FILENAME
FASTQ_Folder='04_fastq'

#REFERENCE GENOME
RG='02_ReferenceGenome/iwgsc_refseqv1.0_all_chromosomes/161010_Chinese_Spring_v1.0_pseudomolecules.fasta'


#Create directory
mkdir ./database ./logs ./tagsForAlign ./tagsAligned ./SNPQualitySummary ./ProductionSNPCallerPluginV2

#Load Java required for Tassel and get the path and add program name (run_pipeline.pl)
#Note the memory provided to TASSEL here should match with the information in the SLURM header
TASSEL='05_Software/tassel-5-standalone/run_pipeline.pl -Xms32g -Xmx32g'

##Step1: GBSSeqToTagDBPlugin
$TASSEL -fork1 -GBSSeqToTagDBPlugin -deleteOldData true -e $E -i $FASTQ_Folder -db database/$STUDY.db -k $KEY -kmerLength 64 -minKmerL 20 -mnQS 20 -mxKmerNum 100000000 -endPlugin -runfork1 > ./logs/GBSSeqToTagDBPlugin_NF.log
#Default quality score changed from 0 to 20

##Step2: TagExportToFastqPlugin
$TASSEL -fork1 -TagExportToFastqPlugin -db database/$STUDY.db -o tagsForAlign/$STUDY\_tagsForAlign.fa.gz -c 5 -endPlugin -runfork1 > ./logs/TagExportToFastqPlugin_NF.log
#Minimum count of reads across taxa for a tag to be output (Default: 1) changed to 5

#Path to BWA plus the name of the program bwa
BWA='05_Software/bwa/bwa'

##Step3: Alignment of tags to Reference Genome
##Note we have mentioned 8 threads are available in the SLURM header (ntasks-per-node=8)
$BWA aln -t 8 $RG tagsForAlign/$STUDY\_tagsForAlign.fa.gz > tagsAligned/tagsAligned.sai
#Using bwa aln invokes BWA-backtrack algorithm. Note for shorter reads this is the preferred algorithm.

##Step4: Convert the *.sai output file from BWA to *.sam format
$BWA samse $RG tagsAligned/tagsAligned.sai tagsForAlign/$STUDY\_tagsForAlign.fa.gz > tagsAligned/tagsAligned.sam

##Step5: SAMToGBSdbPlugin
$TASSEL -fork1 -SAMToGBSdbPlugin -i tagsAligned/tagsAligned.sam -db database/$STUDY.db -endPlugin -runfork1 > ./logs/SAMToGBSdbPlugin_NF.log

##Step6: DiscoverySNPCallerPluginV2
$TASSEL -fork1 -DiscoverySNPCallerPluginV2 -db database/$STUDY.db -sC "chr1A" -eC "chrUn" -mnLCov 0.01 -mnMAF 0.01 -deleteOldData true -endPlugin -runfork1 > ./logs/DiscoverySNPCallerPluginV2_NF.log

##Step7: SNPQualityProfilerPlugin
$TASSEL -fork1 -SNPQualityProfilerPlugin -db database/$STUDY.db -statFile SNPQualitySummary/outputStats.txt -deleteOldData true -endPlugin -runfork1 > ./logs/SNPQualityProfilerPlugin_NF.log

##Step8: ProductionSNPCallerPluginV2
$TASSEL -fork1 -ProductionSNPCallerPluginV2 -ko false -db database/$STUDY.db -e $E -i $FASTQ_Folder -k $KEY -kmerLength 64 -o ProductionSNPCallerPluginV2/${project}.h5 -endPlugin -runfork1 > ./logs/ProductionSNPCallerPluginV2_NF.log

##Step9: Convert the h5 file into vcf and hapmap formats
$TASSEL -h5 ProductionSNPCallerPluginV2/DiscoveryBuild_NF.h5 -export ProductionSNPCallerPluginV2/${project} -exportType Hapmap
$TASSEL -h5 ProductionSNPCallerPluginV2/DiscoveryBuild_NF.h5 -export ProductionSNPCallerPluginV2/${project} -exportType VCF
