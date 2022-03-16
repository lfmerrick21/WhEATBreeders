#!/bin/sh
#SBATCH --partition=cahnrs_bigmem,kamiak,cahnrs	# Partition/Queue to use
#SBATCH --job-name=fastqc	# Job name
#SBATCH --output=output/fastqc_%j.out	# Output file (stdout)
#SBATCH --error=output/fastqc_%j.err	# Error file (stderr)
#SBATCH --time=7-00:00:00	# Wall clock time limit Days-HH:MM:SS
#SBATCH --mail-type=ALL	# Email notification: BEGIN,END,FAIL,ALL
#SBATCH --mail-user=lance.merrick@wsu.edu	# Email address for notifications
#SBATCH --nodes=1		# Number of nodes (min-max)
#SBATCH --ntasks-per-node=1	# Number of tasks per node (max)
#SBATCH --ntasks=1		# Number of tasks (processes)
#SBATCH --cpus-per-task=1	# Number of cores per task (threads)
##SBATCH --mem-per-cpu=112G	# Memory per core (gigabytes)

#module load fastqc/0.11
module load fastqc
module load java
mkdir fastqc
#fastqc -o fastqc /work/demo/SHARED/04_fastq/GBS0551xDUP1_C4RF0ACXX_s_1_fastq.txt.gz
fastqc -o fastqc 04_fastq/HKFMKDRXY_2_fastq.gz
