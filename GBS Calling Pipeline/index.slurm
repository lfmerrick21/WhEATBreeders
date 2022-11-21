#!/bin/sh
#SBATCH --partition=cahnrs_bigmem,kamiak,cahnrs	# Partition/Queue to use
#SBATCH --job-name=IA	# Job name
#SBATCH --output=output/IA_%j.out	# Output file (stdout)
#SBATCH --error=output/IA_%j.err	# Error file (stderr)
#SBATCH --time=7-00:00:00	# Wall clock time limit Days-HH:MM:SS
#SBATCH --mail-type=ALL	# Email notification: BEGIN,END,FAIL,ALL
#SBATCH --mail-user=lance.merrick@wsu.edu	# Email address for notifications
#SBATCH --nodes=1		# Number of nodes (min-max)
#SBATCH --ntasks-per-node=1	# Number of tasks per node (max)
#SBATCH --ntasks=1		# Number of tasks (processes)
#SBATCH --cpus-per-task=1	# Number of cores per task (threads)
#SBATCH --mem-per-cpu=250G	# Memory per core (gigabytes)

#Note path to BWA program
BWA='05_Software/bwa/bwa'

#Define the reference genome
RG='02_ReferenceGenome/iwgsc_refseqv1.0_all_chromosomes/161010_Chinese_Spring_v1.0_pseudomolecules.fasta'

$BWA index -a bwtsw $RG
