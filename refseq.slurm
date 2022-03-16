#!/bin/sh
#SBATCH --partition=cahnrs_bigmem,kamiak,cahnrs	# Partition/Queue to use
#SBATCH --job-name=RF	# Job name
#SBATCH --output=output/RF_%j.out	# Output file (stdout)
#SBATCH --error=output/RF_%j.err	# Error file (stderr)
#SBATCH --time=7-00:00:00	# Wall clock time limit Days-HH:MM:SS
#SBATCH --mail-type=ALL	# Email notification: BEGIN,END,FAIL,ALL
#SBATCH --mail-user=lance.merrick@wsu.edu	# Email address for notifications
#SBATCH --nodes=1		# Number of nodes (min-max)
#SBATCH --ntasks-per-node=1	# Number of tasks per node (max)
#SBATCH --ntasks=1		# Number of tasks (processes)
#SBATCH --cpus-per-task=1	# Number of cores per task (threads)
#SBATCH --mem-per-cpu=90G	# Memory per core (gigabytes)


wget -c https://urgi.versailles.inra.fr/download/iwgsc/IWGSC_RefSeq_Assemblies/v1.0/iwgsc_refseqv1.0_all_chromosomes.zip
