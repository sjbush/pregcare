#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=guppy
#SBATCH --ntasks=10
#SBATCH --mem=10G
#SBATCH --time=15:00:00
#SBATCH --output=%j_%x.out
#SBATCH --error=%j_%x.err
#SBATCH --mail-user=stephen.bush@ndcls.ox.ac.uk
#SBATCH --mail-type=end,fail

/home/s/sbush/programs/ont-guppy-cpu/bin/guppy_basecaller -i /t1-data/project/pregcare/sbush/MinION8/1.fast5/$1 -s /t1-data/project/pregcare/sbush/MinION8/2.fast5_to_fastq/$1 --num_callers 10 --cpu_threads_per_caller 1 --compress_fastq --config dna_r9.4.1_450bps_hac.cfg
