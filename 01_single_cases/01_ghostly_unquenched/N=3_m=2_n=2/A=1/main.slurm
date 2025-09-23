#!/bin/bash
#SBATCH --job-name=singl					# create a short name for your job
#SBATCH --nodes=1									# node count
#SBATCH --ntasks=1								# total number of tasks across all nodes
#SBATCH --cpus-per-task=8					# cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem=256G								# memory in total
#SBATCH --time=96:00:00						# total run time limit (HH:MM:SS)
#SBATCH --mail-type=begin					# send email when job begins
#SBATCH --mail-type=end						# send email when job ends
#SBATCH --mail-user=aaron.held@phys.ens.fr


# run main jobs

julia --threads 8 main.jl
