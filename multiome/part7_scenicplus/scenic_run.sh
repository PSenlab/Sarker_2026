#!/bin/bash
#SBATCH --job-name=scenicplus
#SBATCH --output=scenicplus_output.log
#SBATCH --error=scenicplus_error.log
#SBATCH --time=148:00:00
#SBATCH --partition=largemem
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=50
#SBATCH --mem=1600G



ulimit -u 4096  # Increase the number of user processes
ulimit -n 8192  # Increase the number of open files

# Optional: Set the maximum file size (in blocks, where 1 block = 512 bytes)
# ulimit -f 10000

source myconda
mamba activate scenicplus5


python scenicplus.py
