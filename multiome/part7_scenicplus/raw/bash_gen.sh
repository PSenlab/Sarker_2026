#!/bin/bash
ulimit -u 4096  # Increase the number of user processes
ulimit -n 8192  # Increase the number of open files



source myconda
mamba activate scenicplus5

#python /data/sarkern2/scenicplus/pycistopic_script.py
#python /data/sarkern2/scenicplus/pycistopic_script_mallet.py
python 03_cistopic_merge_dbl.py


# sbatch --cpus-per-task=30 --mem=280g --time=10:00:00 solTE_bash.sh # example of command line code
# sbatch --cpus-per-task=30 --mem=280g --time=10:00:00 solTE_bash.sh # example of command line code
# sbatch --partition=gpu --gres=gpu:v100x:2 --cpus-per-task=50 --mem=600g --time=10:00:00 solTE_bash.sh
#sbatch --partition=gpu --cpus-per-task=50 --gres=gpu:k80:4 bash_gen.sh