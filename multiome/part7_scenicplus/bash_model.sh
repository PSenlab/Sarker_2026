#!/bin/bash
ulimit -u 4096  # Increase the number of user processes
ulimit -n 8192  # Increase the number of open files

# Optional: Set the maximum file size (in blocks, where 1 block = 512 bytes)
# ulimit -f 10000

source myconda
mamba activate scenicplus5


# Set JAVA_HOME based on your discovered Java directory
export JAVA_HOME=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.412.b08-2.el8.x86_64
export PATH=$JAVA_HOME/bin:$PATH

# Set MALLET_HOME to the MALLET directory
export MALLET_HOME=$PWD/Mallet-202108

# Add MALLET_HOME to PATH
export PATH=$MALLET_HOME/bin:$PATH

# Set MALLET_MEMORY for MALLET
export MALLET_MEMORY=500G

# Set Java tool options for heap size
export JAVA_TOOL_OPTIONS="-Xms500g -Xmx500g -XX:+UseG1GC"


# Print to verify
echo "JAVA_HOME: $JAVA_HOME"
echo "JAVA_TOOL_OPTIONS: $JAVA_TOOL_OPTIONS"





#python /data/sarkern2/scenicplus/pycistopic_script.py
#python /data/sarkern2/scenicplus/pycistopic_script_mallet.py
python mallet_model.py

# sbatch --cpus-per-task=30 --mem=280g --time=10:00:00 solTE_bash.sh # example of command line code
# sbatch --cpus-per-task=30 --mem=280g --time=10:00:00 solTE_bash.sh # example of command line code
# sbatch --partition=gpu --gres=gpu:v100x:2 --cpus-per-task=40 --mem=280g --time=10:00:00 solTE_bash.sh