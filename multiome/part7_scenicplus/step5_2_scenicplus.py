import os
import scenicplus
import scanpy as sc
import subprocess


os.system("snakemake --unlock")


import resource
from joblib import Parallel, delayed




# Initialize Snakemake for SCENIC+
#os.system("mkdir -p scplus_pipeline")
#os.system("scenicplus init_snakemake --out_dir scplus_pipeline")

# Display the directory structure to verify initialization
#os.system("tree scplus_pipeline/")
#os.system("tree outs/scenicplus/region_sets")

# Display the content of the configuration file to ensure it is set correctly
#os.system("cat scplus_pipeline/Snakemake/config/config.yaml")

# Change directory to the Snakemake folder
#os.chdir("scplus_pipeline/Snakemake/")

# Run Snakemake with 10 cores (adjust as needed)
os.system("snakemake --cores 50")
