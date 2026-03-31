import os
import pandas as pd
import numpy as np
import pickle
from pycisTopic.lda_models import *

# Define paths
outDir = "outs_trial/"
tmp_path = os.path.join(outDir, "mal_temp")
mal_result_path = os.path.join(outDir, "mal_result")

path_to_mallet_binary = "/data/sarkern2/scenicplus1/Mallet-202108/bin/mallet"

# Ensure output directories exist
os.makedirs(tmp_path, exist_ok=True)
os.makedirs(mal_result_path, exist_ok=True)


# Load cisTopic object
with open(os.path.join(outDir, 'cisTopicObject_filtered_annotated.pkl'), 'rb') as infile:
    cistopic_obj = pickle.load(infile)

# Configure MALLET environment
os.environ["JAVA_TOOL_OPTIONS"] = "-Xms500g -Xmx500g -XX:+UseG1GC"
os.environ['MALLET_MEMORY'] = '500G'  # Allocate 600G memory for MALLET

# Run models
try:
    models = run_cgs_models_mallet(
        mallet_path=path_to_mallet_binary,  # specify the mallet path as mallet_path
        cistopic_obj=cistopic_obj,
        n_topics=[300, 350, 400],
        n_cpu=20,
        n_iter=500,
        random_state=555,
        alpha=50,
        alpha_by_topic=True,
        eta=0.1,
        eta_by_topic=False,
        tmp_path=tmp_path,
        save_path=mal_result_path
    )
except Exception as e:
    print(f"Error running models: {e}")
    models = None

# Save individual models if they were successfully created
if models is not None:
    try:
        for model in models:
            topic_number = model.topic_no  # Assuming the model has an attribute `topic_no` for the number of topics
            model_path = os.path.join(models_dir, f'mallet_model_{topic_number}.pkl')
            with open(model_path, 'wb') as f:
                pickle.dump(model, f)
        print("Models saved successfully.")
    except Exception as e:
        print(f"Error saving models: {e}")
else:
    print("No models to save.")
