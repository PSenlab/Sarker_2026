#!/usr/bin/env python3
# ==============================================================================
# SCENIC+ Step 2: MALLET LDA Model Fitting on cisTopic Object
# ==============================================================================
#
# Description:
#   Fits LDA (Latent Dirichlet Allocation) topic models to the filtered
#   cisTopic object using MALLET. Multiple topic counts are tested; the
#   best model is selected downstream in Step 3 using model-quality metrics.
#
# Input:
#   - cisTopicObject_filtered_annotated.pkl
#     (output of Step 1: step1_pycistopic_preprocessing.py)
#
# Output:
#   - outs/mal_result/mallet_model_<N>.pkl  (one per topic count)
#
# Prerequisites:
#   - MALLET (http://mallet.cs.umass.edu/) installed and accessible
#   - Java 8 available on PATH (JAVA_HOME set appropriately)
#   - Sufficient memory (~500G) for large cisTopic objects
#
#
# ==============================================================================

import os
import pickle
from pycisTopic.lda_models import run_cgs_models_mallet


# ==============================================================================
# CONFIGURATION - UPDATE THESE PATHS
# ==============================================================================

OUT_DIR      = "outs/"
TMP_PATH     = os.path.join(OUT_DIR, "mal_temp")
RESULT_PATH  = os.path.join(OUT_DIR, "mal_result")

# Path to the MALLET binary
# Download from http://mallet.cs.umass.edu/download.php and unpack locally
MALLET_BIN   = "path/to/Mallet-202108/bin/mallet"

# Topic counts to fit (one model per value)
N_TOPICS     = [25, 50, 75, 100, 300, 350, 400]
N_CPU        = 20
N_ITER       = 500
RANDOM_STATE = 555
ALPHA        = 50
ETA          = 0.1


# ==============================================================================
# SETUP
# ==============================================================================

os.makedirs(TMP_PATH, exist_ok=True)
os.makedirs(RESULT_PATH, exist_ok=True)

# MALLET / JVM memory - scale to your machine
os.environ["JAVA_TOOL_OPTIONS"] = "-Xms500g -Xmx500g -XX:+UseG1GC"
os.environ["MALLET_MEMORY"]     = "500G"


# ==============================================================================
# LOAD cisTopic OBJECT
# ==============================================================================

cistopic_pkl = os.path.join(OUT_DIR, "cisTopicObject_filtered_annotated.pkl")
print(f"Loading cisTopic object: {cistopic_pkl}")
with open(cistopic_pkl, "rb") as infile:
    cistopic_obj = pickle.load(infile)
print("  [OK] Loaded cisTopic object")


# ==============================================================================
# FIT LDA MODELS
# ==============================================================================

print(f"\nFitting MALLET LDA models for n_topics = {N_TOPICS}")
print(f"  CPUs:  {N_CPU}")
print(f"  Iters: {N_ITER}")
print(f"  Alpha: {ALPHA} (by_topic=True)")
print(f"  Eta:   {ETA} (by_topic=False)")
print(f"  Seed:  {RANDOM_STATE}")

try:
    models = run_cgs_models_mallet(
        mallet_path    = MALLET_BIN,
        cistopic_obj   = cistopic_obj,
        n_topics       = N_TOPICS,
        n_cpu          = N_CPU,
        n_iter         = N_ITER,
        random_state   = RANDOM_STATE,
        alpha          = ALPHA,
        alpha_by_topic = True,
        eta            = ETA,
        eta_by_topic   = False,
        tmp_path       = TMP_PATH,
        save_path      = RESULT_PATH,
    )
except Exception as e:
    print(f"[ERROR] MALLET run failed: {e}")
    models = None


# ==============================================================================
# SAVE INDIVIDUAL MODELS
# ==============================================================================

if models is not None:
    for model in models:
        topic_number = model.topic_no
        model_path = os.path.join(RESULT_PATH, f"mallet_model_{topic_number}.pkl")
        with open(model_path, "wb") as f:
            pickle.dump(model, f)
        print(f"  [OK] {model_path}")
    print("\nAll models saved successfully.")
else:
    print("\n[ERROR] No models to save.")
