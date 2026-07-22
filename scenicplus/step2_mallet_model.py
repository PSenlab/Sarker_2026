#!/usr/bin/env python3
# ==============================================================================
# SCENIC+ Step 2: MALLET LDA Model Fitting — ALL COMPARTMENTS
# ==============================================================================
#
# Description:
#   Fits LDA topic models with MALLET to each compartment's filtered cisTopic
#   object. Multiple topic counts are fitted per compartment; the best model is
#   chosen in Step 3 from the evaluation metrics (Arun / Cao-Juan crossover,
#   Mimno coherence, log-likelihood).
#
# Input  (per compartment, from Step 1):
#   <out_dir>/cisTopicObject_filtered_annotated_<suffix>.pkl
#
# Output (per compartment):
#   <out_dir>/mal_result/mallet_model_<N>.pkl        (one per topic count)
#
# Prerequisites:
#   - MALLET (http://mallet.cs.umass.edu/) unpacked locally
#   - Java 8 on PATH (JAVA_HOME set appropriately)
#   - Memory scaled to the object: see JVM_MEMORY / per-compartment override
#
# Usage:
#   python step2_mallet_lda_all.py                      # all compartments
#   python step2_mallet_lda_all.py --run T_ILC          # one
#   python step2_mallet_lda_all.py --run myeloid Heps
#   python step2_mallet_lda_all.py --force              # refit even if models exist
#
# NOTE on duplicate model files
#   run_cgs_models_mallet(save_path=...) ALREADY writes each model into
#   save_path. Re-pickling the returned models under a second filename in the
#   same directory produces two files per topic count, which then distorts the
#   Step 3 evaluation curves (the same model counted twice) and is why a
#   dedup-by-topic-count step was needed downstream. This script writes the
#   models EXACTLY ONCE -- see RESAVE_MODELS below.
# ==============================================================================

import argparse
import os
import pickle
import time
import traceback


# ==============================================================================
# SHARED CONFIGURATION
# ==============================================================================

# Path to the MALLET binary
MALLET_BIN = "/data/sarkern2/scenicplus1/Mallet-202108/bin/mallet"

# Sampling parameters (shared)
N_CPU        = 20
N_ITER       = 500
RANDOM_STATE = 555
ALPHA        = 50
ETA          = 0.1

# JVM / MALLET memory. Override per compartment for the big atlas run.
JVM_MEMORY = "500g"

# Default topic grid. Compartments (10k-36k cells) rarely need >100 topics;
# the grid is kept wide so the Step 3 metric curves show a clear crossover.
DEFAULT_N_TOPICS = [5, 10, 15, 20, 25, 30, 40, 50, 75, 100]

# run_cgs_models_mallet already saves into save_path. Leave False so models are
# written once. Set True only if your pycisTopic version does NOT auto-save.
RESAVE_MODELS = False


# ==============================================================================
# PER-COMPARTMENT REGISTRY  (out_dir / suffix must match Step 1)
# ==============================================================================

COMPARTMENTS = {

    # full atlas: many more cells -> wider topic grid and more memory
    "all_celltypes": dict(
        out_dir="outs_all_celltypes",
        suffix="all_celltypes",
        n_topics=[25, 50, 75, 100, 150, 200, 250, 300, 350, 400],
        jvm_memory="500g",
    ),

    "Hepatocyte": dict(
        out_dir="outs_Heps",
        suffix="Heps",
        n_topics=DEFAULT_N_TOPICS,
    ),

    "endothelial_Kupffer02": dict(
        out_dir="outs_endothelial_Kupffer02",
        suffix="endothelial_Kupffer02",
        n_topics=DEFAULT_N_TOPICS,
    ),

    "myeloid": dict(
        out_dir="outs_myeloid",
        suffix="myeloid",
        n_topics=DEFAULT_N_TOPICS,
    ),

    "T_ILC": dict(
        out_dir="outs_T_ILC",
        suffix="T_ILC",
        n_topics=DEFAULT_N_TOPICS,
    ),
}


# ==============================================================================
# HELPERS
# ==============================================================================

def set_jvm_memory(mem):
    """Set JVM / MALLET memory. Must happen before the JVM is launched."""
    os.environ["JAVA_TOOL_OPTIONS"] = f"-Xms{mem} -Xmx{mem} -XX:+UseG1GC"
    os.environ["MALLET_MEMORY"] = mem.upper()
    print(f"  JVM memory: {mem}")


def model_topic_count(model):
    """Topic count attribute name differs across pycisTopic versions."""
    for attr in ("n_topic", "topic_no", "n_topics"):
        if hasattr(model, attr):
            return getattr(model, attr)
    return None


def existing_model_topics(result_path):
    """Topic counts already fitted in result_path (for resume)."""
    found = set()
    if not os.path.isdir(result_path):
        return found
    for fn in os.listdir(result_path):
        if not fn.endswith(".pkl"):
            continue
        try:
            with open(os.path.join(result_path, fn), "rb") as fh:
                n = model_topic_count(pickle.load(fh))
            if n is not None:
                found.add(int(n))
        except Exception:
            print(f"  [WARN] could not read existing model {fn}; ignoring")
    return found


# ==============================================================================
# PER-COMPARTMENT RUNNER
# ==============================================================================

def run_compartment(name, spec, force=False):
    from pycisTopic.lda_models import run_cgs_models_mallet

    out_dir = spec["out_dir"]
    suffix = spec["suffix"]
    n_topics = list(spec.get("n_topics", DEFAULT_N_TOPICS))
    tmp_path = os.path.join(out_dir, "mal_temp")
    result_path = os.path.join(out_dir, "mal_result")

    print("\n" + "#" * 70)
    print(f"#  COMPARTMENT: {name}  ->  {result_path}")
    print("#" * 70)

    os.makedirs(tmp_path, exist_ok=True)
    os.makedirs(result_path, exist_ok=True)
    set_jvm_memory(spec.get("jvm_memory", JVM_MEMORY))

    # --- input ---
    cistopic_pkl = os.path.join(out_dir,
                                f"cisTopicObject_filtered_annotated_{suffix}.pkl")
    if not os.path.exists(cistopic_pkl):
        raise FileNotFoundError(
            f"[ERROR] Step 1 output not found: {cistopic_pkl}\n"
            f"        Run the preprocessing pipeline for '{name}' first.")

    print(f"Loading cisTopic object: {cistopic_pkl}")
    with open(cistopic_pkl, "rb") as infile:
        cistopic_obj = pickle.load(infile)
    n_cells = len(cistopic_obj.cell_names)
    n_regions = len(cistopic_obj.region_names)
    print(f"  [OK] {n_cells} cells x {n_regions} regions")

    # --- resume: skip topic counts already fitted ---
    if not force:
        done = existing_model_topics(result_path)
        if done:
            skip = sorted(set(n_topics) & done)
            n_topics = [n for n in n_topics if n not in done]
            if skip:
                print(f"  [INFO] already fitted, skipping: {skip}")
    if not n_topics:
        print("  [OK] nothing to fit (all topic counts present). "
              "Use --force to refit.")
        return {"name": name, "fitted": 0, "cells": n_cells}

    # --- fit ---
    print(f"\nFitting MALLET LDA models for n_topics = {n_topics}")
    print(f"  CPUs:  {N_CPU}")
    print(f"  Iters: {N_ITER}")
    print(f"  Alpha: {ALPHA} (by_topic=True)")
    print(f"  Eta:   {ETA} (by_topic=False)")
    print(f"  Seed:  {RANDOM_STATE}")

    t0 = time.time()
    models = run_cgs_models_mallet(
        mallet_path=MALLET_BIN,
        cistopic_obj=cistopic_obj,
        n_topics=n_topics,
        n_cpu=N_CPU,
        n_iter=N_ITER,
        random_state=RANDOM_STATE,
        alpha=ALPHA,
        alpha_by_topic=True,
        eta=ETA,
        eta_by_topic=False,
        tmp_path=tmp_path,
        save_path=result_path,      # <-- models are written here by the call itself
    )
    elapsed = (time.time() - t0) / 60
    print(f"  [OK] fitting finished in {elapsed:.1f} min")

    # --- optional re-save (OFF by default; see header note on duplicates) ---
    if RESAVE_MODELS and models is not None:
        for model in models:
            n = model_topic_count(model)
            path = os.path.join(result_path, f"mallet_model_{n}.pkl")
            with open(path, "wb") as f:
                pickle.dump(model, f)
            print(f"  [OK] {path}")

    # --- report what is on disk, deduplicated by topic count ---
    on_disk = sorted(existing_model_topics(result_path))
    n_files = len([f for f in os.listdir(result_path) if f.endswith(".pkl")])
    print(f"\n  models on disk: {len(on_disk)} unique topic counts {on_disk} "
          f"({n_files} .pkl files)")
    if n_files > len(on_disk):
        print("  [WARN] more .pkl files than unique topic counts -> duplicate "
              "model files present. Step 3 must dedup by topic count, or clear "
              "mal_result and refit with RESAVE_MODELS = False.")

    return {"name": name, "fitted": len(n_topics), "cells": n_cells}


# ==============================================================================
# MAIN
# ==============================================================================

def main():
    p = argparse.ArgumentParser(
        description="MALLET LDA model fitting for all compartments (SCENIC+ Step 2)")
    p.add_argument("--run", nargs="+", default=None, choices=sorted(COMPARTMENTS),
                   help="Which compartments to fit (default: all)")
    p.add_argument("--force", action="store_true",
                   help="Refit topic counts even if models already exist")
    p.add_argument("--stop-on-error", action="store_true",
                   help="Abort on first failure (default: continue)")
    args = p.parse_args()

    to_run = args.run or list(COMPARTMENTS)
    print("=" * 70)
    print("SCENIC+ Step 2: MALLET LDA fitting")
    print(f"Compartments: {to_run}")
    print(f"MALLET binary: {MALLET_BIN}")
    print("=" * 70)

    if not os.path.exists(MALLET_BIN):
        raise FileNotFoundError(f"[ERROR] MALLET binary not found: {MALLET_BIN}")

    results = []
    for name in to_run:
        try:
            results.append(run_compartment(name, COMPARTMENTS[name], force=args.force))
        except Exception as e:
            print(f"\n[ERROR] {name} failed: {e}")
            traceback.print_exc()
            if args.stop_on_error:
                raise
            print("[INFO] Continuing to next compartment...\n")
            results.append({"name": name, "fitted": None, "cells": None})

    print("\n" + "=" * 70)
    print("STEP 2 FINISHED")
    for r in results:
        if r["fitted"] is None:
            print(f"   {r['name']}: FAILED")
        else:
            print(f"   {r['name']}: {r['fitted']} models fitted "
                  f"({r['cells']} cells)")
    print("=" * 70)
    print("\nNext: Step 3 -- inspect each compartment's evaluation plot and record "
          "the chosen topic count before running downstream analysis.")


if __name__ == "__main__":
    main()
