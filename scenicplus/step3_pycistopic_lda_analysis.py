#!/usr/bin/env python3
"""
SCENIC+ Step 3: LDA model selection, topic analysis, DARs, and SCENIC+ input prep
=================================================================================
One script for all compartments. Consumes Step 1 (annotated cisTopic objects) and
Step 2 (MALLET models), and produces everything SCENIC+ needs.

TWO PHASES -- topic count cannot be automated
---------------------------------------------
The number of topics has to be chosen by eye from the model-quality curves
(Arun / Cao-Juan crossover, Mimno coherence, log-likelihood), and it differs per
compartment. So:

  PHASE 1   python step3_lda_all.py --evaluate
            -> writes <out_dir>/plots/model_evaluation.pdf for each compartment
               and prints the available topic counts. Nothing else runs.

  (you look at each plot, then set n_topics in the COMPARTMENTS registry below)

  PHASE 2   python step3_lda_all.py --run
            -> select model, cluster/UMAP/t-SNE, binarize topics, topic QC +
               annotation, imputation, DARs (age / celltype / sex), region-set
               BEDs, and the final cisTopic object.

Outputs per compartment
-----------------------
  <out_dir>/plots/...                        evaluation, UMAPs, topic QC
  <out_dir>/topic_qc_metrics.csv
  <out_dir>/topic_annotation_by_<var>.csv
  <out_dir>/region_sets/Topics_otsu/*.bed
  <out_dir>/region_sets/Topics_top_3k/*.bed
  <out_dir>/region_sets/DARs_{age,celltype,sex}/*.bed
  <out_dir>/cisTopicObject_lda_complete.pkl

NEXT STEP
---------
RNA-ATAC integration and the SCENIC+ input bundle are NOT done here -- they are
Step 4 (step4_rna_atac_prep_all.py), which consumes
<out_dir>/cisTopicObject_lda_complete.pkl and <out_dir>/region_sets/. Keeping it
separate means the integration can be redone (different RNA object, different
filter, barcode fixes) without repeating LDA selection and DAR calling.

Usage
-----
  python step3_lda_all.py --evaluate
  python step3_lda_all.py --evaluate --run T_ILC
  python step3_lda_all.py --run                      # all compartments
  python step3_lda_all.py --run myeloid --topics 50  # one-off topic override
  python step3_lda_all.py --run T_ILC --skip 7       # skip a step

Author: Nishat Sarker
License: MIT
"""

import argparse
import os
import pickle
import traceback
import warnings

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")


# =============================================================================
# SHARED PARAMETERS
# =============================================================================

N_TOP_REGIONS = 3000                  # binarization: top-N regions per topic
CLUSTER_K = 10                        # k for the topic-space kNN graph
CLUSTER_RESOLUTIONS = [0.6, 1.2, 3]

# find_diff_features parallelism. n_cpu=1 is deliberate: higher values leak Ray
# workers and hang on large imputed matrices. ray.shutdown() is called before
# every invocation for the same reason.
DIFF_N_CPU = 1

MIN_DISP = 0.05
MIN_MEAN = 0.0125
MAX_MEAN = 3
ADJPVAL_THR = 0.05
LOG2FC_THR = np.log2(1.5)

# Variables to compute DARs / topic annotations for, if present in cell_data
DIFF_VARIABLES = ["age", "celltype", "sex"]


# =============================================================================
# PER-COMPARTMENT REGISTRY
#   out_dir / suffix : must match Step 1 and Step 2
#   celltype_col     : annotation column name in the object ('cell_type' for the
#                      full-atlas run); normalized to 'celltype' internally
#   n_topics         : CHOSEN topic count. None until you have looked at the
#                      evaluation plot from PHASE 1.
# (rna_h5ad lives in Step 4's registry, not here)
# =============================================================================

COMPARTMENTS = {

    "all_celltypes": dict(
        out_dir="outs_all_celltypes",
        suffix="all_celltypes",
        celltype_col="cell_type",
        n_topics=None,                       # <-- set after PHASE 1
    ),

    "Hepatocyte": dict(
        out_dir="outs_Heps",
        suffix="Heps",
        celltype_col="celltype",
        n_topics=None,
    ),

    "endothelial_Kupffer02": dict(
        out_dir="outs_endothelial_Kupffer02",
        suffix="endothelial_Kupffer02",
        celltype_col="celltype",
        n_topics=None,
    ),

    "myeloid": dict(
        out_dir="outs_myeloid",
        suffix="myeloid",
        celltype_col="celltype",
        n_topics=None,
    ),

    "T_ILC": dict(
        out_dir="outs_T_ILC",
        suffix="T_ILC",
        celltype_col="celltype",
        n_topics=None,                       # notebook landed on 50 for T cells
    ),
}


# =============================================================================
# CONFIG OBJECT
# =============================================================================

class Cfg:
    def __init__(self, name, spec):
        self.name = name
        self.outDir = spec["out_dir"]
        self.suffix = spec["suffix"]
        self.celltype_col = spec.get("celltype_col", "celltype")
        self.n_topics = spec.get("n_topics")

        self.models_dir = os.path.join(self.outDir, "mal_result")
        self.plots_dir = os.path.join(self.outDir, "plots")
        self.region_sets_dir = os.path.join(self.outDir, "region_sets")

        self.annotated_pkl = os.path.join(
            self.outDir, f"cisTopicObject_filtered_annotated_{self.suffix}.pkl")
        self.lda_added_pkl = os.path.join(self.outDir, "cisTopicObject_lda_added.pkl")
        self.complete_pkl = os.path.join(self.outDir, "cisTopicObject_lda_complete.pkl")

    def makedirs(self):
        os.makedirs(self.plots_dir, exist_ok=True)
        for sub in ["Topics_otsu", "Topics_top_3k"] + \
                   [f"DARs_{v}" for v in DIFF_VARIABLES]:
            os.makedirs(os.path.join(self.region_sets_dir, sub), exist_ok=True)


# =============================================================================
# HELPERS
# =============================================================================

def load_pickle(path, what="object"):
    print(f" Loading {what} from {path}...")
    with open(path, "rb") as f:
        obj = pickle.load(f)
    return obj


def save_pickle(obj, path, what="object"):
    print(f"[SAVE] Saving {what} to {path}...")
    with open(path, "wb") as f:
        pickle.dump(obj, f)
    print("[OK] Saved.")


def model_topic_count(model):
    """Topic-count attribute name differs across pycisTopic versions."""
    for attr in ("n_topic", "topic_no", "n_topics"):
        if hasattr(model, attr):
            return int(getattr(model, attr))
    raise AttributeError("cannot determine topic count from model")


def load_lda_models(models_dir):
    """Load models, DEDUPLICATED BY TOPIC COUNT.

    Step 2's MALLET call saves each model itself; older runs additionally
    re-pickled them under a second name, so a directory can hold two files per
    topic count. Passing both to evaluate_models double-counts the model and
    distorts the metric curves, so collapse by topic count here.
    """
    if not os.path.isdir(models_dir):
        raise FileNotFoundError(f"Models directory not found: {models_dir}")

    files = [f for f in os.listdir(models_dir) if f.endswith(".pkl")]
    if not files:
        raise FileNotFoundError(f"No .pkl model files in {models_dir}")

    uniq = {}
    for fn in sorted(files):
        with open(os.path.join(models_dir, fn), "rb") as fh:
            m = pickle.load(fh)
        uniq[model_topic_count(m)] = m

    models = [uniq[k] for k in sorted(uniq)]
    print(f"[STATS] {len(files)} model files -> {len(models)} unique topic counts: "
          f"{[model_topic_count(m) for m in models]}")
    if len(files) > len(models):
        print("  [INFO] duplicate model files collapsed by topic count")
    return models


def normalize_celltype_col(cistopic_obj, celltype_col, context=""):
    """Rename the annotation column to 'celltype' so downstream code is uniform."""
    cd = cistopic_obj.cell_data
    if celltype_col in cd.columns:
        if celltype_col != "celltype":
            cistopic_obj.cell_data = cd.rename(columns={celltype_col: "celltype"})
            print(f"[OK] {context}renamed '{celltype_col}' -> 'celltype'")
        return cistopic_obj
    for alias in ("celltype", "cell_type", "celltype2"):
        if alias in cd.columns:
            if alias != "celltype":
                cistopic_obj.cell_data = cd.rename(columns={alias: "celltype"})
            print(f"[WARN] {context}celltype_col='{celltype_col}' not found; "
                  f"using '{alias}'")
            return cistopic_obj
    print(f"[WARN] {context}no celltype column found. "
          f"Available: {cd.columns.tolist()}")
    return cistopic_obj


def transform_barcode_index(index):
    """RNA obs_name -> cisTopic cell name.

    Handles both 'ACGT-geriatric_3' and 'ACGT-1-geriatric_3' (the naive
    rsplit('-', 1) used previously turned the latter into 'ACGT-1-1-...' which
    matched nothing). Applies SAMPLE_RENAME and zero-pads the replicate number:

        ACGT-geriatric_3      -> ACGT-1-geriatric_03___geriatric_03
        ACGT-1-pre_geriatric_3-> ACGT-1-pre_ger_03___pre_ger_03
    """
    s = str(index).split("___")[0]
    parts = s.split("-")
    dna = parts[0]
    rest = [p for p in parts[1:] if p != "1"]
    if not rest:
        raise ValueError(f"cannot parse sample from index: {index}")
    sample = "-".join(rest)
    for old, new in SAMPLE_RENAME.items():
        sample = re.sub(rf"^{old}", new, sample)
    sample = re.sub(r"_(\d+)$", lambda m: f"_{int(m.group(1)):02d}", sample)
    return f"{dna}-1-{sample}___{sample}"


def dump_region_sets(region_dict, out_folder):
    """Write each region set to a BED file."""
    from pycisTopic.utils import region_names_to_coordinates
    os.makedirs(out_folder, exist_ok=True)
    n = 0
    for name, regions in region_dict.items():
        idx = regions.index if hasattr(regions, "index") else regions
        if len(idx) == 0:
            continue
        region_names_to_coordinates(idx).sort_values(
            ["Chromosome", "Start", "End"]).to_csv(
            os.path.join(out_folder, f"{name}.bed"),
            sep="\t", header=False, index=False)
        n += 1
    print(f"   [OK] wrote {n} BED files to {out_folder}")


def ray_reset():
    """Clear any leaked Ray instance before find_diff_features."""
    try:
        import ray
        ray.shutdown()
    except Exception:
        pass


# =============================================================================
# PHASE 1: EVALUATE MODELS
# =============================================================================

def phase_evaluate(cfg):
    """Draw the 4-metric plot so the topic count can be chosen by eye."""
    from pycisTopic.lda_models import evaluate_models

    print("\n" + "#" * 70)
    print(f"#  EVALUATE: {cfg.name}")
    print("#" * 70)
    cfg.makedirs()

    if not os.path.exists(cfg.annotated_pkl):
        raise FileNotFoundError(
            f"[ERROR] Step 1 output missing: {cfg.annotated_pkl}")

    cistopic_obj = load_pickle(cfg.annotated_pkl, "cisTopic object")
    print(f"  cells: {len(cistopic_obj.cell_names)}  "
          f"regions: {len(cistopic_obj.region_names)}")

    models = load_lda_models(cfg.models_dir)

    save_path = os.path.join(cfg.plots_dir, "model_evaluation.pdf")
    evaluate_models(models, select_model=None, return_model=False, save=save_path)
    plt.close("all")
    print(f"[OK] evaluation plot -> {save_path}")
    print(f"[ACTION] inspect the plot, then set n_topics for '{cfg.name}' "
          f"in COMPARTMENTS (available: "
          f"{[model_topic_count(m) for m in models]})")


# =============================================================================
# PHASE 2 STEPS
# =============================================================================

def step1_load(cfg):
    """Load the compartment object. No fallback -- the merged pre-Step-13 object
    contains every QC-passing barcode, not just this compartment."""
    print("\n" + "=" * 60)
    print("STEP 1: Loading cisTopic object")
    print("=" * 60)
    if not os.path.exists(cfg.annotated_pkl):
        raise FileNotFoundError(
            f"[ERROR] Step 1 output missing: {cfg.annotated_pkl}\n"
            f"        (Do NOT substitute cisTopicObject_merged_dbl_filtered.pkl: "
            f"it still holds all cell types.)")
    obj = load_pickle(cfg.annotated_pkl, "cisTopic object")
    obj = normalize_celltype_col(obj, cfg.celltype_col)
    print(f"   Cells: {len(obj.cell_names)}")
    print(f"   Regions: {len(obj.region_names)}")
    print(f"   Metadata: {obj.cell_data.columns.tolist()}")
    return obj


def step2_select_model(cfg, cistopic_obj, n_topics):
    print("\n" + "=" * 60)
    print(f"STEP 2: Selecting LDA model ({n_topics} topics)")
    print("=" * 60)
    from pycisTopic.lda_models import evaluate_models

    models = load_lda_models(cfg.models_dir)
    available = [model_topic_count(m) for m in models]
    if n_topics not in available:
        raise ValueError(f"[ERROR] n_topics={n_topics} not fitted for {cfg.name}. "
                         f"Available: {available}")

    selected = evaluate_models(models, select_model=n_topics, return_model=True)
    if selected is None:                       # some versions return None
        selected = models[available.index(n_topics)]
    plt.close("all")

    cistopic_obj.add_LDA_model(selected)
    print(f"[OK] added model with {model_topic_count(selected)} topics")
    save_pickle(cistopic_obj, cfg.lda_added_pkl, "object with LDA")
    return cistopic_obj


def step3_cluster(cfg, cistopic_obj):
    print("\n" + "=" * 60)
    print("STEP 3: Clustering and dimensionality reduction")
    print("=" * 60)
    from pycisTopic.clust_vis import (find_clusters, run_umap, run_tsne,
                                      plot_metadata, plot_topic, cell_topic_heatmap)

    find_clusters(cistopic_obj, target="cell", k=CLUSTER_K, res=CLUSTER_RESOLUTIONS,
                  prefix="pycisTopic_", scale=True, split_pattern="-")
    run_umap(cistopic_obj, target="cell", scale=True)
    run_tsne(cistopic_obj, target="cell", scale=True)
    print("[OK] clustering + UMAP + t-SNE complete")

    meta_vars = [v for v in ["age", "sex", "sample_id", "celltype"]
                 if v in cistopic_obj.cell_data.columns]
    cluster_cols = [c for c in cistopic_obj.cell_data.columns
                    if c.startswith("pycisTopic_")]

    for var in meta_vars + cluster_cols:
        try:
            plot_metadata(cistopic_obj, reduction_name="UMAP", variables=[var],
                          target="cell", num_columns=1, text_size=10, dot_size=5)
            plt.savefig(os.path.join(cfg.plots_dir, f"umap_{var}.png"),
                        dpi=150, bbox_inches="tight")
            plt.savefig(os.path.join(cfg.plots_dir, f"umap_{var}.pdf"),
                        bbox_inches="tight")
            plt.close("all")
            print(f"   [OK] umap_{var}")
        except Exception as e:
            print(f"   [WARN] could not plot {var}: {e}")

    # topic contributions
    try:
        n_topic = model_topic_count(cistopic_obj.selected_model)
        plot_topic(cistopic_obj, reduction_name="UMAP", target="cell",
                   num_columns=5, topics=list(range(1, min(20, n_topic) + 1)))
        plt.savefig(os.path.join(cfg.plots_dir, "umap_topics.png"),
                    dpi=150, bbox_inches="tight")
        plt.close("all")
        print("   [OK] umap_topics")
    except Exception as e:
        print(f"   [WARN] could not plot topics: {e}")

    # cell-topic heatmap
    try:
        hm_vars = [v for v in ["celltype", "age"]
                   if v in cistopic_obj.cell_data.columns]
        cell_topic_heatmap(cistopic_obj, variables=hm_vars or None, scale=True,
                           legend_loc_x=1.05, legend_loc_y=-0.5, legend_dist_y=-1,
                           figsize=(20, 10))
        plt.savefig(os.path.join(cfg.plots_dir, "cell_topic_heatmap.png"),
                    dpi=150, bbox_inches="tight")
        plt.savefig(os.path.join(cfg.plots_dir, "cell_topic_heatmap.pdf"),
                    bbox_inches="tight")
        plt.close("all")
        print("   [OK] cell_topic_heatmap")
    except Exception as e:
        print(f"   [WARN] could not generate heatmap: {e}")

    return cistopic_obj


def step4_binarize(cfg, cistopic_obj):
    print("\n" + "=" * 60)
    print("STEP 4: Topic binarization")
    print("=" * 60)
    from pycisTopic.topic_binarization import binarize_topics

    top3k = binarize_topics(cistopic_obj, method="ntop", ntop=N_TOP_REGIONS,
                            plot=True, num_columns=5)
    plt.savefig(os.path.join(cfg.plots_dir, "topic_binarization_top3k.png"),
                dpi=150, bbox_inches="tight")
    plt.close("all")
    print(f"   [OK] {len(top3k)} topics (top {N_TOP_REGIONS})")

    otsu = binarize_topics(cistopic_obj, method="otsu", plot=True, num_columns=5)
    plt.savefig(os.path.join(cfg.plots_dir, "topic_binarization_otsu.png"),
                dpi=150, bbox_inches="tight")
    plt.close("all")
    print(f"   [OK] {len(otsu)} topics (Otsu)")

    cell_topic = binarize_topics(cistopic_obj, target="cell", method="li",
                                 plot=True, num_columns=5, nbins=100)
    plt.savefig(os.path.join(cfg.plots_dir, "cell_topic_binarization_li.png"),
                dpi=150, bbox_inches="tight")
    plt.close("all")
    print("   [OK] cell-topic matrix (Li)")

    dump_region_sets(top3k, os.path.join(cfg.region_sets_dir, "Topics_top_3k"))
    dump_region_sets(otsu, os.path.join(cfg.region_sets_dir, "Topics_otsu"))
    return top3k, otsu, cell_topic


def step5_topic_qc(cfg, cistopic_obj, binarized_cell_topic):
    print("\n" + "=" * 60)
    print("STEP 5: Topic QC and annotation")
    print("=" * 60)
    from pycisTopic.topic_qc import (compute_topic_metrics, plot_topic_qc,
                                     topic_annotation)

    metrics = compute_topic_metrics(cistopic_obj)
    metrics.to_csv(os.path.join(cfg.outDir, "topic_qc_metrics.csv"))
    print("   [OK] topic_qc_metrics.csv")

    try:
        plot_topic_qc(metrics, num_columns=4)
        plt.savefig(os.path.join(cfg.plots_dir, "topic_qc.png"),
                    dpi=150, bbox_inches="tight")
        plt.close("all")
        print("   [OK] topic_qc.png")
    except Exception as e:
        print(f"   [WARN] could not plot topic QC: {e}")

    for var in DIFF_VARIABLES:
        if var not in cistopic_obj.cell_data.columns:
            print(f"   [WARN] '{var}' not in cell_data - skipping annotation")
            continue
        try:
            topic_annotation(cistopic_obj, annot_var=var,
                             binarized_cell_topic=binarized_cell_topic,
                             general_topic_thr=0.2).to_csv(
                os.path.join(cfg.outDir, f"topic_annotation_by_{var}.csv"))
            print(f"   [OK] topic_annotation_by_{var}.csv")
        except Exception as e:
            print(f"   [WARN] could not annotate by {var}: {e}")

    return metrics


def step6_dars(cfg, cistopic_obj):
    print("\n" + "=" * 60)
    print("STEP 6: Imputation and differential accessibility")
    print("=" * 60)
    from pycisTopic.diff_features import (impute_accessibility, normalize_scores,
                                          find_highly_variable_features,
                                          find_diff_features)
    from pycisTopic.utils import region_names_to_coordinates

    imputed = impute_accessibility(cistopic_obj, selected_cells=None,
                                   selected_regions=None, scale_factor=10 ** 6)
    print(f"   [OK] imputed: {imputed.mtx.shape if hasattr(imputed,'mtx') else ''}")

    normalized = normalize_scores(imputed, scale_factor=10 ** 4)
    print("   [OK] normalized")

    variable_regions = find_highly_variable_features(
        normalized, min_disp=MIN_DISP, min_mean=MIN_MEAN, max_mean=MAX_MEAN,
        max_disp=np.inf, n_bins=20, n_top_features=None, plot=True)
    plt.savefig(os.path.join(cfg.plots_dir, "highly_variable_regions.png"),
                dpi=150, bbox_inches="tight")
    plt.close("all")
    print(f"   [OK] {len(variable_regions)} highly variable regions")

    region_names_to_coordinates(variable_regions).sort_values(
        ["Chromosome", "Start", "End"]).to_csv(
        os.path.join(cfg.region_sets_dir, "highly_variable_regions.bed"),
        sep="\t", header=False, index=False)

    all_markers = {}
    for var in DIFF_VARIABLES:
        if var not in cistopic_obj.cell_data.columns:
            print(f"\n   [WARN] '{var}' not in cell_data - skipping DARs")
            continue
        print(f"\n[STATS] DARs by {var}...")
        try:
            ray_reset()   # leaked Ray instances hang the next call
            markers = find_diff_features(
                cistopic_obj, imputed, variable=var,
                var_features=variable_regions, contrasts=None,
                adjpval_thr=ADJPVAL_THR, log2fc_thr=LOG2FC_THR,
                n_cpu=DIFF_N_CPU, _temp_dir=None, split_pattern="-")
            all_markers[var] = markers

            for group, df in markers.items():
                df.to_csv(os.path.join(cfg.outDir, f"DARs_{var}_{group}.csv"))
            dump_region_sets(markers, os.path.join(cfg.region_sets_dir, f"DARs_{var}"))

            print(f"   DAR summary ({var}):")
            for group, df in markers.items():
                if "Log2FC" in df.columns:
                    up, dn = int((df["Log2FC"] > 0).sum()), int((df["Log2FC"] < 0).sum())
                    print(f"     {group}: {len(df)} DARs (up {up}, down {dn})")
                else:
                    print(f"     {group}: {len(df)} DARs")
        except Exception as e:
            print(f"   [WARN] error finding {var} DARs: {e}")

    return imputed, variable_regions, all_markers


def step7_plot_features(cfg, cistopic_obj, imputed, all_markers):
    print("\n" + "=" * 60)
    print("STEP 7: Plotting top imputed features")
    print("=" * 60)
    from pycisTopic.clust_vis import plot_imputed_features

    if not all_markers:
        print("[WARN] no markers - skipping")
        return

    for var, markers in all_markers.items():
        for group, df in markers.items():
            if len(df) == 0:
                continue
            if "Log2FC" in df.columns:
                top = df.reindex(df["Log2FC"].abs().sort_values(
                    ascending=False).index).head(6).index.tolist()
            else:
                top = df.head(6).index.tolist()
            try:
                plot_imputed_features(cistopic_obj, reduction_name="UMAP",
                                      imputed_acc_obj=imputed, features=top,
                                      scale=True, num_columns=3)
                plt.savefig(os.path.join(
                    cfg.plots_dir, f"imputed_features_{var}_{group}.png"),
                    dpi=150, bbox_inches="tight")
                plt.close("all")
                print(f"   [OK] imputed_features_{var}_{group}")
            except Exception as e:
                print(f"   [WARN] {var}/{group}: {e}")


def step8_save(cfg, cistopic_obj):
    print("\n" + "=" * 60)
    print("STEP 8: Saving final cisTopic object")
    print("=" * 60)
    save_pickle(cistopic_obj, cfg.complete_pkl, "final object")
    print(f"Cells: {len(cistopic_obj.cell_names)}")
    print(f"Regions: {len(cistopic_obj.region_names)}")
    print(f"Topics: {model_topic_count(cistopic_obj.selected_model)}")
    return cfg.complete_pkl


# =============================================================================
# PHASE 2 DRIVER
# =============================================================================

def phase_run(cfg, n_topics_override=None, skip_steps=()):
    print("\n" + "#" * 70)
    print(f"#  RUN: {cfg.name}  ->  {cfg.outDir}")
    print("#" * 70)
    cfg.makedirs()

    n_topics = n_topics_override or cfg.n_topics
    if n_topics is None:
        raise ValueError(
            f"[ERROR] no topic count set for '{cfg.name}'.\n"
            f"        Run --evaluate, inspect {cfg.plots_dir}/model_evaluation.pdf, "
            f"then set n_topics in COMPARTMENTS (or pass --topics).")

    cistopic_obj = step1_load(cfg)
    cistopic_obj = step2_select_model(cfg, cistopic_obj, n_topics)

    if 3 not in skip_steps:
        cistopic_obj = step3_cluster(cfg, cistopic_obj)

    binarized_cell_topic = None
    if 4 not in skip_steps:
        _, _, binarized_cell_topic = step4_binarize(cfg, cistopic_obj)

    if 5 not in skip_steps and binarized_cell_topic is not None:
        step5_topic_qc(cfg, cistopic_obj, binarized_cell_topic)

    imputed, all_markers = None, {}
    if 6 not in skip_steps:
        imputed, _, all_markers = step6_dars(cfg, cistopic_obj)

    if 7 not in skip_steps and imputed is not None:
        step7_plot_features(cfg, cistopic_obj, imputed, all_markers)

    step8_save(cfg, cistopic_obj)

    print("\n" + "=" * 60)
    print(f"[OK] {cfg.name} COMPLETE ({n_topics} topics)")
    print(f"     next: step4_rna_atac_prep_all.py --run {cfg.name}")
    print("=" * 60)
    return cistopic_obj


# =============================================================================
# MAIN
# =============================================================================

def main():
    p = argparse.ArgumentParser(
        description="SCENIC+ Step 3: LDA model selection, topic analysis, and DARs")
    p.add_argument("--evaluate", action="store_true",
                   help="PHASE 1: draw model-evaluation plots and exit")
    p.add_argument("--run", nargs="*", default=None, metavar="COMPARTMENT",
                   help="PHASE 2: run downstream analysis "
                        "(no names = all compartments)")
    p.add_argument("--topics", type=int, default=None,
                   help="Override n_topics (only sensible with a single --run)")
    p.add_argument("--skip", type=int, nargs="+", default=[],
                   help="Step numbers to skip (3-7), e.g. --skip 7")
    p.add_argument("--stop-on-error", action="store_true",
                   help="Abort on first failure (default: continue)")
    args = p.parse_args()

    if not args.evaluate and args.run is None:
        p.error("nothing to do: pass --evaluate and/or --run")

    names = args.run if args.run else list(COMPARTMENTS)
    unknown = [n for n in names if n not in COMPARTMENTS]
    if unknown:
        p.error(f"unknown compartment(s): {unknown}. "
                f"Choose from {sorted(COMPARTMENTS)}")

    if args.topics is not None and len(names) > 1:
        p.error("--topics applies to a single compartment; "
                "set n_topics in COMPARTMENTS for multi-compartment runs")

    if args.evaluate:
        for name in names:
            try:
                phase_evaluate(Cfg(name, COMPARTMENTS[name]))
            except Exception as e:
                print(f"\n[ERROR] evaluate {name} failed: {e}")
                traceback.print_exc()
                if args.stop_on_error:
                    raise
        if args.run is None:
            print("\n[NEXT] set n_topics per compartment, then re-run with --run")
            return

    if args.run is not None:
        results = {}
        for name in names:
            try:
                phase_run(Cfg(name, COMPARTMENTS[name]),
                          n_topics_override=args.topics,
                          skip_steps=set(args.skip))
                results[name] = "OK"
            except Exception as e:
                print(f"\n[ERROR] {name} failed: {e}")
                traceback.print_exc()
                if args.stop_on_error:
                    raise
                results[name] = "FAILED"
        print("\n" + "=" * 60)
        print("STEP 3 FINISHED")
        for k, v in results.items():
            print(f"   {k}: {v}")
        print("=" * 60)


if __name__ == "__main__":
    main()
