#!/usr/bin/env python3
"""
pycisTopic LDA Topic Modeling and Downstream Analysis Pipeline
This script performs topic modeling, clustering, binarization, and differential
accessibility analysis on a preprocessed cisTopic object.

Author: Nishat Sarker
Date: 2025

NOTE: Do not run in exclusive mode - this script benefits from parallel processing
"""

# =============================================================================
# IMPORTS
# =============================================================================

# === Standard Library ===
import os
import pickle
import warnings

# === Scientific Computing ===
import numpy as np
import pandas as pd
import polars as pl

# === Plotting ===
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

# === Scanpy / AnnData ===
import scanpy as sc
import anndata as ad

# === pycisTopic Modules ===
from pycisTopic.lda_models import evaluate_models, run_cgs_models
from pycisTopic.clust_vis import (
    find_clusters,
    run_umap,
    run_tsne,
    plot_metadata,
    plot_topic,
    cell_topic_heatmap
)
from pycisTopic.topic_binarization import binarize_topics
from pycisTopic.topic_qc import compute_topic_metrics, plot_topic_qc, topic_annotation
from pycisTopic.diff_features import (
    impute_accessibility,
    normalize_scores,
    find_highly_variable_features,
    find_diff_features
)
from pycisTopic.clust_vis import plot_imputed_features
from pycisTopic.utils import region_names_to_coordinates

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')


# =============================================================================
# CONFIGURATION
# =============================================================================

# Directory settings
out_dir = "outs_trial"
models_dir = os.path.join(out_dir, "mal_result")
plots_dir = os.path.join(out_dir, "plots")
region_sets_dir = os.path.join(out_dir, "region_sets")

# Create output directories
os.makedirs(plots_dir, exist_ok=True)
os.makedirs(os.path.join(region_sets_dir, "Topics_otsu"), exist_ok=True)
os.makedirs(os.path.join(region_sets_dir, "Topics_top_3k"), exist_ok=True)
os.makedirs(os.path.join(region_sets_dir, "DARs_age"), exist_ok=True)

# Analysis parameters
SELECTED_TOPICS = 350  # Number of topics to select
N_TOP_REGIONS = 3000   # Number of top regions for binarization
CLUSTER_K = 10         # k for clustering
CLUSTER_RESOLUTIONS = [0.6, 1.2, 3]  # Clustering resolutions
N_CPU = 5              # Number of CPUs for parallel processing

# Differential accessibility parameters
MIN_DISP = 0.05
MIN_MEAN = 0.0125
MAX_MEAN = 3
ADJPVAL_THR = 0.05
LOG2FC_THR = np.log2(1.5)


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def load_cistopic_object(filepath):
    """Load a cisTopic object from pickle file."""
    print(f"📂 Loading cisTopic object from {filepath}...")
    with open(filepath, 'rb') as f:
        cistopic_obj = pickle.load(f)
    print(f"✅ Loaded successfully. Shape: {cistopic_obj.fragment_matrix.shape}")
    return cistopic_obj


def save_cistopic_object(cistopic_obj, filepath):
    """Save a cisTopic object to pickle file."""
    print(f"💾 Saving cisTopic object to {filepath}...")
    with open(filepath, 'wb') as f:
        pickle.dump(cistopic_obj, f)
    print(f"✅ Saved successfully.")


def load_lda_models(models_dir):
    """Load all LDA models from a directory."""
    print(f"\n📂 Loading LDA models from {models_dir}...")
    
    if not os.path.exists(models_dir):
        raise FileNotFoundError(f"Models directory not found: {models_dir}")
    
    model_files = [f for f in os.listdir(models_dir) if f.endswith('.pkl')]
    
    if not model_files:
        raise FileNotFoundError(f"No .pkl model files found in {models_dir}")
    
    models = {}
    for file_name in model_files:
        label = os.path.splitext(file_name)[0]
        model_path = os.path.join(models_dir, file_name)
        with open(model_path, "rb") as f:
            models[label] = pickle.load(f)
        print(f"   ✅ Loaded model: {label}")
    
    print(f"📊 Total models loaded: {len(models)}")
    return models


def save_region_sets_to_bed(region_dict, output_dir, prefix=""):
    """Save region sets to BED files."""
    for name, regions in region_dict.items():
        if hasattr(regions, 'index'):
            region_coords = region_names_to_coordinates(regions.index)
        else:
            region_coords = region_names_to_coordinates(regions)
        
        region_coords.sort_values(
            ["Chromosome", "Start", "End"]
        ).to_csv(
            os.path.join(output_dir, f"{prefix}{name}.bed"),
            sep="\t",
            header=False,
            index=False
        )
    print(f"✅ Saved {len(region_dict)} region sets to {output_dir}")


# =============================================================================
# STEP 1: LOAD CISTOPIC OBJECT
# =============================================================================

def step1_load_data():
    """Load the preprocessed cisTopic object."""
    print("\n" + "="*60)
    print("STEP 1: Loading cisTopic Object")
    print("="*60)
    
    # Try different possible input files
    possible_inputs = [
        os.path.join(out_dir, 'cisTopicObject_filtered_annotated.pkl'),
        os.path.join(out_dir, 'cisTopicObject_filtered_annotated_Heps.pkl'),
        os.path.join(out_dir, 'cisTopicObject_merged_dbl_filtered.pkl')
    ]
    
    cistopic_obj = None
    for filepath in possible_inputs:
        if os.path.exists(filepath):
            cistopic_obj = load_cistopic_object(filepath)
            break
    
    if cistopic_obj is None:
        raise FileNotFoundError(f"No cisTopic object found. Tried: {possible_inputs}")
    
    # Print summary
    print(f"\n📊 cisTopic Object Summary:")
    print(f"   Cells: {len(cistopic_obj.cell_names)}")
    print(f"   Regions: {len(cistopic_obj.region_names)}")
    print(f"   Metadata columns: {cistopic_obj.cell_data.columns.tolist()}")
    
    return cistopic_obj


# =============================================================================
# STEP 2: LOAD AND EVALUATE LDA MODELS
# =============================================================================

def step2_evaluate_models(cistopic_obj):
    """Load pre-trained LDA models and select the best one."""
    print("\n" + "="*60)
    print("STEP 2: Loading and Evaluating LDA Models")
    print("="*60)
    
    # Load models
    models = load_lda_models(models_dir)
    model_list = list(models.values())
    
    # Verify model attributes
    print("\n🔍 Model verification:")
    for label, model in list(models.items())[:3]:  # Show first 3
        print(f"   {label}: {model.n_topic} topics")
    
    # Evaluate models and select best
    print(f"\n📈 Evaluating models (selecting {SELECTED_TOPICS} topics)...")
    
    selected_model = evaluate_models(
        model_list,
        select_model=SELECTED_TOPICS,
        return_model=True,
        save=os.path.join(plots_dir, 'model_evaluation.pdf')
    )
    
    if selected_model is not None:
        print(f"✅ Selected model with {selected_model.n_topic} topics")
    else:
        # If no model returned, find the one with closest topic count
        print(f"⚠️ evaluate_models did not return a model. Selecting manually...")
        topic_counts = {label: model.n_topic for label, model in models.items()}
        closest_label = min(topic_counts, key=lambda x: abs(topic_counts[x] - SELECTED_TOPICS))
        selected_model = models[closest_label]
        print(f"✅ Selected model '{closest_label}' with {selected_model.n_topic} topics")
    
    # Add model to cisTopic object
    print("\n🔗 Adding LDA model to cisTopic object...")
    cistopic_obj.add_LDA_model(selected_model)
    
    # Save intermediate result
    save_cistopic_object(cistopic_obj, os.path.join(out_dir, 'cisTopicObject_lda_added.pkl'))
    
    return cistopic_obj, selected_model


# =============================================================================
# STEP 3: CLUSTERING AND DIMENSIONALITY REDUCTION
# =============================================================================

def step3_clustering_and_visualization(cistopic_obj):
    """Perform clustering and dimensionality reduction."""
    print("\n" + "="*60)
    print("STEP 3: Clustering and Dimensionality Reduction")
    print("="*60)
    
    # Find clusters
    print(f"\n🔍 Finding clusters (k={CLUSTER_K}, resolutions={CLUSTER_RESOLUTIONS})...")
    find_clusters(
        cistopic_obj,
        target='cell',
        k=CLUSTER_K,
        res=CLUSTER_RESOLUTIONS,
        prefix='pycisTopic_',
        scale=True,
        split_pattern='-'
    )
    print("✅ Clustering complete")
    
    # Run UMAP
    print("\n🗺️ Running UMAP...")
    run_umap(
        cistopic_obj,
        target='cell',
        scale=True
    )
    print("✅ UMAP complete")
    
    # Run t-SNE
    print("\n🗺️ Running t-SNE...")
    run_tsne(
        cistopic_obj,
        target='cell',
        scale=True
    )
    print("✅ t-SNE complete")
    
    # Generate visualization plots
    print("\n📊 Generating visualization plots...")
    
    # Plot UMAP colored by metadata
    metadata_vars = ['age', 'sex', 'sample_id', 'celltype']
    available_vars = [v for v in metadata_vars if v in cistopic_obj.cell_data.columns]
    
    for var in available_vars:
        try:
            fig = plot_metadata(
                cistopic_obj,
                reduction_name='UMAP',
                variables=[var],
                target='cell',
                num_columns=1,
                text_size=10,
                dot_size=5
            )
            plt.savefig(os.path.join(plots_dir, f'umap_{var}.png'), dpi=150, bbox_inches='tight')
            plt.savefig(os.path.join(plots_dir, f'umap_{var}.pdf'), bbox_inches='tight')
            plt.close()
            print(f"   ✅ Saved UMAP colored by {var}")
        except Exception as e:
            print(f"   ⚠️ Could not plot {var}: {e}")
    
    # Plot UMAP colored by cluster
    cluster_cols = [c for c in cistopic_obj.cell_data.columns if c.startswith('pycisTopic_')]
    for cluster_col in cluster_cols:
        try:
            fig = plot_metadata(
                cistopic_obj,
                reduction_name='UMAP',
                variables=[cluster_col],
                target='cell',
                num_columns=1,
                text_size=10,
                dot_size=5
            )
            plt.savefig(os.path.join(plots_dir, f'umap_{cluster_col}.png'), dpi=150, bbox_inches='tight')
            plt.close()
            print(f"   ✅ Saved UMAP colored by {cluster_col}")
        except Exception as e:
            print(f"   ⚠️ Could not plot {cluster_col}: {e}")
    
    # Plot topic contributions on UMAP
    print("\n📊 Plotting topic contributions...")
    n_topics = cistopic_obj.selected_model.n_topic
    topics_to_plot = min(20, n_topics)  # Plot first 20 topics
    
    try:
        fig = plot_topic(
            cistopic_obj,
            reduction_name='UMAP',
            target='cell',
            num_columns=5,
            topics=list(range(1, topics_to_plot + 1))
        )
        plt.savefig(os.path.join(plots_dir, 'umap_topics.png'), dpi=150, bbox_inches='tight')
        plt.close()
        print(f"   ✅ Saved topic UMAP plots")
    except Exception as e:
        print(f"   ⚠️ Could not plot topics: {e}")
    
    # Cell-topic heatmap
    print("\n📊 Generating cell-topic heatmap...")
    try:
        fig = cell_topic_heatmap(
            cistopic_obj,
            variables=['age'] if 'age' in cistopic_obj.cell_data.columns else None,
            scale=True,
            legend_loc_x=1.05,
            legend_loc_y=-0.5,
            legend_dist_y=-1,
            figsize=(20, 10)
        )
        plt.savefig(os.path.join(plots_dir, 'cell_topic_heatmap.png'), dpi=150, bbox_inches='tight')
        plt.savefig(os.path.join(plots_dir, 'cell_topic_heatmap.pdf'), bbox_inches='tight')
        plt.close()
        print("   ✅ Saved cell-topic heatmap")
    except Exception as e:
        print(f"   ⚠️ Could not generate heatmap: {e}")
    
    # Save intermediate result
    save_cistopic_object(cistopic_obj, os.path.join(out_dir, 'cisTopicObject_clustered.pkl'))
    
    return cistopic_obj


# =============================================================================
# STEP 4: TOPIC BINARIZATION
# =============================================================================

def step4_topic_binarization(cistopic_obj):
    """Binarize topics using multiple methods."""
    print("\n" + "="*60)
    print("STEP 4: Topic Binarization")
    print("="*60)
    
    # Method 1: Top N regions
    print(f"\n📊 Binarizing topics (top {N_TOP_REGIONS} regions)...")
    region_bin_topics_top_3k = binarize_topics(
        cistopic_obj,
        method='ntop',
        ntop=N_TOP_REGIONS,
        plot=True,
        num_columns=5
    )
    plt.savefig(os.path.join(plots_dir, 'topic_binarization_top3k.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print(f"   ✅ Binarized {len(region_bin_topics_top_3k)} topics (top {N_TOP_REGIONS})")
    
    # Method 2: Otsu thresholding
    print("\n📊 Binarizing topics (Otsu method)...")
    region_bin_topics_otsu = binarize_topics(
        cistopic_obj,
        method='otsu',
        plot=True,
        num_columns=5
    )
    plt.savefig(os.path.join(plots_dir, 'topic_binarization_otsu.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print(f"   ✅ Binarized {len(region_bin_topics_otsu)} topics (Otsu)")
    
    # Method 3: Cell-topic binarization (Li method)
    print("\n📊 Binarizing cell-topic matrix (Li method)...")
    binarized_cell_topic = binarize_topics(
        cistopic_obj,
        target='cell',
        method='li',
        plot=True,
        num_columns=5,
        nbins=100
    )
    plt.savefig(os.path.join(plots_dir, 'cell_topic_binarization_li.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print(f"   ✅ Binarized cell-topic matrix")
    
    # Save region sets to BED files
    print("\n💾 Saving binarized region sets...")
    save_region_sets_to_bed(
        region_bin_topics_top_3k,
        os.path.join(region_sets_dir, "Topics_top_3k")
    )
    save_region_sets_to_bed(
        region_bin_topics_otsu,
        os.path.join(region_sets_dir, "Topics_otsu")
    )
    
    return region_bin_topics_top_3k, region_bin_topics_otsu, binarized_cell_topic


# =============================================================================
# STEP 5: TOPIC QC AND ANNOTATION
# =============================================================================

def step5_topic_qc(cistopic_obj, binarized_cell_topic):
    """Compute topic QC metrics and annotations."""
    print("\n" + "="*60)
    print("STEP 5: Topic QC and Annotation")
    print("="*60)
    
    # Compute topic metrics
    print("\n📊 Computing topic QC metrics...")
    topic_qc_metrics = compute_topic_metrics(cistopic_obj)
    
    # Save metrics
    topic_qc_metrics.to_csv(os.path.join(out_dir, 'topic_qc_metrics.csv'))
    print(f"   ✅ Saved topic QC metrics")
    
    # Plot topic QC
    try:
        fig = plot_topic_qc(topic_qc_metrics, num_columns=4)
        plt.savefig(os.path.join(plots_dir, 'topic_qc.png'), dpi=150, bbox_inches='tight')
        plt.close()
        print(f"   ✅ Saved topic QC plots")
    except Exception as e:
        print(f"   ⚠️ Could not plot topic QC: {e}")
    
    # Topic annotation by age
    print("\n📊 Annotating topics by age...")
    topic_annot = None
    if 'age' in cistopic_obj.cell_data.columns:
        try:
            topic_annot = topic_annotation(
                cistopic_obj,
                annot_var='age',
                binarized_cell_topic=binarized_cell_topic,
                general_topic_thr=0.2
            )
            topic_annot.to_csv(os.path.join(out_dir, 'topic_annotation_by_age.csv'))
            print(f"   ✅ Saved topic annotations")
        except Exception as e:
            print(f"   ⚠️ Could not annotate topics: {e}")
    else:
        print("   ⚠️ 'age' column not found in cell_data")
    
    # Topic annotation by celltype (if available)
    if 'celltype' in cistopic_obj.cell_data.columns:
        print("\n📊 Annotating topics by celltype...")
        try:
            topic_annot_celltype = topic_annotation(
                cistopic_obj,
                annot_var='celltype',
                binarized_cell_topic=binarized_cell_topic,
                general_topic_thr=0.2
            )
            topic_annot_celltype.to_csv(os.path.join(out_dir, 'topic_annotation_by_celltype.csv'))
            print(f"   ✅ Saved celltype topic annotations")
        except Exception as e:
            print(f"   ⚠️ Could not annotate topics by celltype: {e}")
    
    return topic_qc_metrics, topic_annot


# =============================================================================
# STEP 6: IMPUTATION AND DIFFERENTIAL ACCESSIBILITY
# =============================================================================

def step6_differential_accessibility(cistopic_obj):
    """Perform accessibility imputation and differential analysis."""
    print("\n" + "="*60)
    print("STEP 6: Imputation and Differential Accessibility Analysis")
    print("="*60)
    
    # Impute accessibility
    print("\n🔬 Imputing accessibility scores...")
    imputed_acc_obj = impute_accessibility(
        cistopic_obj,
        selected_cells=None,
        selected_regions=None,
        scale_factor=10**6
    )
    print(f"   ✅ Imputation complete. Shape: {imputed_acc_obj.X.shape}")
    
    # Normalize scores
    print("\n📊 Normalizing imputed scores...")
    normalized_imputed_acc_obj = normalize_scores(
        imputed_acc_obj,
        scale_factor=10**4
    )
    print(f"   ✅ Normalization complete")
    
    # Find highly variable features
    print("\n🔍 Finding highly variable regions...")
    variable_regions = find_highly_variable_features(
        normalized_imputed_acc_obj,
        min_disp=MIN_DISP,
        min_mean=MIN_MEAN,
        max_mean=MAX_MEAN,
        max_disp=np.inf,
        n_bins=20,
        n_top_features=None,
        plot=True
    )
    plt.savefig(os.path.join(plots_dir, 'highly_variable_regions.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print(f"   ✅ Found {len(variable_regions)} highly variable regions")
    
    # Save variable regions
    variable_regions_df = region_names_to_coordinates(variable_regions)
    variable_regions_df.sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(region_sets_dir, "highly_variable_regions.bed"),
        sep="\t",
        header=False,
        index=False
    )
    print(f"   ✅ Saved highly variable regions to BED file")
    
    # Find differentially accessible regions by age
    markers_dict = None
    if 'age' in cistopic_obj.cell_data.columns:
        print("\n📊 Finding differentially accessible regions by age...")
        try:
            markers_dict = find_diff_features(
                cistopic_obj,
                imputed_acc_obj,
                variable='age',
                var_features=variable_regions,
                contrasts=None,
                adjpval_thr=ADJPVAL_THR,
                log2fc_thr=LOG2FC_THR,
                n_cpu=N_CPU,
                _temp_dir=None,
                split_pattern='-'
            )
            
            # Save DAR results
            print("\n💾 Saving differentially accessible regions...")
            for group_name, dar_df in markers_dict.items():
                # Save full results as CSV
                dar_df.to_csv(
                    os.path.join(out_dir, f'DARs_age_{group_name}.csv')
                )
                
                # Save as BED file
                if len(dar_df) > 0:
                    region_names_to_coordinates(
                        dar_df.index
                    ).sort_values(
                        ["Chromosome", "Start", "End"]
                    ).to_csv(
                        os.path.join(region_sets_dir, "DARs_age", f"{group_name}.bed"),
                        sep="\t",
                        header=False,
                        index=False
                    )
            
            print(f"   ✅ Saved DARs for {len(markers_dict)} age groups")
            
            # Print summary
            print("\n📊 DAR Summary by Age Group:")
            for group_name, dar_df in markers_dict.items():
                n_up = sum(dar_df['Log2FC'] > 0) if 'Log2FC' in dar_df.columns else 'N/A'
                n_down = sum(dar_df['Log2FC'] < 0) if 'Log2FC' in dar_df.columns else 'N/A'
                print(f"   {group_name}: {len(dar_df)} DARs (↑{n_up}, ↓{n_down})")
                
        except Exception as e:
            print(f"   ⚠️ Error finding DARs: {e}")
    else:
        print("   ⚠️ 'age' column not found - skipping DAR analysis")
    
    # Find DARs by celltype (if available)
    if 'celltype' in cistopic_obj.cell_data.columns:
        print("\n📊 Finding differentially accessible regions by celltype...")
        try:
            markers_dict_celltype = find_diff_features(
                cistopic_obj,
                imputed_acc_obj,
                variable='celltype',
                var_features=variable_regions,
                contrasts=None,
                adjpval_thr=ADJPVAL_THR,
                log2fc_thr=LOG2FC_THR,
                n_cpu=N_CPU,
                _temp_dir=None,
                split_pattern='-'
            )
            
            # Save DAR results
            os.makedirs(os.path.join(region_sets_dir, "DARs_celltype"), exist_ok=True)
            for group_name, dar_df in markers_dict_celltype.items():
                dar_df.to_csv(
                    os.path.join(out_dir, f'DARs_celltype_{group_name}.csv')
                )
                if len(dar_df) > 0:
                    region_names_to_coordinates(
                        dar_df.index
                    ).sort_values(
                        ["Chromosome", "Start", "End"]
                    ).to_csv(
                        os.path.join(region_sets_dir, "DARs_celltype", f"{group_name}.bed"),
                        sep="\t",
                        header=False,
                        index=False
                    )
            
            print(f"   ✅ Saved DARs for {len(markers_dict_celltype)} cell types")
            
        except Exception as e:
            print(f"   ⚠️ Error finding celltype DARs: {e}")
    
    return imputed_acc_obj, normalized_imputed_acc_obj, variable_regions, markers_dict


# =============================================================================
# STEP 7: PLOT IMPUTED FEATURES (OPTIONAL)
# =============================================================================

def step7_plot_imputed_features(cistopic_obj, imputed_acc_obj, markers_dict):
    """Plot top differentially accessible regions on UMAP."""
    print("\n" + "="*60)
    print("STEP 7: Plotting Imputed Features")
    print("="*60)
    
    if markers_dict is None:
        print("⚠️ No markers dict available - skipping feature plots")
        return
    
    # Plot top DARs for each age group
    for group_name, dar_df in markers_dict.items():
        if len(dar_df) == 0:
            continue
        
        # Get top 6 regions by absolute log2FC
        if 'Log2FC' in dar_df.columns:
            top_regions = dar_df.reindex(
                dar_df['Log2FC'].abs().sort_values(ascending=False).index
            ).head(6).index.tolist()
        else:
            top_regions = dar_df.head(6).index.tolist()
        
        if len(top_regions) > 0:
            try:
                fig = plot_imputed_features(
                    cistopic_obj,
                    reduction_name='UMAP',
                    imputed_acc_obj=imputed_acc_obj,
                    features=top_regions,
                    scale=True,
                    num_columns=3
                )
                plt.savefig(
                    os.path.join(plots_dir, f'imputed_features_{group_name}.png'),
                    dpi=150,
                    bbox_inches='tight'
                )
                plt.close()
                print(f"   ✅ Saved imputed feature plots for {group_name}")
            except Exception as e:
                print(f"   ⚠️ Could not plot features for {group_name}: {e}")


# =============================================================================
# STEP 8: SAVE FINAL OBJECT
# =============================================================================

def step8_save_final(cistopic_obj):
    """Save the final cisTopic object with all analyses."""
    print("\n" + "="*60)
    print("STEP 8: Saving Final cisTopic Object")
    print("="*60)
    
    final_path = os.path.join(out_dir, 'cisTopicObject_lda_complete.pkl')
    save_cistopic_object(cistopic_obj, final_path)
    
    # Print final summary
    print("\n" + "="*60)
    print("📊 FINAL SUMMARY")
    print("="*60)
    print(f"Cells: {len(cistopic_obj.cell_names)}")
    print(f"Regions: {len(cistopic_obj.region_names)}")
    print(f"Topics: {cistopic_obj.selected_model.n_topic}")
    print(f"Metadata columns: {cistopic_obj.cell_data.columns.tolist()}")
    
    # List output files
    print("\n📁 Output Files:")
    print(f"   Main object: {final_path}")
    print(f"   Plots: {plots_dir}/")
    print(f"   Region sets: {region_sets_dir}/")
    
    return final_path


# =============================================================================
# MAIN EXECUTION
# =============================================================================

def main(skip_steps=None):
    """
    Run the complete LDA analysis pipeline.
    
    Parameters
    ----------
    skip_steps : list, optional
        List of step numbers to skip (e.g., [3, 7] to skip clustering and feature plots)
    """
    if skip_steps is None:
        skip_steps = []
    
    print("\n" + "="*60)
    print("pycisTopic LDA ANALYSIS PIPELINE")
    print("="*60)
    
    # Step 1: Load data
    cistopic_obj = step1_load_data()
    
    # Step 2: Evaluate and select LDA model
    if 2 not in skip_steps:
        cistopic_obj, selected_model = step2_evaluate_models(cistopic_obj)
    else:
        # Load pre-existing object with LDA
        cistopic_obj = load_cistopic_object(os.path.join(out_dir, 'cisTopicObject_lda_added.pkl'))
    
    # Step 3: Clustering and visualization
    if 3 not in skip_steps:
        cistopic_obj = step3_clustering_and_visualization(cistopic_obj)
    
    # Step 4: Topic binarization
    if 4 not in skip_steps:
        region_bin_topics_top_3k, region_bin_topics_otsu, binarized_cell_topic = step4_topic_binarization(cistopic_obj)
    else:
        # Create dummy binarized_cell_topic if skipped
        binarized_cell_topic = None
    
    # Step 5: Topic QC
    if 5 not in skip_steps and binarized_cell_topic is not None:
        topic_qc_metrics, topic_annot = step5_topic_qc(cistopic_obj, binarized_cell_topic)
    
    # Step 6: Differential accessibility
    if 6 not in skip_steps:
        imputed_acc_obj, normalized_imputed_acc_obj, variable_regions, markers_dict = step6_differential_accessibility(cistopic_obj)
    else:
        markers_dict = None
        imputed_acc_obj = None
    
    # Step 7: Plot imputed features
    if 7 not in skip_steps and markers_dict is not None:
        step7_plot_imputed_features(cistopic_obj, imputed_acc_obj, markers_dict)
    
    # Step 8: Save final object
    final_path = step8_save_final(cistopic_obj)
    
    print("\n" + "="*60)
    print("✅ PIPELINE COMPLETED SUCCESSFULLY")
    print("="*60)
    
    return cistopic_obj


# =============================================================================
# ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='pycisTopic LDA analysis pipeline')
    parser.add_argument('--skip', type=int, nargs='+', default=[],
                        help='Step numbers to skip (e.g., --skip 3 7)')
    parser.add_argument('--topics', type=int, default=350,
                        help='Number of topics to select')
    parser.add_argument('--n-cpu', type=int, default=5,
                        help='Number of CPUs for parallel processing')
    
    args = parser.parse_args()
    
    # Update configuration
    SELECTED_TOPICS = args.topics
    N_CPU = args.n_cpu
    
    # Run pipeline
    cistopic_obj = main(skip_steps=args.skip)
