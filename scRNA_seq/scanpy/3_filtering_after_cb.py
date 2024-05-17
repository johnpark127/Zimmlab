import os
import scanpy as sc
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


# Getting information from ids.csv file
sample_info = pd.read_csv("../ids.csv", encoding="latin-1")

# Getting GSM numbers from ids.csv and make list called "samp_name"
samp_name = sample_info["GSM"].astype(str).tolist()


# Load metrics for cellbender ambient RNA filteration
metrics = []
for file in [x for x in os.listdir('cellbender/') if x.endswith('metrics.csv')]:
    _ = pd.read_csv('cellbender/' + file,
                header = None, names = ['Metric', 'Value']).set_index('Metric').T
    _ ['File'] = file
    metrics.append(_)

metrics = pd.concat(metrics).reset_index()

metrics.hist('fraction_counts_removed')

# adatas = [sc.read_10x_h5(x) for x in os.listdir('cellbender/') if x.endswith('filtered.h5')]


def qc(adata_single):
    
    # Load adata
    adata = sc.read_10x_h5(adata_single)

    # Checking quality befor filteration
    sns.displot(adata.obs['total_counts'], bins=100, kde=False, save = samp + "_total_counts_cb.png")
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True, save= samp + "_cb.png")
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', color='pct_counts_mt', save= samp + "_scatter_cb.png")

    # Basic filtering: filtering out mitochondrial genes and cells that have gene counts over 5k
    # Keep in mind that it's already filtered out by min_genes = 200 and min_cells = 3 BEFORE the cellbender
    # Consider adding additional filtering steps such as MAD(median absolute deviations) if downstream analysis doesn't look good (check below)
    adata = adata[adata.obs.n_genes_by_counts < 5000, :]
    adata = adata[adata.obs.pct_counts_mt < 50, :]


#    # MAD definition
#    def is_outlier(adata, metric: str, nmads: int):
#        M = adata.obs[metric]
#        outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
#            np.median(M) + nmads * median_abs_deviation(M) < M
#        )
#        return outlier
#    data.obs["outlier"] = (
#        is_outlier(adata, "log1p_total_counts", 5)
#        | is_outlier(adata, "log1p_n_genes_by_counts", 5)
#        | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
#    )
#    # Check the outlier:
#    adata.obs.outlier.value_counts()
#
#    # Filtering out using MAD
#    print(f"Total number of cells: {adata.n_obs}")
#    adata = adata[(~adata.obs.outlier)].copy()
#    print(f"Number of cells after filtering of low quality cells: {adata.n_obs}")
#    # MAD is not applied to mitochondrial genes; so manually set threshold
#    adata = adata[adata.obs.pct_counts_mt < 50, :]


    # Checking quailtiy after filteration
    sns.displot(adata.obs['total_counts'], bins=100, kde=False, save = samp + "_total_counts_filtered.png")
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True, save= samp + "_filtered.png")
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', color='pct_counts_mt', save= samp + "_scatter_filtered.png")
    
    # Save filtered data
    adata.write_h5ad(f"processed{samp}_cb_filtered.h5ad")
    print(f"Filteration completed: {samp}")

# Looping qc() function for each h5 files 
for samp in samp_name:
    # Defining adatas path
    adatas = [os.path.join('cellbender/', x) for x in os.listdir('cellbender/') if x.startswith(samp) and x.endswith('filtered.h5')]

    for adata_single in adatas:
        qc(adata_single)

