### This script is generating raw h5ad files which is not filtered at all for CellBender tools

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


# Defining function for reading, filtering matrix and adding metadata
def pp(matrix_path):
    adata = sc.read_10x_mtx(matrix_path, var_names = 'gene_symbols', cache = True)

    # Specifying which sample (GSM) is being processed in this loop
    sample_row = sample_info[sample_info["GSM"] == samp]

    # Adding metadata
    adata.obs['sample'] = sample_row['GSM'].values[0]
    adata.obs['group'] = sample_row['group'].values[0]
    adata.obs['control_type'] = sample_row['control_type'].values[0]
    adata.obs['type_of_injury'] = sample_row['type_of_injury'].values[0]
    adata.obs['injury_dose'] = sample_row['injury_dose'].values[0]
    adata.obs['surgery_time'] = sample_row['surgery_time'].values[0]
    adata.obs['surgery_temp'] = sample_row['surgery_temp'].values[0]
    adata.obs['nephrectomy'] = sample_row['nephrectomy'].values[0]
    adata.obs['treatment'] = sample_row['treatment'].values[0]
    adata.obs['anesthetics'] = sample_row['anesthetics'].values[0]
    adata.obs['age_at_injury'] = sample_row['age_at_injury'].values[0]
    adata.obs['days_post_injury'] = sample_row['days_post_injury'].values[0]
    adata.obs['species'] = sample_row['species'].values[0]
    adata.obs['strain'] = sample_row['strain'].values[0]
    adata.obs['sex'] = sample_row['sex'].values[0]
    adata.obs['sorting'] = sample_row['sorting'].values[0]
    adata.obs['unique_ann'] = sample_row['unique_ann'].values[0]
    adata.obs['GSE_numbers'] = sample_row['GSE_numbers'].values[0]

    print("Metadata added: ", samp, '\n')  

    # Define mitochondrial genes
    # Mitochondrial genes are annotated with 'mt-' in mouse samples (in human, 'MT-')
    adata.var['mt'] = adata.var_names.str.startswith('mt-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    # Save raw data
    adata.write_h5ad('unfiltered_for_cb/' + samp + '_raw.h5ad')

    # Saving gene counts, total counts, and mt gene counts after filtering
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True, save= samp + "_raw.png")
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', save= samp + "_total_mt.png")
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', save= samp + "_total_gene.png")

# Looping pp() function for each GSM samples
for samp in samp_name:
    # Defining matrix path
    matrix_path = os.path.join("../outs_cellranger", samp, "outs", "filtered_feature_bc_matrix")
    pp(matrix_path)

print("Filtering and preprocessing completed for all samples.")


