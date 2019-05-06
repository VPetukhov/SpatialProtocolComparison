import pandas as pd
from pandas import DataFrame
import h5py
import numpy as np

from tqdm import tqdm_notebook
import scanpy as sc

import matplotlib.pyplot as plt


def process_scanpy(adata: sc.AnnData, n_neighbors: int = 10, n_pcs: int = 50, metric: str = 'cosine', cl_resolution: float = 0.5):
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, metric=metric)
    sc.tl.umap(adata)
    sc.tl.louvain(adata, resolution=cl_resolution)


def plot_clustering(adata: sc.AnnData, clustering_type: str = "louvain", ax=None, **kwargs):
    if ax is None:
        ax = plt.gca()

    for cl in sorted(adata.obs[clustering_type].unique()):
        c_coord = adata.obs[adata.obs[clustering_type].values == cl]
        ax.scatter(c_coord["x"].values, c_coord["y"].values, label=str(cl), **kwargs)

    ax.legend()


def plot_clusters_spatial(adata: sc.AnnData, clustering_type: str = "louvain", figsize=(15, 6), s=0.5, **kwargs):
    fig, axes = plt.subplots(ncols=2, figsize=figsize)
    sc.pl.umap(adata, color=clustering_type, legend_loc='on data', ax=axes[0], show=False)
    plot_clustering(adata, s=s, ax=axes[1], **kwargs)


# Slide-seq

def load_slide_seq(mat_path: str = "bijectivemapping.mat", gene_path: str = "Genes.csv", location_path="Locations.csv",
                   data_path: str = None, min_cells: int = 5, min_genes: int = 5):
    if data_path is not None:
        mat_path = data_path + "/" + mat_path
        gene_path = data_path + "/" + gene_path
        location_path = data_path + "/" + location_path

    with h5py.File(mat_path, "r") as data:
        unique_mapped_dge = data["UniqueMappedDGE"][:, :]

    gene_names = pd.read_csv(gene_path, header=None).values.flatten()
    coordinates = pd.read_csv(location_path, header=None).values

    adata = sc.AnnData(DataFrame(unique_mapped_dge, columns=gene_names))
    adata.obs_names = adata.obs_names.map(str)

    mit_gene_mask = adata.var_names.map(lambda x: x[:3] == 'mt-')
    adata.obs["mit_frac"] = adata[:, mit_gene_mask].X.sum(axis=1) / adata.X.sum(axis=1)
    adata.obs["x"] = coordinates[:, 0]
    adata.obs["y"] = coordinates[:, 1]

    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    return adata


def merge_slide_seq_beads(adata: sc.AnnData, grid_size: int):
    df = adata.to_df().copy()
    df["x"] = adata.obs["x"].values
    df["y"] = adata.obs["y"].values
    x_borders, y_borders = [np.linspace(df[s].min(), df[s].max(), grid_size) for s in ["x", "y"]]

    cells_collapsed = []
    sizes = []
    x_vals, y_vals = [], []
    for xs, xe in tqdm_notebook(zip(x_borders[:-1], x_borders[1:]), total=len(x_borders) - 1):
        for ys, ye in zip(y_borders[:-1], y_borders[1:]):
            mask = (df["x"].values > xs) & (df["x"].values <= xe) & (df["y"].values > ys) & (df["y"].values <= ye)
            if mask.any():
                df_cur = df[mask]
                cells_collapsed.append(df_cur.iloc[:, :-2].sum())
                sizes.append(df_cur.shape[0])
                x_vals.append(df_cur["x"].mean())
                y_vals.append(df_cur["y"].mean())

    adata_collapsed = sc.AnnData(DataFrame(cells_collapsed))
    adata_collapsed.obs_names = adata_collapsed.obs_names.map(str)
    adata_collapsed.obs["n_merged"] = sizes
    adata_collapsed.obs["x"] = x_vals
    adata_collapsed.obs["y"] = y_vals

    mit_gene_mask = adata_collapsed.var_names.map(lambda x: x[:3] == 'mt-')
    adata_collapsed.obs["mit_frac"] = adata_collapsed[:, mit_gene_mask].X.sum(axis=1) / adata_collapsed.X.sum(axis=1)

    return adata_collapsed


# SpatialTranscriptmics

def load_spatial_transcriptomics(data_path: str, min_cells: int = 5, min_genes: int = 5):
    df = pd.read_table(data_path)
    coordinates = np.array(list(df.iloc[:, 0].map(lambda x: [float(v) for v in x.split('x')]).values))

    adata = sc.AnnData(df.iloc[:, 1:])
    adata.obs_names = adata.obs_names.map(str)

    adata.obs["x"] = coordinates[:, 0]
    adata.obs["y"] = coordinates[:, 1]

    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    return adata


# STARmap

def read_starmap(data_path: str):
    expression_df = np.load(data_path + "/cell_barcode_count.npy")
    gene_names = pd.read_csv(data_path + "/genes.csv", header=None)[0].values
    return sc.AnnData(DataFrame(expression_df, columns=gene_names))


# Seq-FISH+

def read_seq_fish_df(path, cell_id):
    with open(path) as f:
        lines = [np.array(l[:-1].split(",")) for l in f]
    
    df = DataFrame(dict(gene=lines[0], x=lines[1].astype(float), y=lines[2].astype(float)))
    df["cell"] = cell_id
    return df
