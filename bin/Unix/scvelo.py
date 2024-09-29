#!/usr/bin/env python

import sys
sys.path.insert(0,'/home/jibzhang/miniconda3/envs/singlecell/lib/python3.9/site-packages')
import anndata
import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib as plt
import cellrank as cr
import scanpy as sc

early_velo=sys.argv[1]
late_velo=sys.argv[2]
cellID=sys.argv[3]
umap=sys.argv[4]
cellclust=sys.argv[5]
out=sys.argv[6]
veloc_sample=sys.argv[7]
veloc_cluster=sys.argv[8]

D40_90 = anndata.read_loom(early_velo)
LTC = anndata.read_loom(late_velo)
D40_90.var_names_make_unique()
LTC.var_names_make_unique()
cellID_obs = pd.read_csv(cellID)
umap_cord = pd.read_csv(umap)
cell_clusters = pd.read_csv(cellclust)

cellID_obs_early = cellID_obs[cellID_obs['x'].str.contains("D40_90_")]
cellID_obs_LTC = cellID_obs[cellID_obs['x'].str.contains("LTC_")]
D40_90.obs.index=D40_90.obs.index.str.replace("37191:","D40_90_")
D40_90.obs.index=D40_90.obs.index.str.replace("x","-1")
LTC.obs.index=LTC.obs.index.str.replace("39550:","LTC_")
LTC.obs.index=LTC.obs.index.str.replace("x","-1")

sample_index = pd.DataFrame(adata.obs.index)
cluster_ordered = sample_index.merge(cell_clusters, on = "CellID")
umap_ordered = sample_index.merge(umap_cord, on="CellID")

umap = umap_ordered.iloc[:,1:3]
cluster = cluster_ordered.iloc[:,1]
color = cluster_ordered.iloc[:,2]
sample_name = umap_ordered.iloc[:,3]
adata.obsm['X_umap'] = umap.values
adata.obs['Sample_name'] = sample_name.values
adata.obs['Clusters'] = cluster.values
adata.uns['Seurat_Cluster'] = cluster.values
adata.uns['Seurat_Cluster_colors'] = color.values
adata.write(out,compression='gzip')

scv.pp.filter_and_normalize(adata, min_shared_counts=3, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=20, n_neighbors=30)
scv.tl.velocity(adata,mode="stochastic")
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_stream(adata, basis='umap', color='Sample_name',legend_loc='right margin',
                                 palette=['blue','red','yellow','green','purple'], save=veloc_sample)
scv.pl.velocity_embedding_stream(adata, basis='umap', color='Clusters', legend_loc='right margin',
                                 palette=['#F8766D','#CD9600','#7CAE00','#00BE67','#00BFC4','#00A9FF','#C77CFF','#FF61CC'],
                                 save=veloc_cluster)
