---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.1'
      jupytext_version: 1.1.1
  kernelspec:
    display_name: Python [conda env:anaconda3]
    language: python
    name: conda-env-anaconda3-py
---

```python init_cell=true run_control={"read_only": false, "frozen": false}
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pandas import DataFrame, Series

import h5py
from tqdm import tqdm_notebook

%matplotlib inline

%run ../src/count_matrix_metrics.py
%run ../src/utils.py
```

## Slide-seq

```python run_control={"read_only": false, "frozen": false}
adata = load_slide_seq(data_path="../data/slide_seq_ob/")
adata = adata[(adata.obs["x"].values > 500) & (adata.obs["x"].values < 5000) & (adata.obs["y"].values > 800)]
```

```python run_control={"read_only": false, "frozen": false}
plot_expression_metrics(adata.X, xlim1=1000, xlim2=1000, xlim3=500);
```

```python run_control={"read_only": false, "frozen": false}
get_expression_value_fracs(adata.X)
```

```python run_control={"read_only": false, "frozen": false}
get_scalar_metrics(adata.X)
```

#### Visualize

```python run_control={"read_only": false, "frozen": false}
plt.scatter(adata.obs.x, adata.obs.y, c=adata.obs.mit_frac, s=0.5);
plt.title("Mirochondrial fraction");
```

```python run_control={"read_only": false, "frozen": false}
adata.obs.mit_frac.hist(bins=30);
plt.xlabel("Mitochondrial fraction");
```

```python run_control={"read_only": false, "frozen": false}
%time adata_processed = process_scanpy(adata, cl_resolution=0.3, n_neighbors=10)
```

```python
plot_clusters_spatial(adata_processed, titles=["Expression UMAP", "Spatial position"])
```

### Merge

```python run_control={"read_only": false, "frozen": false}
adata_collapsed = merge_slide_seq_beads(adata, grid_size=101)
```

```python run_control={"read_only": false, "frozen": false}
adata_collapsed.obs.mit_frac.hist(bins=30);
```

```python run_control={"read_only": false, "frozen": false}
plot_expression_metrics(adata_collapsed.X, xlim1=2000, xlim2=2000, xlim3=500);
```

```python run_control={"read_only": false, "frozen": false}
get_expression_value_fracs(adata_collapsed.X)
```

```python run_control={"read_only": false, "frozen": false}
get_scalar_metrics(adata_collapsed.X)
```

#### Visualize

```python run_control={"read_only": false, "frozen": false}
%time adata_collapsed_processed = process_scanpy(adata_collapsed, cl_resolution=0.4)
```

```python run_control={"read_only": false, "frozen": false}
plot_clusters_spatial(adata_collapsed_processed, s=10, figsize=(20, 10), titles=["Expression UMAP", "Spatial position"])
```

```python run_control={"read_only": false, "frozen": false}
plt.scatter(adata_collapsed.obs.x, adata_collapsed.obs.y, c=adata_collapsed[:,"Pbx1"].X, s=2);
# plt.title("Mirochondrial fraction");
```

<!-- #region {"run_control": {"read_only": false, "frozen": false}} -->
### Merge hard
<!-- #endregion -->

```python run_control={"read_only": false, "frozen": false}
adata_collapsed_hard = merge_slide_seq_beads(adata, grid_size=20)
```

```python run_control={"read_only": false, "frozen": false}
adata_collapsed_hard.obs["n_merged"].hist(bins=20);
```

```python run_control={"read_only": false, "frozen": false}
adata_collapsed_hard = adata_collapsed_hard[adata_collapsed_hard.obs["n_merged"].values > 25,:]
plot_expression_metrics(adata_collapsed_hard.X, xlim1=80000, xlim2=8000, xlim3=300);
```

```python run_control={"read_only": false, "frozen": false}
get_expression_value_fracs(adata_collapsed_hard.X)
```

```python run_control={"read_only": false, "frozen": false}
get_scalar_metrics(adata_collapsed_hard.X)
```

#### Visualization

```python run_control={"read_only": false, "frozen": false}
fig = plt.figure(figsize=(5, 5))
adata_collapsed_hard_processed = process_scanpy(adata_collapsed_hard, n_neighbors=10, do_log=True, n_od_genes=1000)
plot_clustering(adata_collapsed_hard_processed)
```

## SpatialTranscriptomics

```python run_control={"read_only": false, "frozen": false}
adata_st = load_spatial_transcriptomics("../data/spatial_transcriptomics_ob/Rep1_MOB_count_matrix-1.tsv")
```

```python run_control={"read_only": false, "frozen": false}
plot_expression_metrics(adata_st.X, xlim1=80000, xlim2=8000, xlim3=300, figsize=(15, 3));
```

```python run_control={"read_only": false, "frozen": false}
get_expression_value_fracs(adata_st.X)
```

```python run_control={"read_only": false, "frozen": false}
get_scalar_metrics(adata_st.X)
```

```python
adata_st_processed = process_scanpy(adata_st, n_neighbors=10, do_log=True, n_od_genes=1000)
plot_clusters_spatial(adata_st_processed, s=40, titles=["Expression UMAP", "Spatial position"])
```
