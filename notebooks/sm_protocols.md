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

```python run_control={"read_only": false, "frozen": false}
import pandas as pd
import os
from tqdm import tqdm_notebook
import scanpy as sc

%run ../src/count_matrix_metrics.py
%run ../src/utils.py
```

<!-- #region {"heading_collapsed": true} -->
### Allen smFISH
<!-- #endregion -->

```python run_control={"read_only": false, "frozen": false} hidden=true
df_spatial = pd.read_csv("../data/allen_smfish.csv")
adata = sc.AnnData(df_spatial.groupby(["gene", "cell"]).size().unstack(fill_value=0).T)
```

```python run_control={"read_only": false, "frozen": false} hidden=true
plot_expression_metrics(adata.X, xlim1=1000, xlim2=25, xlim3=5000);
```

```python run_control={"read_only": false, "frozen": false} hidden=true
get_expression_value_fracs(adata.X)
```

```python run_control={"read_only": false, "frozen": false} hidden=true
get_scalar_metrics(adata.X)
```

<!-- #region {"heading_collapsed": true} -->
### osmFISH
<!-- #endregion -->

```python run_control={"read_only": false, "frozen": false} hidden=true
adata = sc.read_loom("../data/osmFISH_SScortex_mouse_all_cells.loom")
adata.X = adata.X.A
```

```python run_control={"read_only": false, "frozen": false} hidden=true
plot_expression_metrics(adata.X, xlim1=1000, xlim2=40, xlim3=8000);
```

```python run_control={"read_only": false, "frozen": false} hidden=true
get_expression_value_fracs(adata.X)
```

```python run_control={"read_only": false, "frozen": false} hidden=true
get_scalar_metrics(adata.X)
```

<!-- #region {"heading_collapsed": true} -->
### STARmap
<!-- #endregion -->

<!-- #region {"heading_collapsed": true, "hidden": true} -->
#### 1020 genes
<!-- #endregion -->

```python run_control={"read_only": false, "frozen": false} hidden=true
adata_1020 = read_starmap("../data/star_map/visual_1020_20180505_BY3_1kgenes/")
```

```python run_control={"read_only": false, "frozen": false} hidden=true
plot_expression_metrics(adata_1020.X, xlim1=1300, xlim2=600, xlim3=1000);
```

```python run_control={"read_only": false, "frozen": false} hidden=true
get_expression_value_fracs(adata_1020.X)
```

```python run_control={"read_only": false, "frozen": false} hidden=true
get_scalar_metrics(adata_1020.X)
```

<!-- #region {"heading_collapsed": true, "hidden": true} -->
#### 160 genes
<!-- #endregion -->

```python run_control={"read_only": false, "frozen": false} hidden=true
adata_160 = read_starmap("../data/star_map/visual_160_20171120_BF4_light/")
```

```python run_control={"read_only": false, "frozen": false} hidden=true
plot_expression_metrics(adata_160.X, xlim1=1200, xlim2=160, xlim3=1000);
```

```python run_control={"read_only": false, "frozen": false} hidden=true
get_expression_value_fracs(adata_160.X)
```

```python run_control={"read_only": false, "frozen": false} hidden=true
get_scalar_metrics(adata_160.X)
```

<!-- #region {"heading_collapsed": true, "hidden": true} -->
#### 1020 genes, subset
<!-- #endregion -->

```python run_control={"read_only": false, "frozen": false} hidden=true
adata_1020s = sc.AnnData(adata_1020.to_df()[np.intersect1d(adata_160.var_names, adata_1020.var_names)])
```

```python run_control={"read_only": false, "frozen": false} hidden=true
plot_expression_metrics(adata_1020s.X, xlim1=200, xlim2=75, xlim3=1000);
```

```python run_control={"read_only": false, "frozen": false} hidden=true
get_expression_value_fracs(adata_1020s.X)
```

```python run_control={"read_only": false, "frozen": false} hidden=true
get_scalar_metrics(adata_1020s.X)
```

### seqFISH+

```python run_control={"read_only": false, "frozen": false}
df_spatial = pd.read_csv("../data/seq_fish_nih3t3.csv")
adata = sc.AnnData(df_spatial.groupby(["gene", "cell"]).size().unstack(fill_value=0).T)
```

```python run_control={"read_only": false, "frozen": false}
plot_expression_metrics(adata.X, xlim1=60000, xlim2=10000, xlim3=110);
```

```python run_control={"read_only": false, "frozen": false}
get_expression_value_fracs(adata.X)
```

```python run_control={"read_only": false, "frozen": false}
get_scalar_metrics(adata.X)
```

## MERFISH

```python
df_spatial = pd.read_csv("../data/merfish_coords_perprocessed.csv")
df_spatial = df_spatial[df_spatial.cell > 0]
adata = sc.AnnData(df_spatial.groupby(["gene", "cell"]).size().unstack(fill_value=0).T)
adata.shape
```

```python
plot_expression_metrics(adata.X, xlim1=1000, xlim2=120, xlim3=6500);
```

```python
get_expression_value_fracs(adata.X)
```

```python
get_scalar_metrics(adata.X)
```
