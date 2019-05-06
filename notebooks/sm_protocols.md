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

### Allen smFISH

```python run_control={"read_only": false, "frozen": false}
df_spatial = pd.read_csv("../data/allen_smfish.csv")
adata = sc.AnnData(df_spatial.groupby(["gene", "cell"]).size().unstack(fill_value=0).T)
```

```python run_control={"read_only": false, "frozen": false}
plot_expression_metrics(adata.X, xlim1=1000, xlim2=25, xlim3=5000);
```

```python run_control={"read_only": false, "frozen": false}
get_expression_value_fracs(adata.X)
```

```python run_control={"read_only": false, "frozen": false}
get_scalar_metrics(adata.X)
```

### osmFISH

```python run_control={"read_only": false, "frozen": false}
adata = sc.read_loom("../data/osmFISH_SScortex_mouse_all_cells.loom")
adata.X = adata.X.A
```

```python run_control={"read_only": false, "frozen": false}
plot_expression_metrics(adata.X, xlim1=1000, xlim2=40, xlim3=8000);
```

```python run_control={"read_only": false, "frozen": false}
get_expression_value_fracs(adata.X)
```

```python run_control={"read_only": false, "frozen": false}
get_scalar_metrics(adata.X)
```

### STARmap


#### 1020 genes

```python run_control={"read_only": false, "frozen": false}
adata_1020 = read_starmap("../data/star_map/visual_1020_20180505_BY3_1kgenes/")
```

```python run_control={"read_only": false, "frozen": false}
plot_expression_metrics(adata_1020.X, xlim1=1300, xlim2=600, xlim3=1000);
```

```python run_control={"read_only": false, "frozen": false}
get_expression_value_fracs(adata_1020.X)
```

```python run_control={"read_only": false, "frozen": false}
get_scalar_metrics(adata_1020.X)
```

#### 160 genes

```python run_control={"read_only": false, "frozen": false}
adata_160 = read_starmap("../data/star_map/visual_160_20171120_BF4_light/")
```

```python run_control={"read_only": false, "frozen": false}
plot_expression_metrics(adata_160.X, xlim1=1200, xlim2=160, xlim3=1000);
```

```python run_control={"read_only": false, "frozen": false}
get_expression_value_fracs(adata_160.X)
```

```python run_control={"read_only": false, "frozen": false}
get_scalar_metrics(adata_160.X)
```

#### 1020 genes, subset

```python run_control={"read_only": false, "frozen": false}
adata_1020s = sc.AnnData(adata_1020.to_df()[np.intersect1d(adata_160.var_names, adata_1020.var_names)])
```

```python run_control={"read_only": false, "frozen": false}
plot_expression_metrics(adata_1020s.X, xlim1=200, xlim2=75, xlim3=1000);
```

```python run_control={"read_only": false, "frozen": false}
get_expression_value_fracs(adata_1020s.X)
```

```python run_control={"read_only": false, "frozen": false}
get_scalar_metrics(adata_1020s.X)
```

### SeqFish+

```python run_control={"read_only": false, "frozen": false}
data_dir = "../data/seq_fish/rna_locations/"
%time dfs = [read_seq_fish_df(data_dir + "/" + f, i) for i,f in enumerate(os.listdir(data_dir))]
df_spatial = pd.concat(dfs, ignore_index=True)
adata = sc.AnnData(df_spatial.groupby(["gene", "cell"]).size().unstack(fill_value=0).T)
```

```python run_control={"read_only": false, "frozen": false}
plot_expression_metrics(adata.X, xlim1=15000, xlim2=6000, xlim3=1000);
```

```python run_control={"read_only": false, "frozen": false}
get_expression_value_fracs(adata.X)
```

```python run_control={"read_only": false, "frozen": false}
get_scalar_metrics(adata.X)
```
