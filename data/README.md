## Download data


### SpatialTranscriptomics Olfactory Bulb


[Link](http://www.spatialtranscriptomicsresearch.org/datasets/doi-10-1126science-aaf2403/) to the dataset.

```python
from os import system

data_dir = "spatial_transcriptomics_ob"
system("mkdir -p " + data_dir)
for i in range(1, 13):
    system("wget http://www.spatialtranscriptomicsresearch.org/wp-content/uploads/2016/07/Rep{}_MOB_count_matrix-1.tsv -O {}/Rep{}_MOB_count_matrix-1.tsv".format(i, data_dir, i))
```

### Slide-seq

<!-- #region -->
[Slide-seq data](https://portals.broadinstitute.org/single_cell/study/SCP354/slide-seq-public) requires login to download. File: *"Puck_180430_3.tar.gz"*.

Then, extract it:
```bash
mkdir slide-seq
tar xvf Puck_180430_3.tar.gz
mv Puck_180430_3/BeadMapping_6-12_1036/bijectivemapping.mat slide-seq/
rm -rf Puck_180430_3
```

And convert necessary fields from Matlab format (was using [Matlab on O2](https://wiki.rc.hms.harvard.edu/display/O2/Using+MATLAB))

```matlab
load("./bijectivemapping.mat")
csvwrite('Locations.csv', [UniqueMappedBeads.Locations]')
fid = fopen('Genes.csv', 'w'); fprintf(fid, '%s\n', GeneNames{1:end}); fclose(fid);
```
<!-- #endregion -->

### osmFISH

```bash
wget http://linnarssonlab.org/osmFISH/osmFISH_SScortex_mouse_all_cells.loom
```

### SeqFISH+


[Data repository](https://github.com/CaiGroup/seqFISH-PLUS)


```bash
wget https://github.com/CaiGroup/seqFISH-PLUS/raw/master/sourcedata.zip
unzip sourcedata.zip
mv sourcedata seq_fish && rm sourcedata.zip
unzip seq_fish/RNA_locations_cell_cortex.zip -d seq_fish/rna_locations
```


### STARmap

[Data](https://www.starmapresources.com/data/). Requires downloading from DropBox.

```bash
mkdir star_map
# Download data there
unzip visual_160_20171120_BF4_light.zip -d visual_160_20171120_BF4_light
unzip visual_1020_20180505_BY3_1kgenes.zip -d visual_1020_20180505_BY3_1kgenes
```
