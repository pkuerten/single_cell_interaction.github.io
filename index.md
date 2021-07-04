Single Cell Workshop - Cell interaction analysis
================
Alexandre Mondaini and Asif Javed 


## Data download

We will use a PBMC dataset made available as part of the `SeuratData` package. This dataset contains two PBMC samples: Stimulated sample 
treated with IFN Beta as well as control samples. Please follow instructions on `SeuratData` download below. On the first attempt to access it, the data 
will be downloaded to your directory. Subsequent attempts would read it from that folder.

For interaction analysis we will be using NicheNet which requires a few curated datasets. While these can be downloaded from within R environment, 
the download speed was too slow for me. I have instead downloaded it to my dropbox. The files can be downloaded from the following links.
[gr_network.rds](https://www.dropbox.com/s/ajv5ldiv4x3r5uf/gr_network.rds?dl=0)
[ligand_target_matrix.rds](https://www.dropbox.com/s/yf2vzjre3i2qvo2/ligand_target_matrix.rds?dl=0)
[lr_network.rds](https://www.dropbox.com/s/muv638ah3n5f0f4/lr_network.rds?dl=0)
[weighted_networks.rds](https://www.dropbox.com/s/x1po2sk500ac6aj/weighted_networks.rds?dl=0)

Please copy this dataset to the working directory you intend to use for the workshop.

## Prerequisite

The workshop content uses `Seurat` <a href="https://satijalab.org/seurat/articles/integration_introduction.html">data integration vignettes</a>
to generate the input data for interaction analysis. It assumes you are familiar with the preprocessing components.

## Install packages

Let’s start by installing the necessary packages one by one

``` r
# I needed to install limma separately
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("limma")
```

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("glmGamPoi")

```

``` r
# install.packages("devtools")
devtools::install_github('satijalab/seurat-data')
```
```r
devtools::install_github("saeyslab/nichenetr")

```

``` r
install.packages("tidyverse")

```

``` r
install.packages("dplyr")
```

And loading them into R.

``` r
# load into your session
library(SeuratData)
library(Seurat)

# load dataset
LoadData("ifnb")
```

## Read in the dataset and Create a Seurat objects

The first step is to read the data from the `SeuratData` package. The input data is split based on sample and individually `SCtranform` is applied. Notice
this replaces the `NormalizeData`, `ScaleData`, and `FindVariableFeatures`. SCTransform is a statistical approach specifically designed for single cell UMI 
count data. It overcomes some of the overfitting limitations of prior bulk designed normalization methods. More details of its advantages can be found in 
[SCtranform manuscript](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1). Notice that we applied no quality control steps. This 
is because the data was precleaned obviating the need to repeat these steps.

``` r
ifnb.list <- SplitObject(ifnb, split.by = "stim")
ifnb.list <- lapply(X = ifnb.list, FUN = SCTransform, method = "glmGamPoi")
features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 3000)
ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)
ifnb.list <- lapply(X = ifnb.list, FUN = RunPCA, features = features)
```

Next we integrate the two datasets. The integration relies on highly variable features which are common to both datasets and is conducted in two steps. The first
step `FindIntegrationAnchors` defines <b> anchors</b> or cell pairs (one member from each sample) which are highly similar and hence can be confidently expected 
to be assigned to the same cluster (cell type and state). The second step `IntegrateData` uses the defined anchor pairs to align the complete datasets. 
More details on the integration method upgrades since the original Seurat paper can be found in their [more recent publication](https://pubmed.ncbi.nlm.nih.gov/34062119/)

``` r
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, normalization.method = "SCT",    anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 20)
immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT", dims = 1:30)
```

The integrated object contains both the batch effect corrected values as well as the original count values as separate assays.

``` r
immune.combined.sct <- RunPCA(immune.combined.sct, verbose = FALSE)
immune.combined.sct <- RunUMAP(immune.combined.sct, reduction = "pca", dims = 1:30)
```
``` r

# Visualization
p1 <- DimPlot(immune.combined.sct, reduction = "umap", group.by = "stim")
p2 <- DimPlot(immune.combined.sct, reduction = "umap", group.by = "seurat_annotations", label = TRUE,
    repel = TRUE)
p1 + p2
```

## Cell interaction analysis

[NicheNet](https://pubmed.ncbi.nlm.nih.gov/31819264/) aims to predicts ligand and target cell interaction links by combining cluster specific single 
cell expression data with prior knowledge of ligand-receptor pairs and gene regulatory networks downstream of the targets. In particular it aims to define
ligands which best explain the differential expression observed in target cluster. 

We begin by loading the neccessary R packages.
:Warning: Some of the `Seurat` integration commands failed to successfully execute after I loaded the NicheNet package. There might be some incompatabilities between the two packages.

``` r
library(nichenetr)
library(tidyverse)
```
