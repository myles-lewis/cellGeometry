# Geometric single cell deconvolution

Ultrafast deconvolution of bulk RNA-Seq datasets using a single-cell RNA-Seq
reference dataset in which cell clusters have been defined.

### Installation

Bioconductor version >=3.20 must be installed first for this package to install
correctly. For full package functionality, particularly with sparse matrices
stored on disc in the h5ad format, we recommend that the Bioconductor packages
zellkonverter, rhdf5 and HDF5Array must also be installed to be able to read
h5ad files. If you are using Seurat, then it needs to be installed. We also
recommend installing AnnotationHub to enable conversion of ensembl gene ids to
symbols.

```
# Bioconductor must be installed +/- updated first
BiocManager::install(version = "3.20")

# minimum necessary Bioconductor packages to install cellGeometry package
BiocManager::install(c("ensembldb", "DelayedArray"))

# packages needed to read h5ad files
BiocManager::install(c("zellkonverter", "rhdf5", "HDF5Array"))

# optional, if you are using Seurat
install.packages("Seurat")

# package needed to convert ensembl gene ids to symbols
BiocManager::install("AnnotationHub")
```

Install from Github
```
devtools::install_github("myles-lewis/cellGeometry")
```

### Algorithm

The algorithm is performed in two stages:

1. Optimal gene markers for each cell subclass are identified. In this part,
each gene is considered as a vector in high dimensions with cell clusters as
dimensions.

2. The bulk RNA-Seq is deconvoluted by calculating the vector projection of each
bulk RNA-Seq sample against a vector representing each cell cluster in high
dimensional gene marker space using the vector dot product. In order to adjust
for spillover in the vector projection between cell clusters, a compensation
matrix is applied.

### Example h5ad file

The following example is based on a the Cell Typist dataset (Global) which
contains 329,762 immune cells and is available on the CZ cellxgene repository
here:
https://cellxgene.cziscience.com/collections/62ef75e4-cbea-454e-a0ce-998ec40223d3

The h5ad file (2.9 Gb) for the example can be downloaded from CZ cellxgene
repository directly using this link:
https://datasets.cellxgene.cziscience.com/2ac906a5-9725-4258-8e36-21a9f6c0302a.h5ad

First we load the file in HDF5 format so that the full data remains on disc and
only subsets of the data are loaded/processed when necessary using the HDF5Array
and DelayedArray packages.

```
library(zellkonverter)
library(SingleCellExperiment)
library(cellGeometry)

typist_h5 <- readH5AD("2ac906a5-9725-4258-8e36-21a9f6c0302a.h5ad",
                      use_hdf5 = TRUE, reader = "R")
```

We extract the main count matrix and cell metadata. cellGeometry needs rownames
on the count matrix.

```
mat <- typist_h5@assays@data$X
rownames(mat) <- rownames(typist_h5)
meta <- typist_h5@colData@listData
```

### Example Seurat file

Some users report difficulties with installing zellkonverter. cellGeometry can
also be used with Seurat files although these become progressively slower with
larger datasets as well as needing substantial amounts of RAM, so for datasets >1M
cells we recommend persevering with zellkonverter and the h5ad format since it
is much faster. We include example code for loading a Seurat file below as an
alternative to h5ad.

At time of writing the rds file (2.9 Gb) in Seurat format can be downloaded from
CZ cellxgene repository directly using this link:
https://datasets.cellxgene.cziscience.com/2ac906a5-9725-4258-8e36-21a9f6c0302a.rds

CZ cellxgene state that Seurat support will end after Dec 2024.

```
library(Seurat)
typist <- readRDS("08f58b32-a01b-4300-8ebc-2b93c18f26f7.rds")  # 15.5 GB in memory

mat <- typist@assays$RNA$counts
meta <- typist@meta.data
```

### Extract cell subclasses and clusters

We first check cell cluster subclasses. Then we extract a vector which contains
the subclass cluster for each cell and a 2nd vector for broader cell groups. We
restrict the dataset to blood so that we can deconvolute blood bulk
RNA-Seq data later (this is optional).

```
table(meta$Majority_voting_CellTypist)

subcl <- meta$Majority_voting_CellTypist
cellgrp <- meta$Majority_voting_CellTypist_high

# reduce dataset to only blood (optional)
subcl[meta$tissue != "blood"] <- NA
cellgrp[meta$tissue != "blood"] <- NA
```

We then run the 1st stage of cellGeometry which generates mean gene expression
for each cell cluster (this is the slowest part). Then the best cell cluster and
cell group gene markers are identified.

```
mk <- cellMarkers(mat, subclass = subcl, cellgroup = cellgrp,
                  dual_mean = TRUE, cores = 2)
```

The `dual_mean` argument only needs to be set for the purpose of the simulation
later. Most users do not need to set this. It calculates both the standard mean
gene expression, which is mean(log2(counts +1)), as well as the arithmetic mean
of the (unlogged) counts.

The derivation of mean gene expression for each cluster and cell group is the
slowest part. If you are on linux or mac, this can be sped up using
parallelisation by setting `cores = 2` or more. Note that this can increase
memory requirements dramatically unless HFD5 is used. For this particular
dataset which is moderate in size, we find significant speed up with 4-8 cores
(64 Gb machine). For very large datasets (>1M cells) if the sc data is kept on
disc via HFD5 then many cores can be used. But if the data or subsets of it have
to be loaded into memory then we typically apportion around 16 Gb per core (e.g.
3 cores on a 64 Gb machine). So the limit on cores depends on the size of the
single-cell data, available RAM and whether HFD5 is used.

Windows users can invoke parallelisation using the future backend and setting up
a multisession plan.

We have not specified a bulk RNA-Seq dataset at this stage as this example is
based on simulation alone. However, if you have a bulk RNA-Seq dataset it is
helpful to specify it during the first call to `cellMarkers()`. It is only used
for its rownames to identify genes that overlap between the 2 datasets. The
marker signature can be updated later for different bulk datasets using
`updateMarkers()` (see below).

We convert the ensembl ids in the cellMarkers object using the built-in function
`gene2symbol()`. This needs an ensembl database to be loaded.

```
library(AnnotationHub)
ah <- AnnotationHub()
ensDb_v110 <- ah[["AH113665"]]
mk <- gene2symbol(mk, ensDb_v110)
```

The signature gene matrix can be displayed as follows.

```
signature_heatmap(mk)
```

The spillover heatmap between cell clusters can also be visualised.

```
spillover_heatmap(mk)
```

This heatmap as well as the signature heatmap reveals that some cell subclasses
'spillover' too strongly into other cell subclasses. In other words some cell
types are too similar - perhaps one is really a closely related subset of the
other. Here we see that Helper T cells are the most affected and their signature
is similar to Tcm/Naive helper T cells.

Below we update the cellMarkers object to remove 2 cell clusters which overlap
with other cell clusters and are therefore likely to be difficult to deconvolute
well if applied to real world bulk RNA-Seq. For the simulation it does not
matter whether these are removed or not.

```
mk <- updateMarkers(mk,
                    remove_subclass = c("Helper T cells", "Cytotoxic T cells"))
```

### Simulated pseudo-bulk RNA-Seq

We can generate pseudo-bulk to test the deconvolution using the following
commands. Here `generate_samples()` makes 25 samples with random cell counts,
`sim_counts`. The simulate_bulk() function operates in 2 modes. In the first
mode, the average gene expression for each cell cluster is extracted from the
cellMarkers object and used to generate the pseudo-bulk totals. In the 2nd
mode (see below) the original single-cell count data is sampled.

```
# simulated bulk
set.seed(3)
sim_counts <- generate_samples(mk, 25)
sim_percent <- sim_counts / rowSums(sim_counts) * 100
sim_pseudo <- simulate_bulk(mk, sim_counts)
```

Deconvolution itself is performed as a 2nd function `deconvolute()`. The
`plot_set()` function can be used to plot the results. The `metric_set()`
function generates a table of results.

```
# mode 1: (perfect deconvolution)
fit <- deconvolute(mk, sim_pseudo,
                   count_space = TRUE, convert_bulk = FALSE, use_filter = FALSE)
plot_set(sim_counts, fit$subclass$output)
plot_set(sim_percent, fit$subclass$percent)

metric_set(sim_percent, fit$subclass$percent)  # table of results
```

In the 2nd mode, the original scRNA-Seq count dataset is sampled. Here we
oversample the actual cell counts in `sim_counts` by 3x by setting `times = 3`.
Cells are sampled with replacement. The desired cell counts are simply
multiplied by `times` prior to sampling. Users will find that increasing `times`
from 1 to 30 or more improves the deconvolution as the sum of the gene counts
per sampled cell approaches the arithmetic mean of gene counts for each cell
cluster.

```
# mode 2: sample from original sc count matrix
sim_sampled <- simulate_bulk(mat, sim_counts, subcl, times = 3)

# fix rownames
rownames(sim_sampled) <- gene2symbol(rownames(sim_sampled), ensDb_v110)

# near optimal deconvolution of counts sampled from the original scRNA-Seq
fit2 <- deconvolute(mk, sim_sampled,
                    count_space = TRUE, convert_bulk = FALSE, use_filter = FALSE,
                    arith_mean = TRUE)

# plot results
plot_set(sim_counts, fit2$subclass$output / 3)  # adjust for 3x oversampling
plot_set(sim_percent, fit2$subclass$percent)

metric_set(sim_percent, fit2$subclass$percent)
```

Note that these settings are mathematically ideal for simulated bulk data. In
reality, we expect the scRNA-Seq signature to differ from real-world bulk
RNA-Seq due to differences in chemistry and the amplification step required by
single-cell sequencing. So we recommend the default settings for real-world
bulk data.

The marker object `mk` can be rapidly updated with new settings, e.g. to alter
the number of genes used per subclass, using the function `updateMarkers()`. If
some signature genes are missing from the bulk data, `deconvolute()` will stop
and warn you these genes are missing. `updateMarkers()` can then be used to
refine the gene signatures using only genes which are also found in the bulk
RNA-Seq dataset.

There is also a powerful function `tune_deconv()` which allows users to tune any
of the parameters available in `updateMarkers()` based on a bulk reference
dataset. The simulated pseudo-bulk data can be used for this purpose, but real
bulk RNA-Seq would be better (more realistic and better for tuning).

Also, 2 scRNA-Seq datasets can be merged using the function `mergeMarkers()`.
This merges the `cellMarkers` objects derived from each single cell dataset. One
dataset is defined as reference, and the 2nd dataset is merged into it after
adjustment for its overall distribution based on quantile mapping.
