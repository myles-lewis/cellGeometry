# Geometric single cell deconvolution

Deconvolution of bulk RNA-Seq datasets using a single-cell RNA-Seq reference 
dataset in which cell clusters have been defined.

### Installation

Bioconductor version >=3.20 must be installed first for this package to install
correctly. For full package functionality, particularly with sparse matrices
stored on disc in the h5ad format, we recommend that the Bioconductor packages
zellkonverter, rhdf5 and HDF5Array must also be installed to be able to read
h5ad files. We also recommend installing AnnotationHub to enable conversion of
ensembl gene ids to symbols.

```
# Bioconductor must be installed +/- updated first
BiocManager::install(version = "3.20")

# minimum necessary Bioconductor packages to install cellGeometry package
BiocManager::install(c("ensembldb", "DelayedArray"))

# packages needed to read h5ad files
BiocManager::install(c("zellkonverter", "rhdf5", "HDF5Array"))

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

### Example

The following example is based on a the Cell Typist dataset which is available
on the CZI cellxgene repository here:
https://cellxgene.cziscience.com/collections/62ef75e4-cbea-454e-a0ce-998ec40223d3

The h5ad file for the example can be downloaded from CZI cellxgene repository
directly using this link:
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

We first check cell cluster subclasses. Then we extract a vector which contains
the subclass cluster for each cell and a 2nd vector for broader cell groups. We
restrict the dataset to blood so that we can deconvolute blood bulk
RNA-Seq data later (this is optional).

```
table(meta$Majority_voting_CellTypist)

subcl <- meta$Majority_voting_CellTypist
subcl[meta$tissue != "blood"] <- NA
cellgrp <- meta$Majority_voting_CellTypist_high
cellgrp[meta$tissue != "blood"] <- NA
```

We then run the 1st stage of cellGeometry which generates mean gene expression
for each cell cluster (this is the slowest part). Then the best cell cluster and
cell group gene markers are identified.

```
mk <- cellMarkers(mat, subclass = subcl, cellgroup = cellgrp,
                  remove_subclass = c("Helper T cells", "Cytotoxic T cells",
                                      "Cycling T cells"),
                  dual_mean = TRUE)
```

Here we have removed 3 cell clusters which overlap too much with other cell
clusters and are therefore difficult to deconvolute.

The `dual_mean` argument only needs to be set for the purpose of the simulation
later. Most users do not need to set this. It calculates both the standard mean
gene expression, which is mean(log2(counts +1)), as well as the arithmetic mean
of the (unlogged) counts.

The derivation of mean gene expression for each cluster and cell group is the
slowest part. If you are on linux or mac, this can be sped up using
parallelisation by setting `cores = 2` or more. Note that this can increase
memory requirements dramatically. We typically use 3 cores on a 64 Gb machine,
but the limit on cores depends on the size of the single-cell data.

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

# fix rownames
rownames(sim_pseudo) <- gene2symbol(rownames(sim_pseudo), ensDb_v110)
```

Deconvolution itself is performed as a 2nd function `deconvolute()`. The
`plot_set()` function can be used to plot the results. The `metric_set()`
function generates a table of results.

```
# mode 1: (perfect deconvolution)
fit <- deconvolute(mk, sim_pseudo,
                   exp_signature = TRUE, convert_bulk = FALSE, use_filter = FALSE)
plot_set(sim_counts, fit$subclass$output)
plot_set(sim_percent, fit$subclass$percent)

metric_set(sim_percent, fit$subclass$percent)  # table of results
```

In the 2nd mode, the original scRNA-Seq count dataset is sampled. Here we
oversample the actual cell counts in `sim_counts` by 3x. Users can see that if
times is varied from 1 to 30 or more that the deconvolution improves as the
sampling approaches the arithmetic mean of the counts for each cluster.

```
# mode 2: sample from original sc count matrix
sim_sampled <- simulate_bulk(mat, sim_counts, subcl, times = 3)

# fix rownames
rownames(sim_sampled) <- gene2symbol(rownames(sim_sampled), ensDb_v110)

# near optimal deconvolution of counts sampled from the original scRNA-Seq
fit2 <- deconvolute(mk, sim_sampled,
                    exp_signature = TRUE, convert_bulk = FALSE, use_filter = FALSE,
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
