# Geometric single cell deconvolution

Tools for identifying best gene markers for single-cell and bulk RNA-Seq datasets.

### Installation

Install from Github
```
devtools::install_github("myles-lewis/cellGeometry")
```

### BLAS libraries

Matrix operations can be substantially sped up by using an optimised BLAS 
library. The standard library installed with R is accurate but not optimised for 
speed. We suggest installing [OpenBLAS](https://www.openblas.net) on 
linux/windows and using Apple's vecLib on mac. CRAN gives instructions for 
switching to vecLib [here](https://cran.r-project.org/bin/macosx/RMacOSX-FAQ.html#Which-BLAS-is-used-and-how-can-it-be-changed_003f)
in the R for macOS FAQ. There is additional information in the R for windows FAQ
[here](https://cran.r-project.org/bin/windows/base/rw-FAQ.html#Can-I-use-a-fast-BLAS_003f).

