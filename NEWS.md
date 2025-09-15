News
=====

# cellGeometry 0.5.5
###### 10/09/2025
* Fix bug in `mergeMarkers()` if cellMarkers object has no cell grouping vector.

# cellGeometry 0.5.4
###### 09/09/2025
* Fix CRAN checks
* Switch from using `cat()` to `message()`

# cellGeometry 0.5.3
###### 26/08/2025
* Massive speed up (4-5x) in compensation optimisation

# cellGeometry 0.5.2
###### 30/07/2025

* Expanded to 3 methods for identifying outlier genes (var.e, Cook's distance,
Studentized residuals)
* Improved (nested) parallelisation of `tune_deconv()`
* Rewrite weights code (faster)

# cellGeometry 0.5.0
###### 14/07/2025

* Calculation of SE
* Detection and multipass removal of problematic genes with extreme residuals

# cellGeometry 0.2.0
###### 24/01/2025

* This is the initial build of cellGeometry
