## R CMD check results

0 errors | 0 warnings | 0 note

* There are no references/citations available for this novel method (it is 
unpublished). But we are going to submit the maths of the method as a new 
manuscript soon.

* All instances of `cat()` have been changed to `message()` throughout. All 
requested functions have been fixed.

* The only exception is the function `diagnose()` which is designed to print 
diagnostics, similar to `summary()`, so this retains `cat()`.

* Nearly all functions have a `verbose` argument which can also be used to 
suppress messages.

* The function `txtProgressBar2()` was redundant and has been removed.

* This is a resubmission of a new release.
