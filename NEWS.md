# drcte 1.00.61
* 2024-02-28: Corrected an error in compCDF, that provoked wrong resampling with nonparametric NPMLE fitting
* 2023-11-06: Corrected an error in 'melt_te()', that prevented the n.subject argument to provoke an error when a column is selected from a tibble. Improved the separate fitting. Now, whenever using one of the internal distributions, the 2-parameter model is fitted when the three parameter model fails and the 1-parameter model is fitted when also the 2-parameter model fails.
* 2023-01-10: Version 1.00.30 ready for upload to CRAN. Bugs edited and check successful.
* 2022-08-19: Several edits have been made, to correct some small bugs. We added the possibility of obtaining robust standard errors for parameter estimates and predictions from 'drcteList' objects.
* 2022-06-16: Functions from package 'kedd' have been included in 'drcte', to remove the dependency on 'keed', that is no longer available at CRAN for R 4.2.0 (vers. 1.0.10)
* 2022-07-04: A bug has been corrected, that, in some instances, scrambled the namings of curveid levels with KDE and NPMLE fits (vers. 1.0.11)
