# Version 1.0.30

This is the first submission of `drcte` to CRAN.

## Fixing the following comments from reviewer and resubmitting

* From: Benjamin Altmann <benjamin.altmann@wu.ac.at>
* Date: 03.05.2023 19:24

  - Please put the title of your reference in your DESCRIPTION in quotes.

    * Changed description.

  - Please write TRUE and FALSE instead of T and F. (Please don't use 'T' or
    'F' as vector names. 'T' and 'F' instead of TRUE and FALSE)

    * Done. I have updated the R files 'quantile.drcte.R', 'reshape.te.R', 
    together with the related 'quantile.drcte.Rd' and 'reshape.drcte.Rd'

  - You are setting options(warn=-1) in your function. This is not allowed. 
  Please rather use suppressWarnings() if really needed.

    * I have found two instances where I was setting options(warn = -1) and, in both
    cases I have changed the code to use suppressWarnings(), instead.
    
  - Please do not modifiy the .GlobalEnv. This is not allowed by the CRAN policies.
    -> R/drm_te.R
    
    * I have updated the file R/drm_te.R and edited a line where a commad was
    evaluated in the global environment
    

## Test environments

* Local:
  - NixOS (Linux), R 4.0.4 (x86_64-pc-linux-gnu)
  - Windows 10, R 4.0.5 (x86_64-w64-mingw32/x64 (64-bit))
* travis-ci:
  - Ubuntu Linux 16.04.6 LTS (xenial) (release and devel)
* win-builder:
  - Windows Server 2008 (release (4.0.5) and devel (4.1.0 alpha))
    - 1 NOTE:
      * New submission
* R-hub builder (https://builder.r-hub.io)
  - Windows Server 2008 R2 SP1, R-devel, 32/64 bit
    - 1 NOTE:
      * New submission
  - Ubuntu Linux 20.04.1 LTS, R-release, GCC
    - 1 NOTE:
      * New submission
  - Fedora Linux, R-devel, clang, gfortran
    - 2 NOTEs:
      * New submission
      * unable to verify current time

## Reverse dependencies

There are currently no downstream dependencies for this package.

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

