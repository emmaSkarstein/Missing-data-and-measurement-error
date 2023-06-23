# Missing-data-and-measurement-error

This repository contains the files both to 

1) reproduce all analysis in the article *A joint Bayesian framework for missing data and measurement error using integrated nested Laplace approximations*, 
2) generate the website https://emmaskarstein.github.io/Missing-data-and-measurement-error/. 

The website is made using Quarto, and the compiled html-files can be found in the folder called `docs`. I recommend opening the file at `docs/index.html` to get started! 

If you prefer looking directly at/running `R`-scripts, you can clone this repository and run the `R`-scripts in the folder called `R`. These `R`-scripts are automatically generated from the Quarto documents using `knitr::purl`.

## How to reproduce figures from article

- Figure 2: Run the file `simulation_study.qmd`
- Figure 3: Run the file `missing_covariate_imputation.qmd`
- Figure 4: Run the file `missing_and_mismeasured_covariate.qmd`

## Data sets

See the respective data description files, `docs/data_bloodpressure.html`, `docs/data_nhanes.html`, `docs/data_simulated.html`.

## Session info

R-INLA version ..........: 22.05.07
Date ....................: Sat May 7 12:43:31 PM +03 2022 (Version_22.05.07)
Maintainers .............: Havard Rue <hrue@r-inla.org>
                         : Finn Lindgren <finn.lindgren@gmail.com>
                         : Elias Teixeira Krainski <elias@r-inla.org>
Main web-page ...........: www.r-inla.org
Download-page ...........: inla.r-inla-download.org
Repository ..............: github.com/hrue/r-inla
Email support ...........: help@r-inla.org
                         : r-inla-discussion-group@googlegroups.com
                         
                         

─ Session info ───────────────────────────────────────────────────────────────
 setting  value
 version  R version 4.2.1 (2022-06-23)
 os       macOS Big Sur ... 10.16
 system   x86_64, darwin17.0
 ui       X11
 language (EN)
 collate  en_US.UTF-8
 ctype    en_US.UTF-8
 tz       Europe/Oslo
 date     2023-06-23
 pandoc   3.1.1 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/ (via rmarkdown)

─ Packages ───────────────────────────────────────────────────────────────────
 package      * version  date (UTC) lib source
 class          7.3-20   2022-01-16 [1] CRAN (R 4.2.1)
 classInt       0.4-8    2022-09-29 [1] CRAN (R 4.2.0)
 cli            3.6.1    2023-03-23 [1] CRAN (R 4.2.0)
 codetools      0.2-18   2020-11-04 [1] CRAN (R 4.2.1)
 colorspace     2.1-0    2023-01-23 [1] CRAN (R 4.2.0)
 DBI            1.1.3    2022-06-18 [1] CRAN (R 4.2.0)
 e1071          1.7-12   2022-10-24 [1] CRAN (R 4.2.0)
 fansi          1.0.4    2023-01-22 [1] CRAN (R 4.2.0)
 farver         2.1.1    2022-07-06 [1] CRAN (R 4.2.0)
 foreach      * 1.5.2    2022-02-02 [1] CRAN (R 4.2.0)
 ggplot2        3.4.2    2023-04-03 [1] CRAN (R 4.2.0)
 glue           1.6.2    2022-02-24 [1] CRAN (R 4.2.0)
 gtable         0.3.3    2023-03-21 [1] CRAN (R 4.2.0)
 INLA         * 22.05.07 2022-05-07 [1] local
 inlabru        2.7.0    2022-12-02 [1] CRAN (R 4.2.0)
 isoband        0.2.7    2022-12-20 [1] CRAN (R 4.2.0)
 iterators      1.0.14   2022-02-05 [1] CRAN (R 4.2.0)
 KernSmooth     2.23-20  2021-05-03 [1] CRAN (R 4.2.1)
 labeling       0.4.2    2020-10-20 [1] CRAN (R 4.2.0)
 lattice        0.20-45  2021-09-22 [1] CRAN (R 4.2.1)
 lifecycle      1.0.3    2022-10-07 [1] CRAN (R 4.2.0)
 magrittr       2.0.3    2022-03-30 [1] CRAN (R 4.2.0)
 MASS           7.3-57   2022-04-22 [1] CRAN (R 4.2.1)
 Matrix       * 1.5-1    2022-09-13 [1] CRAN (R 4.2.0)
 MatrixModels   0.5-1    2022-09-11 [1] CRAN (R 4.2.0)
 mgcv           1.8-40   2022-03-29 [1] CRAN (R 4.2.1)
 munsell        0.5.0    2018-06-12 [1] CRAN (R 4.2.0)
 nlme           3.1-157  2022-03-25 [1] CRAN (R 4.2.1)
 patchwork      1.1.2    2022-08-19 [1] CRAN (R 4.2.0)
 pillar         1.9.0    2023-03-22 [1] CRAN (R 4.2.0)
 pkgconfig      2.0.3    2019-09-22 [1] CRAN (R 4.2.0)
 plyr           1.8.7    2022-03-24 [1] CRAN (R 4.2.0)
 proxy          0.4-27   2022-06-09 [1] CRAN (R 4.2.0)
 R6             2.5.1    2021-08-19 [1] CRAN (R 4.2.0)
 RColorBrewer   1.1-3    2022-04-03 [1] CRAN (R 4.2.0)
 Rcpp           1.0.10   2023-01-22 [1] CRAN (R 4.2.0)
 rgdal          1.5-32   2022-05-09 [1] CRAN (R 4.2.0)
 rlang          1.1.1    2023-04-28 [1] CRAN (R 4.2.0)
 s2             1.1.0    2022-07-18 [1] CRAN (R 4.2.0)
 scales         1.2.1    2022-08-20 [1] CRAN (R 4.2.0)
 sf             1.0-12   2023-03-19 [1] CRAN (R 4.2.0)
 sp           * 1.6-0    2023-01-19 [1] CRAN (R 4.2.0)
 tibble         3.2.1    2023-03-20 [1] CRAN (R 4.2.0)
 units          0.8-0    2022-02-05 [1] CRAN (R 4.2.0)
 utf8           1.2.3    2023-01-31 [1] CRAN (R 4.2.0)
 vctrs          0.6.2    2023-04-19 [1] CRAN (R 4.2.0)
 viridisLite    0.4.1    2022-08-22 [1] CRAN (R 4.2.0)
 withr          2.5.0    2022-03-03 [1] CRAN (R 4.2.0)
 wk             0.7.0    2022-10-13 [1] CRAN (R 4.2.0)

 [1] /Library/Frameworks/R.framework/Versions/4.2/Resources/library

──────────────────────────────────────────────────────────────────────────────