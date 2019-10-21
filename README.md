
<!-- README.md is generated from README.Rmd. Please edit that file -->
entbal - An alternative implementation of entropy balancing weights
===================================================================

The goal of entbal is to create an easy to use implementation of the entropy balancing algorithm outlined in Hainmueller (2012) for applied researchers to use.

Hainmueller, Jens. "Entropy balancing for causal effects: A multivariate reweighting method to produce balanced samples in observational studies." Political Analysis 20.1 (2012): 25-46.

Installation
------------

You can install entbal from github with:

``` r
# install.packages("devtools")
devtools::install_github("bvegetabile/entbal")
```

Example
-------

This is a basic example which shows you how to solve a common problem:

``` r
library(entbal)
library(survey, quietly = T, warn.conflicts = F, verbose = F)
#> Warning: package 'survey' was built under R version 3.5.2
dset <- dw99cps1
estwts <- entbal(TA ~ age + education + black + hispanic + married + nodegree + RE74 + RE75,
                      data = dset,
                      n_moments = 3,
                      verbose = 1,
                      estimand = "ATT")
#> iter   10 value -3.307619
#> iter   20 value -3.463395
#> iter   30 value -3.470933
#> iter   40 value -3.485851
#> iter   50 value -3.581389
#> iter   60 value -3.595018
#> iter   70 value -3.611668
#> iter   80 value -3.612052
#> iter   90 value -3.612056
#> final  value -3.612062 
#> converged
summary(estwts)
#> ------------------------------------------------------------------------------------------
#> Unweighted Balance Statistics:
#> ------------------------------------------------------------------------------------------
#>           MeanGroup1 SEGroup1 MeanGroup0 SEGroup0 StdDiffMeans LogRatioSE
#> age            25.82     7.14      33.23    11.04        -0.80      -0.44
#> education      10.35     2.01      12.03     2.87        -0.68      -0.36
#> black           0.84     0.36       0.07     0.26         2.43       0.13
#> hispanic        0.06     0.24       0.07     0.26        -0.05       0.06
#> married         0.19     0.39       0.71     0.45        -1.23       0.15
#> nodegree        0.71     0.45       0.30     0.46         0.91       0.21
#> RE74         2095.57  4873.40   14016.80  9569.50        -1.57      -0.67
#> RE75         1532.06  3210.54   13650.80  9270.11        -1.75      -1.06
#>           MaxKS
#> age        0.34
#> education  0.41
#> black      0.77
#> hispanic   0.01
#> married    0.52
#> nodegree   0.41
#> RE74       0.60
#> RE75       0.65
#> ------------------------------------------------------------------------------------------
#> Weighted Balance Statistics:
#> ------------------------------------------------------------------------------------------
#>           MeanGroup1 SEGroup1 MeanGroup0 SEGroup0 StdDiffMeans LogRatioSE
#> age            25.82     7.14      25.82     7.14            0          0
#> education      10.35     2.01      10.35     2.01            0          0
#> black           0.84     0.36       0.84     0.36            0          0
#> hispanic        0.06     0.24       0.06     0.24            0          0
#> married         0.19     0.39       0.19     0.39            0          0
#> nodegree        0.71     0.45       0.71     0.45            0          0
#> RE74         2095.57  4873.40    2095.57  4873.40            0          0
#> RE75         1532.06  3210.54    1532.05  3210.48            0          0
#>           MaxKS
#> age        0.07
#> education  0.04
#> black      0.00
#> hispanic   0.00
#> married    0.00
#> nodegree   0.00
#> RE74       0.16
#> RE75       0.10
#> ------------------------------------------------------------------------------------------
#> TA: 1,   Original N = 185
#>        Weighted ESS = 185
#> TA: 0,   Original N = 15992
#>        Weighted ESS = 116.46
#> ------------------------------------------------------------------------------------------
dset$wts <- estwts$wts
design <- svydesign(ids=~1, weights=~wts, data = dset)
resp <- svyglm('RE78 ~ TA', design = design)
summary(resp)
#> 
#> Call:
#> svyglm(formula = "RE78 ~ TA", design = design)
#> 
#> Survey design:
#> svydesign(ids = ~1, weights = ~wts, data = dset)
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)   4957.0      434.9  11.397   <2e-16 ***
#> TA            1392.2      722.5   1.927    0.054 .  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for gaussian family taken to be 61441391)
#> 
#> Number of Fisher Scoring iterations: 2
```
