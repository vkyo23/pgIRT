
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pgIRT

`pgIRT` is an **R** package that implements Item Response Theory (IRT)
model with Polya-Gamma data augmentation and EM Algorithm. The function
is available for the data containing K &gt;= 2 response category (binary
\~ K-multinomial) and having different categories across items. In
addition, the implementation includes the dynamic IRT model. The
algorithm here is based on Goplerud (2019). This package also utilizes
the parametric bootstrap method proposed by Lewis and Poole (2004) to
estimate confidence interval. 

<!-- badges: start -->
<!-- badges: end -->

## Installation

You can install the development version of pgIRT from
[GitHub](https://github.com) with:

``` r
remotes::install_github("vkyo23/pgIRT")
```
