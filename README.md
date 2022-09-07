
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pgIRT

<!-- badges: start -->
<!-- badges: end -->

`pgIRT` is an **R** package that implements Item Response Theory (IRT)
model with Polya-Gamma data augmentation and EM Algorithm.
Implementation includes binary, multinomial along with dynamic IRT
models. The algorithm here is based on the procedure proposed by
Goplerud (2019). This package also utilizes the parametric bootstrap
method proposed by Lewis and Poole (2004) to estimate statistical
uncertainty.

## Installation

You can install the development version of pgIRT from
[GitHub](https://github.com) with:

``` r
remotes::install_github("vkyo23/pgIRT")
```

## Usage

### Binary IRT

Using `Senate` data from `MCMCpack`.

``` r
library(pgIRT)
require(dplyr)
require(tidyr)

data(Senate, package = "MCMCpack")
mat <- as.matrix(Senate[, 6:ncol(Senate)])
colnames(mat) <- 1:ncol(mat)
rownames(mat) <- 1:nrow(mat)

# `make_init` is an auxiliary function which computes initial values for `pgIRT`.
init <- make_init(mat, model = "bin")
fit <- pgIRT(mat, 
             model = "bin",
             init = init,
             verbose = 20)
## =========================================================================
## Polya-Gamma data augmentation Item Response Theory Model via EM Algorithm
## Model = Binomial 
## =========================================================================
## Iteration 20 : eval = alpha 2.480436e-05 elapsed 0.1 sec
## Iteration 40 : eval = alpha 4.468147e-06 elapsed 0.2 sec
## Iteration 60 : eval = alpha 1.615617e-06 elapsed 0.2 sec
## Model converged at iteration 72 : 0.3 sec

summary(fit, parameter = "theta")
## =============== Parameter = theta =============== 
## # A tibble: 102 x 3
##    variable unit_id estimate
##    <chr>      <dbl>    <dbl>
##  1 theta          1    1.84 
##  2 theta          2    1.48 
##  3 theta          3    1.80 
##  4 theta          4    1.02 
##  5 theta          5    2.05 
##  6 theta          6    0.833
##  7 theta          7    1.84 
##  8 theta          8   -1.04 
##  9 theta          9   -3.28 
## 10 theta         10   -1.35 
## # ... with 92 more rows
```

### Dynamic binary IRT

Using `Rehnquist` data from `MCMCpack`. To introduce some auxiliary
function of `pgIRT`, I fisrt construct long-format dataframe.

``` r
# Convert into long-format dataframe
data(Rehnquist, package = "MCMCpack")
long_df <- Rehnquist %>% 
  mutate(rcid = row_number()) %>% 
  select(-term) %>% 
  mutate(time = time - min(time)) %>% 
  pivot_longer(-c(time, rcid))
jname <- tibble(name = unique(long_df$name)) %>% 
  mutate(judge_id = row_number())
long_df <- long_df %>% 
  left_join(jname, by = "name")
head(long_df) %>% 
  knitr::kable()
```

| time | rcid | name      | value | judge\_id |
|-----:|-----:|:----------|------:|----------:|
|    0 |    1 | Rehnquist |     0 |         1 |
|    0 |    1 | Stevens   |     1 |         2 |
|    0 |    1 | O.Connor  |     0 |         3 |
|    0 |    1 | Scalia    |     0 |         4 |
|    0 |    1 | Kennedy   |     1 |         5 |
|    0 |    1 | Souter    |     1 |         6 |

``` r
# Generating roll-call matrix from long-format dataframe
mat_dyn <- make_rollcall(long_df, 
                         unit_id = "judge_id", 
                         bill_id = "rcid",
                         vote_col = "value")

# Initial values
constraint_dyn <- jname$judge_id[jname$name == "Thomas"]
init_dyn <- make_init(mat_dyn,
                      model = "bin_dyn",
                      T = length(unique(long_df$time)),
                      constraint = constraint_dyn)

# Generating options for dynamic estimation
dyn_ops <- make_dyn_options(long_df,
                            unit_id = "judge_id",
                            bill_id = "rcid",
                            time_id = "time",
                            vote_col = "value")

# IRT
fit_dyn <- pgIRT(mat_dyn,
             model = "bin_dyn",
             init = init_dyn,
             constraint = constraint_dyn,
             dyn_options = dyn_ops,
             verbose = 20)
## =========================================================================
## Polya-Gamma data augmentation Item Response Theory Model via EM Algorithm
## Model = Dynamic Binomial 
## =========================================================================
## Iteration 20 : eval = beta 1.025311e-05 elapsed 0 sec
## Model converged at iteration 34 : 0 sec
```

To compute confidence interval, run bootstrap via `pgIRT_boot`:

``` r
# Bootstrap
boot <- pgIRT_boot(fit_dyn, boot = 100, verbose = 20)
## ================================================================
## Parametric Bootstrap for pgIRT ( Dynamic Binomial )
## ================================================================
## Boostrap 20 DONE : 0.7 sec
## Boostrap 40 DONE : 1.6 sec
## Boostrap 60 DONE : 2.3 sec
## Boostrap 80 DONE : 3 sec
## Boostrap 100 DONE : 3.8 sec

summary(boot, parameter = "theta", ci = .95)
## ==================== Parameter = theta ==================== 
## # A tibble: 99 x 7
##    variable unit_id session ci    estimate   lwr   upr
##    <chr>      <dbl>   <int> <chr>    <dbl> <dbl> <dbl>
##  1 theta          1       1 95%      0.987 0.837  1.13
##  2 theta          1       2 95%      0.973 0.846  1.11
##  3 theta          1       3 95%      0.943 0.807  1.08
##  4 theta          1       4 95%      0.928 0.783  1.07
##  5 theta          1       5 95%      0.962 0.824  1.10
##  6 theta          1       6 95%      0.947 0.821  1.08
##  7 theta          1       7 95%      0.957 0.826  1.10
##  8 theta          1       8 95%      0.927 0.788  1.07
##  9 theta          1       9 95%      0.902 0.752  1.05
## 10 theta          1      10 95%      0.934 0.797  1.08
## # ... with 89 more rows
```

### Multinomial IRT

Using a simulated voting data `m_data`.

``` r
data("m_data")

init_mlt <- make_init(m_data,
                      model = "multi")
fit_mlt <- pgIRT(m_data,
                 model = "multi",
                 init = init_mlt, 
                 constraint = 1,
                 verbose = 20)
## =========================================================================
## Polya-Gamma data augmentation Item Response Theory Model via EM Algorithm
## Model = Multinomial (Stick-Breaking) 
## =========================================================================
## Iteration 20 : eval = beta2 0.000277376 elapsed 0 sec
## Iteration 40 : eval = beta2 9.464876e-06 elapsed 0 sec
## Model converged at iteration 55 : 0 sec

summary(fit_mlt, parameter = c("alpha", "beta"))
## =============== Parameter = alpha =============== 
## # A tibble: 30 x 3
##    variable bill_id estimate
##    <chr>      <dbl>    <dbl>
##  1 alpha1         1   -1.93 
##  2 alpha1         2    1.10 
##  3 alpha1         3    0.627
##  4 alpha1         4   -2.31 
##  5 alpha1         5   -1.85 
##  6 alpha1         6    3.10 
##  7 alpha1         7    4.66 
##  8 alpha1         8    1.61 
##  9 alpha1         9    1.21 
## 10 alpha1        10   -0.286
## # ... with 20 more rows
## =============== Parameter = beta =============== 
## # A tibble: 30 x 3
##    variable bill_id estimate
##    <chr>      <dbl>    <dbl>
##  1 beta1          1    -5.75
##  2 beta1          2    -4.15
##  3 beta1          3    -4.94
##  4 beta1          4    -5.98
##  5 beta1          5    -5.58
##  6 beta1          6    -1.10
##  7 beta1          7    -2.12
##  8 beta1          8    -4.90
##  9 beta1          9    -5.18
## 10 beta1         10    -5.38
## # ... with 20 more rows
```

### Dynamic multinomial IRT

Using a simulated voting data (long-format dataframe) `m_data_dyn` and a
simulated matching bill data `sim_match` for dynamic estimation (See
more detail of across time estimation in Bailey (2007)).

``` r
# multinomial dynamic
data("m_data_dyn")
data("sim_match")

m_mlt_d <- make_rollcall(m_data_dyn,
                         unit_id = "unit",
                         bill_id = "bill",
                         vote_col = "vote",
                         drop_unanimous = TRUE)
## Remove some bills because they are unanimous votings: 4 9 14 20 32 53 58 77 79 96 100

# Generating matching bill indicator for across time estimation
bill_match <- make_bill_match(m_mlt_d, sim_match)
init_mlt_d <- make_init(m_mlt_d, 
                        model = "multi_dyn",
                        T = length(unique(m_data_dyn$time)),
                        bill_match = bill_match,
                        constraint = 1)

dyn_ops_mlt <- make_dyn_options(m_data_dyn,
                                unit_id = "unit",
                                bill_id = "bill",
                                time_id = "time",
                                vote_col = "vote",
                                bill_match = bill_match,
                                drop_unanimous = TRUE)
## Remove some bills because they are unanimous votings: 4 9 14 20 32 53 58 77 79 96 100

fit_mlt_d <- pgIRT(m_mlt_d,
                   mode = "multi_dyn",
                   init = init_mlt_d,
                   constraint = 1,
                   dyn_options = dyn_ops_mlt,
                   verbose = 20)
## =========================================================================
## Polya-Gamma data augmentation Item Response Theory Model via EM Algorithm
## Model = Dynamic Multinomial (Stick-Breaking) 
## =========================================================================
## Iteration 20 : eval = alpha2 0.0001816278 elapsed 0.1 sec
## Iteration 40 : eval = alpha2 2.649164e-06 elapsed 0.2 sec
## Model converged at iteration 46 : 0.2 sec
```

Returning 99% confidence interval:

``` r
boot_mlt_d <- pgIRT_boot(fit_mlt_d, boot = 100, verbose = 20)
## ================================================================
## Parametric Bootstrap for pgIRT ( Dynamic Multinomial )
## ================================================================
## Boostrap 20 DONE : 4.6 sec
## Boostrap 40 DONE : 8.7 sec
## Boostrap 60 DONE : 13 sec
## Boostrap 80 DONE : 17.1 sec
## Boostrap 100 DONE : 21.1 sec

summary(boot_mlt_d, parameter = "theta", ci = .99)
## ==================== Parameter = theta ==================== 
## # A tibble: 1,000 x 7
##    variable unit_id session ci    estimate   lwr   upr
##    <chr>      <dbl>   <int> <chr>    <dbl> <dbl> <dbl>
##  1 theta          1       1 99%       2.34  2.09  2.58
##  2 theta          1       2 99%       2.33  2.07  2.58
##  3 theta          1       3 99%       2.34  2.08  2.58
##  4 theta          1       4 99%       2.33  2.10  2.57
##  5 theta          1       5 99%       2.33  2.13  2.58
##  6 theta          1       6 99%       2.33  2.13  2.57
##  7 theta          1       7 99%       2.32  2.10  2.57
##  8 theta          1       8 99%       2.31  2.12  2.56
##  9 theta          1       9 99%       2.32  2.12  2.56
## 10 theta          1      10 99%       2.33  2.12  2.59
## # ... with 990 more rows
```

## References

-   Bailey, M. A. (2007). “Comparable preference estimates across time
    and institutions for the court, congress, and presidency”. *American
    Journal of Political Science*, 51(3), 433-448.
-   Goplerud, M. (2019). “A Multinomial Framework for Ideal Point
    Estimation”. *Political Analysis*, 27(1), 69-89.
-   Lewis, J. B., & Poole, K. T. (2004). “Measuring bias and uncertainty
    in ideal point estimates via the parametric bootstrap”. *Political
    Analysis*, 12(2), 105-127.
-   Martin A.D., Quinn K.M. & Park J.H. (2011). “MCMCpack: Markov Chain
    Monte Carlo in R.” *Journal of Statistical Software*, 42(9), 22.
