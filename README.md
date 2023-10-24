
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

## Usage

### Binary IRT

Monte Carlo simulation:

``` r
library(pgIRT)
require(dplyr)
require(tidyr)
set.seed(1)

# Number of individuals and items
I <- 100 
J <- 1000

# DGP
theta_true <- seq(-2, 2, length = I)
alpha_true <- rnorm(J)
beta_true <- rnorm(J)
Y_star <- cbind(1, theta_true) %*% rbind(alpha_true, beta_true)
Y <- matrix(rbinom(I * J, 1, plogis(Y_star)), I, J)
```

Function `pgIRT` implements IRT. Please note that the elements of an
input matrix should start from 1, not 0. Then, I convert 0 into 2.

``` r
Y[Y == 0] <- 2 
pgfit <- pgIRT(Y,
               model = "default",
               constraint = which.max(theta_true))
#> =========================================================================
#> Polya-Gamma data augmentation Item Response Theory Model via EM Algorithm
#> =========================================================================
#> * Format ------> Binary
#> * Model -------> Default pgIRT 
#> *
#> === Expectation-Maximization ===
#> Model converged at iteration 28 : 0.4 sec.

cor(pgfit$parameter$theta, theta_true)
#>              [,1]
#> theta_1 0.9957472
```

### Dynamic binary IRT

Using `Rehnquist` data from `MCMCpack`. To introduce some auxiliary
function of `pgIRT`, I first construct long-format dataframe.

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
#> * Created 9 x 485 matrix.
rownames(mat_dyn) <- jname$name
mat_dyn[mat_dyn == 0] <- 2
constraint_dyn <- jname$judge_id[jname$name == "Thomas"]

# Options for dynamic estimation
dyn_ops <- make_dyn_options(long_df,
                            unit_id = "judge_id",
                            bill_id = "rcid",
                            time_id = "time",
                            vote_col = "value")
fit_dyn <- pgIRT(mat_dyn,
                 model = "dynamic",
                 constraint = constraint_dyn,
                 dyn_options = dyn_ops,
                 verbose = 10)
#> =========================================================================
#> Polya-Gamma data augmentation Item Response Theory Model via EM Algorithm
#> =========================================================================
#> * Format ------> Binary
#> * Model -------> Dynamic pgIRT 
#> *
#> === Expectation-Maximization ===
#> Iteration 10: eval = 6.69532e-05
#> Iteration 20: eval = 6.36456e-06
#> Iteration 30: eval = 1.31581e-06
#> Model converged at iteration 33 : 0 sec.
```

To compute confidence interval, run bootstrap via `pgIRT_boot`:

``` r
# Bootstrap
boot <- pgIRT_boot(fit_dyn, boot = 100, verbose = 20)
#> ================================================================
#> Parametric Bootstrap for pgIRT ( Dynamic model )
#> ================================================================
#> Boostrap 20 DONE : 1.1 sec
#> Boostrap 40 DONE : 2.1 sec
#> Boostrap 60 DONE : 2.9 sec
#> Boostrap 80 DONE : 3.6 sec
#> Boostrap 100 DONE : 4.7 sec

summary(boot, parameter = "theta", ci = .95)
#> ==================== Parameter = theta ==================== 
#> # A tibble: 99 x 7
#>    unit      variable session ci       lwr estimate    upr
#>    <chr>     <chr>      <int> <chr>  <dbl>    <dbl>  <dbl>
#>  1 Rehnquist theta          1 95%    0.838    0.988  1.13 
#>  2 Stevens   theta          1 95%   -2.45    -2.33  -2.23 
#>  3 O.Connor  theta          1 95%    0.156    0.303  0.448
#>  4 Scalia    theta          1 95%    1.49     1.71   1.90 
#>  5 Kennedy   theta          1 95%    0.162    0.314  0.480
#>  6 Souter    theta          1 95%   -0.932   -0.713 -0.543
#>  7 Thomas    theta          1 95%    1.72     1.89   2.06 
#>  8 Ginsburg  theta          1 95%   -1.08    -0.931 -0.725
#>  9 Breyer    theta          1 95%   -0.971   -0.740 -0.547
#> 10 Rehnquist theta          2 95%    0.841    0.975  1.10 
#> # ... with 89 more rows
```

### Multinomial IRT

Using a simulated multinomial response data (`m_data_dyn`) data.

``` r
data(m_data_dyn)

m_mlt_d <- make_rollcall(m_data_dyn,
                         unit_id = "unit",
                         bill_id = "bill",
                         vote_col = "vote") %>% 
  clean_rollcall()
#> * Created 100 x 120 matrix.
#> * Removed unanimous items: 4 9 14 20 32 53 58 77 79 96 100 
#> * Remaining: 100 x 109

fit_mlt <- pgIRT(m_mlt_d,
                 model = "default",
                 constraint = 1,
                 verbose = 20)
#> =========================================================================
#> Polya-Gamma data augmentation Item Response Theory Model via EM Algorithm
#> =========================================================================
#> * Format ------> Multinomial ( # of categories: Min = 2 / Max = 3 )
#> * Model -------> Default pgIRT 
#> *
#> === Expectation-Maximization ===
#> Iteration 20: eval = 5.92125e-05
#> Iteration 40: eval = 1.22543e-05
#> Iteration 60: eval = 3.6913e-06
#> Iteration 80: eval = 1.31633e-06
#> Model converged at iteration 87 : 0.2 sec.

summary(fit_mlt)$theta
#> # A tibble: 100 x 3
#>    unit  variable estimate
#>    <chr> <chr>       <dbl>
#>  1 1     theta       5.46 
#>  2 2     theta       1.17 
#>  3 3     theta      -1.45 
#>  4 4     theta      -2.99 
#>  5 5     theta      -1.01 
#>  6 6     theta      -1.45 
#>  7 7     theta      -1.46 
#>  8 8     theta      -1.44 
#>  9 9     theta      -1.45 
#> 10 10    theta      -0.747
#> # ... with 90 more rows
```

### Dynamic multinomial IRT

Using the simulated response data and a simulated matching bill data
`sim_match` for dynamic estimation (See more detail of across time
estimation in Bailey (2007)).

``` r
# multinomial dynamic
data(sim_match)

# Generating matching bill indicator for across time estimation
bill_match <- make_bill_match(m_mlt_d, sim_match)

dyn_ops_mlt <- make_dyn_options(m_data_dyn,
                                unit_id = "unit",
                                bill_id = "bill",
                                time_id = "time",
                                vote_col = "vote",
                                add_matched_bill = bill_match,
                                clean = TRUE)
#> * Removed unanimous items: 4 9 14 20 32 53 58 77 79 96 100

fit_mlt_d <- pgIRT(m_mlt_d,
                   mode = "dynamic",
                   constraint = 1,
                   dyn_options = dyn_ops_mlt,
                   verbose = 20)
#> =========================================================================
#> Polya-Gamma data augmentation Item Response Theory Model via EM Algorithm
#> =========================================================================
#> * Format ------> Multinomial ( # of categories: Min = 2 / Max = 3 )
#> * Model -------> Dynamic pgIRT 
#> *
#> === Expectation-Maximization ===
#> Iteration 20: eval = 8.52446e-05
#> Iteration 40: eval = 1.75347e-05
#> Iteration 60: eval = 5.34718e-06
#> Iteration 80: eval = 1.76179e-06
#> Model converged at iteration 92 : 0.3 sec.
```

Returning 99% confidence interval:

``` r
boot_mlt_d <- pgIRT_boot(fit_mlt_d, boot = 100, verbose = 20)
#> ================================================================
#> Parametric Bootstrap for pgIRT ( Dynamic model )
#> ================================================================
#> Boostrap 20 DONE : 8.4 sec
#> Boostrap 40 DONE : 15.5 sec
#> Boostrap 60 DONE : 22.9 sec
#> Boostrap 80 DONE : 30.9 sec
#> Boostrap 100 DONE : 38.9 sec

summary(boot_mlt_d, parameter = "theta", ci = .99)
#> ==================== Parameter = theta ==================== 
#> # A tibble: 985 x 7
#>    unit  variable session ci       lwr estimate     upr
#>    <chr> <chr>      <int> <chr>  <dbl>    <dbl>   <dbl>
#>  1 1     theta          1 99%    1.76     1.99   2.22  
#>  2 2     theta          1 99%    0.454    0.620  0.756 
#>  3 3     theta          1 99%   -0.469   -0.313 -0.198 
#>  4 4     theta          1 99%   -1.09    -0.894 -0.753 
#>  5 5     theta          1 99%   -0.377   -0.245 -0.0884
#>  6 6     theta          1 99%   -0.436   -0.314 -0.191 
#>  7 7     theta          1 99%   -0.517   -0.314 -0.189 
#>  8 8     theta          1 99%   -0.428   -0.305 -0.232 
#>  9 9     theta          1 99%   -0.442   -0.312 -0.186 
#> 10 10    theta          1 99%    0.212    0.362  0.446 
#> # ... with 975 more rows
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
