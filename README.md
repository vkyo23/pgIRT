
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pgIRT

<!-- badges: start -->
<!-- badges: end -->

`pgIRT` is an **R** package that implements Item Response Theory (IRT)
model with Polya-Gamma data augmentation and EM Algorithm. The function
is available for the data containing K \>= 2 response category (binary
\~ K-multinomial) and having different categories across items. The
algorithm here is based on Goplerud (2019). This package also utilizes
the parametric bootstrap method proposed by Lewis and Poole (2004) to
estimate confidence interval.

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
#> = Format ------> Binary
#> = Model -------> Default pgIRT 
#> =
#> ---------- Implementing EM ----------
#> Model converged at iteration 28 : 0.9 sec.

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

| time | rcid | name      | value | judge_id |
|-----:|-----:|:----------|------:|---------:|
|    0 |    1 | Rehnquist |     0 |        1 |
|    0 |    1 | Stevens   |     1 |        2 |
|    0 |    1 | O.Connor  |     0 |        3 |
|    0 |    1 | Scalia    |     0 |        4 |
|    0 |    1 | Kennedy   |     1 |        5 |
|    0 |    1 | Souter    |     1 |        6 |

``` r
# Generating roll-call matrix from long-format dataframe
mat_dyn <- make_rollcall(long_df, 
                         unit_id = "judge_id", 
                         bill_id = "rcid",
                         vote_col = "value")
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
#> = Format ------> Binary
#> = Model -------> Dynamic pgIRT 
#> =
#> ---------- Implementing EM ----------
#> Iteration 10: eval = 5.72897e-05
#> Iteration 20: eval = 1.91845e-06
#> Model converged at iteration 22 : 0 sec.
```

To compute confidence interval, run bootstrap via `pgIRT_boot`:

``` r
# Bootstrap
boot <- pgIRT_boot(fit_dyn, boot = 100, verbose = 20)
#> ================================================================
#> Parametric Bootstrap for pgIRT ( Dynamic model )
#> ================================================================
#> Boostrap 20 DONE : 1.2 sec
#> Boostrap 40 DONE : 2.4 sec
#> Boostrap 60 DONE : 3.6 sec
#> Boostrap 80 DONE : 4.7 sec
#> Boostrap 100 DONE : 6 sec

summary(boot, parameter = "theta", ci = .95)
#> ==================== Parameter = theta ==================== 
#> # A tibble: 99 × 7
#>    unit      variable session ci       lwr estimate    upr
#>    <chr>     <chr>      <int> <chr>  <dbl>    <dbl>  <dbl>
#>  1 Rehnquist theta          1 95%    0.937    1.08   1.20 
#>  2 Stevens   theta          1 95%   -1.67    -1.55  -1.45 
#>  3 O.Connor  theta          1 95%    0.285    0.450  0.622
#>  4 Scalia    theta          1 95%    1.12     1.24   1.35 
#>  5 Kennedy   theta          1 95%    0.294    0.460  0.605
#>  6 Souter    theta          1 95%   -1.20    -1.08  -0.937
#>  7 Thomas    theta          1 95%    1.25     1.38   1.49 
#>  8 Ginsburg  theta          1 95%   -1.27    -1.11  -0.970
#>  9 Breyer    theta          1 95%   -0.937   -0.813 -0.696
#> 10 Rehnquist theta          2 95%    0.939    1.08   1.20 
#> # … with 89 more rows
```

### Multinomial IRT

Using simulated a multinomial response data (`m_data_dyn`) data.

``` r
data(m_data_dyn)

m_mlt_d <- make_rollcall(m_data_dyn,
                         unit_id = "unit",
                         bill_id = "bill",
                         vote_col = "vote") %>% 
  clean_rollcall()
#> Remove some bills because they are unanimous votings: 4 9 14 20 32 53 58 77 79 96 100

fit_mlt <- pgIRT(m_mlt_d,
                 model = "default",
                 constraint = 1,
                 verbose = 20)
#> =========================================================================
#> Polya-Gamma data augmentation Item Response Theory Model via EM Algorithm
#> =========================================================================
#> = Format ------> Multinomial (# of categories: Min = 2 / Max = 3 )
#> = Model -------> Default pgIRT 
#> =
#> ---------- Implementing EM ----------
#> Iteration 20: eval = 5.92125e-05
#> Iteration 40: eval = 1.22543e-05
#> Iteration 60: eval = 3.6913e-06
#> Iteration 80: eval = 1.31633e-06
#> Model converged at iteration 87 : 0.6 sec.

summary(fit_mlt)$theta
#> # A tibble: 100 × 3
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
#> # … with 90 more rows
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
#> Remove some bills because they are unanimous votings: 4 9 14 20 32 53 58 77 79 96 100

fit_mlt_d <- pgIRT(m_mlt_d,
                   mode = "dynamic",
                   constraint = 1,
                   dyn_options = dyn_ops_mlt,
                   verbose = 20)
#> =========================================================================
#> Polya-Gamma data augmentation Item Response Theory Model via EM Algorithm
#> =========================================================================
#> = Format ------> Multinomial (# of categories: Min = 2 / Max = 3 )
#> = Model -------> Dynamic pgIRT 
#> =
#> ---------- Implementing EM ----------
#> Iteration 20: eval = 1.28497e-05
#> Iteration 40: eval = 2.30141e-06
#> Iteration 60: eval = 1.35263e-06
#> Model converged at iteration 74 : 0.8 sec.
```

Returning 99% confidence interval:

``` r
boot_mlt_d <- pgIRT_boot(fit_mlt_d, boot = 100, verbose = 20)
#> ================================================================
#> Parametric Bootstrap for pgIRT ( Dynamic model )
#> ================================================================
#> Boostrap 20 DONE : 11.5 sec
#> Boostrap 40 DONE : 25.4 sec
#> Boostrap 60 DONE : 37.4 sec
#> Boostrap 80 DONE : 50 sec
#> Boostrap 100 DONE : 65.5 sec

summary(boot_mlt_d, parameter = "theta", ci = .99)
#> ==================== Parameter = theta ==================== 
#> # A tibble: 985 × 7
#>    unit  variable session ci       lwr estimate     upr
#>    <chr> <chr>      <int> <chr>  <dbl>    <dbl>   <dbl>
#>  1 1     theta          1 99%    0.921    1.01   1.09  
#>  2 2     theta          1 99%    0.380    0.467  0.560 
#>  3 3     theta          1 99%   -0.358   -0.317 -0.262 
#>  4 4     theta          1 99%   -0.276   -0.217 -0.153 
#>  5 5     theta          1 99%   -0.241   -0.177 -0.0931
#>  6 6     theta          1 99%   -0.365   -0.318 -0.261 
#>  7 7     theta          1 99%   -0.378   -0.323 -0.272 
#>  8 8     theta          1 99%   -0.353   -0.299 -0.236 
#>  9 9     theta          1 99%   -0.361   -0.317 -0.269 
#> 10 10    theta          1 99%   -0.237   -0.179 -0.119 
#> # … with 975 more rows
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
