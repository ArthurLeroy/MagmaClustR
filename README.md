
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MagmaClustR <img src="man/figures/logo.png" align="right" width="120" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/ArthurLeroy/MagmaClustR/workflows/R-CMD-check/badge.svg)](https://github.com/ArthurLeroy/MagmaClustR/actions)
[![codecov](https://codecov.io/gh/ArthurLeroy/MagmaClustR/branch/master/graph/badge.svg?token=KH7SOKOLKX)](https://codecov.io/gh/ArthurLeroy/MagmaClustR)
<!-- badges: end -->

The *MagmaClustR* package implements two main algorithms, called Magma
(*Leroy et al., 2022*) and MagmaClust (*Leroy et al., 2020*), using a
multi-task Gaussian processes (GP) model to perform predictions for
supervised learning problems. Applications involving functional data,
such as multiple time series, are particularly well-handled. Theses
approaches leverage the learning of cluster-specific mean processes,
which are common across similar tasks, to provide enhanced prediction
performances (even far from data points) at a linear computational cost
(in the number of tasks). MagmaClust is a generalisation of Magma where
the tasks are simultaneously clustered into groups, each being
associated to a specific mean process. User-oriented functions in the
package are decomposed into training, prediction and plotting functions.
Some basic features of standard GPs are also implemented.

Leroy, A., Latouche, P., Guedj, B., Gey, S. MAGMA: inference and
prediction using multi-task Gaussian processes with common mean. *Mach
Learn* **111**, 1821–1849 (2022).
<https://doi.org/10.1007/s10994-022-06172-1>

Leroy, A., Latouche, P., Guedj, B., & Gey, S. Cluster-Specific
Predictions with Multi-Task Gaussian Processes. *arXiv preprint* (2020).
<https://arxiv.org/abs/2011.07866>

## Installation

You can install the released version of MagmaClustR from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("MagmaClustR")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ArthurLeroy/MagmaClustR")
```

## Example: Magma

Here is a basic example on how to simulate a dataset with the adequate
format, then train a Magma model and use it to perform predictions.

### Data generation

``` r
library(MagmaClustR)
## Simulate a dataset with 11 individuals, each observed at 10 input locations
set.seed(2)
data_magma <- simu_db(M = 11, N = 10, common_input = FALSE)
## Split individuals into training and prediction sets, and define test points
magma_train <- data_magma %>% subset(ID %in% 1:10)
magma_pred <- data_magma %>% subset(ID == 11) %>% head(5)
magma_test <- data_magma %>% subset(ID == 11) %>% tail(5)

data_magma
#> # A tibble: 110 x 3
#>    ID    Output Input
#>    <chr>  <dbl> <dbl>
#>  1 1       1.52  1.9 
#>  2 1      -5.76  3.3 
#>  3 1      -3.78  3.35
#>  4 1       7.23  5.25
#>  5 1      14.9   6.2 
#>  6 1       7.48  7.4 
#>  7 1       7.21  8.15
#>  8 1      11.5   8.2 
#>  9 1      13.1   8.9 
#> 10 1       8.72  9.5 
#> # ... with 100 more rows
```

### Training and prediction with Magma

``` r
model <- train_magma(data = magma_train)
#> The 'prior_mean' argument has not been specified. The hyper_prior mean function is thus set to be 0 everywhere.
#>  
#> The 'ini_hp_0' argument has not been specified. Random values of hyper-parameters for the mean process are used as initialisation.
#>  
#> Called from: train_magma(data = magma_train)
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#238: if ("ID" %in% names(hp_0)) {
#>     hp_0 = hp_0[names(hp_0) != "ID"]
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#243: if (kern_i %>% is.function()) {
#>     if (ini_hp_i %>% is.null()) {
#>         stop("When using a custom kernel function the 'ini_hp_i' argument is ", 
#>             "mandatory, in order to provide the name of the hyper-parameters. ", 
#>             "You can use the function 'hp()' to easily generate a tibble of random", 
#>             " hyper-parameters with the desired format for initialisation.")
#>     }
#> } else {
#>     if (ini_hp_i %>% is.null()) {
#>         hp_i <- hp(kern_i, list_ID = list_ID, common_hp = common_hp, 
#>             noise = TRUE)
#>         cat("The 'ini_hp_i' argument has not been specified. Random values of", 
#>             "hyper-parameters for the individal processes are used as", 
#>             "initialisation.\n \n")
#>     }
#>     else if (!("ID" %in% names(ini_hp_i))) {
#>         hp_i <- tibble::tibble(ID = list_ID, dplyr::bind_rows(ini_hp_i))
#>     }
#>     else if (!(all(as.character(ini_hp_i$ID) %in% as.character(list_ID)) & 
#>         all(as.character(list_ID) %in% as.character(ini_hp_i$ID)))) {
#>         stop("The 'ID' column in 'ini_hp_i' is different from the 'ID' of the ", 
#>             "'data'.")
#>     }
#>     else {
#>         hp_i <- ini_hp_i
#>     }
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#253: if (ini_hp_i %>% is.null()) {
#>     hp_i <- hp(kern_i, list_ID = list_ID, common_hp = common_hp, 
#>         noise = TRUE)
#>     cat("The 'ini_hp_i' argument has not been specified. Random values of", 
#>         "hyper-parameters for the individal processes are used as", 
#>         "initialisation.\n \n")
#> } else if (!("ID" %in% names(ini_hp_i))) {
#>     hp_i <- tibble::tibble(ID = list_ID, dplyr::bind_rows(ini_hp_i))
#> } else if (!(all(as.character(ini_hp_i$ID) %in% as.character(list_ID)) & 
#>     all(as.character(list_ID) %in% as.character(ini_hp_i$ID)))) {
#>     stop("The 'ID' column in 'ini_hp_i' is different from the 'ID' of the ", 
#>         "'data'.")
#> } else {
#>     hp_i <- ini_hp_i
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#254: hp_i <- hp(kern_i, list_ID = list_ID, common_hp = common_hp, 
#>     noise = TRUE)
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#255: cat("The 'ini_hp_i' argument has not been specified. Random values of", 
#>     "hyper-parameters for the individal processes are used as", 
#>     "initialisation.\n \n")
#> The 'ini_hp_i' argument has not been specified. Random values of hyper-parameters for the individal processes are used as initialisation.
#>  
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#278: if (!("noise" %in% names(hp_i))) {
#>     if (common_hp) {
#>         hp_i <- hp_i %>% dplyr::mutate(hp(NULL, noise = T))
#>     }
#>     else {
#>         hp_i <- hp_i %>% dplyr::left_join(hp(NULL, list_ID = hp_i$ID, 
#>             noise = T), by = "ID")
#>     }
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#288: hp_i_ini <- hp_i
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#289: hp_0_ini <- hp_0
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#291: cv <- FALSE
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#292: logL_monitoring <- -Inf
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#293: seq_loglikelihood <- c()
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#296: for (i in 1:n_iter_max) {
#>     t_i_1 <- Sys.time()
#>     post <- e_step(db = data, m_0 = m_0, kern_0 = kern_0, kern_i = kern_i, 
#>         hp_0 = hp_0, hp_i = hp_i, pen_diag = pen_diag)
#>     if (fast_approx) {
#>         seq_loglikelihood <- logL_monitoring(hp_0 = hp_0, hp_i = hp_i, 
#>             db = data, m_0 = m_0, kern_0 = kern_0, kern_i = kern_i, 
#>             post_mean = post$mean, post_cov = post$cov, pen_diag = pen_diag)
#>         cv <- FALSE
#>         break
#>     }
#>     new_hp <- m_step(db = data, m_0 = m_0, kern_0 = kern_0, kern_i = kern_i, 
#>         old_hp_0 = hp_0, old_hp_i = hp_i, post_mean = post$mean, 
#>         post_cov = post$cov, common_hp = common_hp, pen_diag = pen_diag)
#>     new_hp_0 <- new_hp$hp_0
#>     new_hp_i <- new_hp$hp_i
#>     if (any(is.na(new_hp_0)) | any(is.na(new_hp_i))) {
#>         warning(paste0("The M-step encountered an error at iteration : ", 
#>             i))
#>         warning("Training has stopped and hyper-parameters values from the ", 
#>             "last valid iteration are returned.")
#>         break
#>     }
#>     new_logL_monitoring <- logL_monitoring(hp_0 = new_hp_0, hp_i = new_hp_i, 
#>         db = data, m_0 = m_0, kern_0 = kern_0, kern_i = kern_i, 
#>         post_mean = post$mean, post_cov = post$cov, pen_diag = pen_diag)
#>     diff_logL <- new_logL_monitoring - logL_monitoring
#>     if (diff_logL %>% is.nan()) {
#>         diff_logL <- -Inf
#>     }
#>     if (diff_logL < 0) {
#>         warning("The likelihood descreased. Possible numerical issues.")
#>     }
#>     hp_0 <- new_hp_0
#>     hp_i <- new_hp_i
#>     logL_monitoring <- new_logL_monitoring
#>     seq_loglikelihood <- c(seq_loglikelihood, logL_monitoring)
#>     eps <- diff_logL/abs(logL_monitoring)
#>     if (eps %>% is.nan()) {
#>         eps <- 1
#>     }
#>     t_i_2 <- Sys.time()
#>     paste0("EM algorithm, step ", i, ": ", difftime(t_i_2, t_i_1, 
#>         units = "secs") %>% round(2), " seconds \n \n") %>% cat()
#>     paste0("Value of the likelihood: ", logL_monitoring %>% round(5), 
#>         " --- Convergence ratio = ", eps %>% round(5), "\n \n") %>% 
#>         cat()
#>     if (abs(eps) < cv_threshold) {
#>         cat("The EM algorithm successfully converged, training is completed.", 
#>             "\n \n")
#>         cv <- TRUE
#>         break
#>     }
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#299: t_i_1 <- Sys.time()
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#302: post <- e_step(db = data, m_0 = m_0, kern_0 = kern_0, kern_i = kern_i, 
#>     hp_0 = hp_0, hp_i = hp_i, pen_diag = pen_diag)
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#313: if (fast_approx) {
#>     seq_loglikelihood <- logL_monitoring(hp_0 = hp_0, hp_i = hp_i, 
#>         db = data, m_0 = m_0, kern_0 = kern_0, kern_i = kern_i, 
#>         post_mean = post$mean, post_cov = post$cov, pen_diag = pen_diag)
#>     cv <- FALSE
#>     break
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#332: new_hp <- m_step(db = data, m_0 = m_0, kern_0 = kern_0, kern_i = kern_i, 
#>     old_hp_0 = hp_0, old_hp_i = hp_i, post_mean = post$mean, 
#>     post_cov = post$cov, common_hp = common_hp, pen_diag = pen_diag)
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#344: new_hp_0 <- new_hp$hp_0
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#345: new_hp_i <- new_hp$hp_i
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#348: if (any(is.na(new_hp_0)) | any(is.na(new_hp_i))) {
#>     warning(paste0("The M-step encountered an error at iteration : ", 
#>         i))
#>     warning("Training has stopped and hyper-parameters values from the ", 
#>         "last valid iteration are returned.")
#>     break
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#358: new_logL_monitoring <- logL_monitoring(hp_0 = new_hp_0, hp_i = new_hp_i, 
#>     db = data, m_0 = m_0, kern_0 = kern_0, kern_i = kern_i, post_mean = post$mean, 
#>     post_cov = post$cov, pen_diag = pen_diag)
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#370: diff_logL <- new_logL_monitoring - logL_monitoring
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#371: if (diff_logL %>% is.nan()) {
#>     diff_logL <- -Inf
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#375: if (diff_logL < 0) {
#>     warning("The likelihood descreased. Possible numerical issues.")
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#380: hp_0 <- new_hp_0
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#381: hp_i <- new_hp_i
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#382: logL_monitoring <- new_logL_monitoring
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#385: seq_loglikelihood <- c(seq_loglikelihood, logL_monitoring)
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#388: eps <- diff_logL/abs(logL_monitoring)
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#389: if (eps %>% is.nan()) {
#>     eps <- 1
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#394: t_i_2 <- Sys.time()
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#395: paste0("EM algorithm, step ", i, ": ", difftime(t_i_2, t_i_1, 
#>     units = "secs") %>% round(2), " seconds \n \n") %>% cat()
#> EM algorithm, step 1: 18.48 seconds 
#>  
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#402: paste0("Value of the likelihood: ", logL_monitoring %>% round(5), 
#>     " --- Convergence ratio = ", eps %>% round(5), "\n \n") %>% 
#>     cat()
#> Value of the likelihood: -379.55987 --- Convergence ratio = Inf
#>  
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#412: if (abs(eps) < cv_threshold) {
#>     cat("The EM algorithm successfully converged, training is completed.", 
#>         "\n \n")
#>     cv <- TRUE
#>     break
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#299: t_i_1 <- Sys.time()
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#302: post <- e_step(db = data, m_0 = m_0, kern_0 = kern_0, kern_i = kern_i, 
#>     hp_0 = hp_0, hp_i = hp_i, pen_diag = pen_diag)
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#313: if (fast_approx) {
#>     seq_loglikelihood <- logL_monitoring(hp_0 = hp_0, hp_i = hp_i, 
#>         db = data, m_0 = m_0, kern_0 = kern_0, kern_i = kern_i, 
#>         post_mean = post$mean, post_cov = post$cov, pen_diag = pen_diag)
#>     cv <- FALSE
#>     break
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#332: new_hp <- m_step(db = data, m_0 = m_0, kern_0 = kern_0, kern_i = kern_i, 
#>     old_hp_0 = hp_0, old_hp_i = hp_i, post_mean = post$mean, 
#>     post_cov = post$cov, common_hp = common_hp, pen_diag = pen_diag)
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#344: new_hp_0 <- new_hp$hp_0
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#345: new_hp_i <- new_hp$hp_i
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#348: if (any(is.na(new_hp_0)) | any(is.na(new_hp_i))) {
#>     warning(paste0("The M-step encountered an error at iteration : ", 
#>         i))
#>     warning("Training has stopped and hyper-parameters values from the ", 
#>         "last valid iteration are returned.")
#>     break
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#358: new_logL_monitoring <- logL_monitoring(hp_0 = new_hp_0, hp_i = new_hp_i, 
#>     db = data, m_0 = m_0, kern_0 = kern_0, kern_i = kern_i, post_mean = post$mean, 
#>     post_cov = post$cov, pen_diag = pen_diag)
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#370: diff_logL <- new_logL_monitoring - logL_monitoring
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#371: if (diff_logL %>% is.nan()) {
#>     diff_logL <- -Inf
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#375: if (diff_logL < 0) {
#>     warning("The likelihood descreased. Possible numerical issues.")
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#380: hp_0 <- new_hp_0
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#381: hp_i <- new_hp_i
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#382: logL_monitoring <- new_logL_monitoring
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#385: seq_loglikelihood <- c(seq_loglikelihood, logL_monitoring)
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#388: eps <- diff_logL/abs(logL_monitoring)
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#389: if (eps %>% is.nan()) {
#>     eps <- 1
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#394: t_i_2 <- Sys.time()
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#395: paste0("EM algorithm, step ", i, ": ", difftime(t_i_2, t_i_1, 
#>     units = "secs") %>% round(2), " seconds \n \n") %>% cat()
#> EM algorithm, step 2: 9.83 seconds 
#>  
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#402: paste0("Value of the likelihood: ", logL_monitoring %>% round(5), 
#>     " --- Convergence ratio = ", eps %>% round(5), "\n \n") %>% 
#>     cat()
#> Value of the likelihood: -374.3461 --- Convergence ratio = 0.01393
#>  
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#412: if (abs(eps) < cv_threshold) {
#>     cat("The EM algorithm successfully converged, training is completed.", 
#>         "\n \n")
#>     cv <- TRUE
#>     break
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#299: t_i_1 <- Sys.time()
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#302: post <- e_step(db = data, m_0 = m_0, kern_0 = kern_0, kern_i = kern_i, 
#>     hp_0 = hp_0, hp_i = hp_i, pen_diag = pen_diag)
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#313: if (fast_approx) {
#>     seq_loglikelihood <- logL_monitoring(hp_0 = hp_0, hp_i = hp_i, 
#>         db = data, m_0 = m_0, kern_0 = kern_0, kern_i = kern_i, 
#>         post_mean = post$mean, post_cov = post$cov, pen_diag = pen_diag)
#>     cv <- FALSE
#>     break
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#332: new_hp <- m_step(db = data, m_0 = m_0, kern_0 = kern_0, kern_i = kern_i, 
#>     old_hp_0 = hp_0, old_hp_i = hp_i, post_mean = post$mean, 
#>     post_cov = post$cov, common_hp = common_hp, pen_diag = pen_diag)
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#344: new_hp_0 <- new_hp$hp_0
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#345: new_hp_i <- new_hp$hp_i
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#348: if (any(is.na(new_hp_0)) | any(is.na(new_hp_i))) {
#>     warning(paste0("The M-step encountered an error at iteration : ", 
#>         i))
#>     warning("Training has stopped and hyper-parameters values from the ", 
#>         "last valid iteration are returned.")
#>     break
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#358: new_logL_monitoring <- logL_monitoring(hp_0 = new_hp_0, hp_i = new_hp_i, 
#>     db = data, m_0 = m_0, kern_0 = kern_0, kern_i = kern_i, post_mean = post$mean, 
#>     post_cov = post$cov, pen_diag = pen_diag)
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#370: diff_logL <- new_logL_monitoring - logL_monitoring
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#371: if (diff_logL %>% is.nan()) {
#>     diff_logL <- -Inf
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#375: if (diff_logL < 0) {
#>     warning("The likelihood descreased. Possible numerical issues.")
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#380: hp_0 <- new_hp_0
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#381: hp_i <- new_hp_i
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#382: logL_monitoring <- new_logL_monitoring
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#385: seq_loglikelihood <- c(seq_loglikelihood, logL_monitoring)
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#388: eps <- diff_logL/abs(logL_monitoring)
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#389: if (eps %>% is.nan()) {
#>     eps <- 1
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#394: t_i_2 <- Sys.time()
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#395: paste0("EM algorithm, step ", i, ": ", difftime(t_i_2, t_i_1, 
#>     units = "secs") %>% round(2), " seconds \n \n") %>% cat()
#> EM algorithm, step 3: 9.7 seconds 
#>  
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#402: paste0("Value of the likelihood: ", logL_monitoring %>% round(5), 
#>     " --- Convergence ratio = ", eps %>% round(5), "\n \n") %>% 
#>     cat()
#> Value of the likelihood: -374.17964 --- Convergence ratio = 0.00044
#>  
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#412: if (abs(eps) < cv_threshold) {
#>     cat("The EM algorithm successfully converged, training is completed.", 
#>         "\n \n")
#>     cv <- TRUE
#>     break
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#413: cat("The EM algorithm successfully converged, training is completed.", 
#>     "\n \n")
#> The EM algorithm successfully converged, training is completed. 
#>  
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#417: cv <- TRUE
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#418: break
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#422: if (!cv & (i == n_iter_max)) {
#>     warning("The EM algorithm has reached the maximum number of iterations ", 
#>         "before convergence, training might be sub-optimal \n \n")
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#430: if (!is.null(grid_inputs)) {
#>     cat("Start evaluating hyper-posterior distribution of the mean process", 
#>         "on the provided grid of inputs... \n \n")
#>     post <- hyperposterior(data = data, hp_0 = hp_0, hp_i = hp_i, 
#>         kern_0 = kern_0, kern_i = kern_i, prior_mean = prior_mean, 
#>         grid_inputs = grid_inputs, pen_diag = pen_diag)
#>     cat("Done!\n \n")
#> } else {
#>     post$pred <- tibble::tibble(Input = post$mean %>% dplyr::pull(.data$Input), 
#>         Mean = post$mean %>% dplyr::pull(.data$Output), Var = post$cov %>% 
#>             diag() %>% as.vector())
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#449: post$pred <- tibble::tibble(Input = post$mean %>% dplyr::pull(.data$Input), 
#>     Mean = post$mean %>% dplyr::pull(.data$Output), Var = post$cov %>% 
#>         diag() %>% as.vector())
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#457: fct_args <- list(data = data, prior_mean = prior_mean, ini_hp_0 = hp_0_ini, 
#>     ini_hp_i = hp_i_ini, kern_0 = kern_0, kern_i = kern_i, common_hp = common_hp, 
#>     grid_inputs = grid_inputs, pen_diag = pen_diag, n_iter_max = n_iter_max, 
#>     cv_threshold = cv_threshold)
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#471: t_2 <- Sys.time()
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#473: list(hp_0 = hp_0, hp_i = hp_i, hyperpost = post, ini_args = fct_args, 
#>     seq_loglikelihood = seq_loglikelihood, converged = cv, training_time = difftime(t_2, 
#>         t_1, units = "secs")) %>% return()

pred  <- pred_magma(data = magma_pred,
                    trained_model = model, 
                    grid_inputs = seq(0,10, 0.01))
#> The hyper-posterior distribution of the mean process provided in 'hyperpost' argument isn't evaluated on the expected inputs.
#>  
#>  Start evaluating the hyper-posterior on the correct inputs...
#>  
#> The 'prior_mean' argument has not been specified. The hyper-prior mean function is thus set to be 0 everywhere.
#>  
#> Done!
#> 
```

<img src="man/figures/README-train_and_predict_Magma-1.png" width="100%" />

Note that the `grid_inputs` argument is optional. It merely allows users
to control the grid of values on which the prediction is performed.

### Display the resulting predictions

Several arguments are available in a specific plotting function to offer
additional control in the display of results. For instance, the GP
prediction can be represented as a heatmap of probabilities:

``` r
plot_gp(pred_gp = pred,
        data = magma_pred,
        data_train = magma_train,
        prior_mean = model$hyperpost$mean,
        heatmap = TRUE) 
```

<img src="man/figures/README-display_Magma-1.png" width="100%" />

Additionally, it is also possible to create animated representations by
using functions that generate GIFs. For instance, below, the true
testing points have been represented as red dots and we can observe how
the prediction evolves as we add more data points to our prediction
dataset.

``` r
pred_gif  <- pred_gif(data = magma_pred,
                      trained_model = model,
                      grid_inputs = seq(0, 10, 0.01))
#>  => 1 => 2 => 3 => 4 => 5

plot_gif(pred_gp = pred_gif,
        data = magma_pred,
        data_train = magma_train,
        prior_mean = model$hyperpost$mean) + 
  ggplot2::geom_point(data = magma_test,
                       ggplot2::aes(x = Input, y = Output),
                       color = 'red')
```

<img src="man/figures/README-gif_Magma-1.gif" width="100%" />

Note that the `grid_inputs` argument is optional. It merely allows users
to control the grid of values on which the prediction is performed.

## Example: MagmaClust

Here is a basic example on how to simulate a dataset with the adequate
format, then train a MagmaClust model and use it to perform simultaneous
clustering and predictions.

### Data generation

``` r
## Simulate a dataset containing 3 clusters of 4 individuals, each observed at 10 input locations
set.seed(2) 
data_magmaclust <- simu_db(M = 4, N = 10, K = 3) 
## Split individuals into training and prediction sets, and define test points
list_ID = unique(data_magmaclust$ID)
magmaclust_train <- data_magmaclust %>% subset(ID %in% list_ID[1:11])
magmaclust_pred <- data_magmaclust %>% subset(ID == list_ID[12]) %>% head(5)
magmaclust_test <- data_magmaclust %>% subset(ID == list_ID[12]) %>% tail(5)

data_magmaclust
#> # A tibble: 120 x 3
#>    ID         Output Input
#>    <chr>       <dbl> <dbl>
#>  1 ID1-Clust1 -11.1   0.25
#>  2 ID1-Clust1  -7.64  0.8 
#>  3 ID1-Clust1  -4.91  2   
#>  4 ID1-Clust1 -13.2   4.2 
#>  5 ID1-Clust1 -14.3   4.6 
#>  6 ID1-Clust1 -13.0   6.2 
#>  7 ID1-Clust1 -14.1   6.75
#>  8 ID1-Clust1 -20.3   7.95
#>  9 ID1-Clust1 -14.5   8.85
#> 10 ID1-Clust1 -12.2   9.85
#> # ... with 110 more rows
```

### Training and prediction with MagmaClust

``` r
model_clust <- train_magmaclust(data = magmaclust_train)
#> The number of cluster argument has not been specified. There will be 3 cluster by default. 
#>  
#> The 'ini_hp_i' argument has not been specified. Random values of hyper-parameters for the individual processes are used as initialisation.
#>  
#> The 'ini_hp_k' argument has not been specified. Random values of hyper-parameters for the mean processes are used as initialisation.
#>  
#> The 'prior_mean' argument has not been specified. The hyper_prior mean function is thus set to be 0 everywhere.
#>  
#> VEM algorithm, step 1: 83.24 seconds 
#>  
#> Value of the elbo: -403.8673 --- Convergence ratio = Inf
#>  
#> VEM algorithm, step 2: 28.76 seconds 
#>  
#> Value of the elbo: -383.34763 --- Convergence ratio = 0.05353
#>  
#> VEM algorithm, step 3: 23.63 seconds 
#>  
#> Value of the elbo: -383.08831 --- Convergence ratio = 0.00068
#>  
#> The EM algorithm successfully converged, training is completed. 
#> 

pred_clust  <- pred_magmaclust(data = magmaclust_pred,
                    trained_model = model_clust,
                    grid_inputs = seq(0, 10, 0.01), 
                    plot = FALSE)
#> The hyper-posterior distribution of the mean process provided in 'hyperpost' argument isn't evaluated on the expected inputs. Start evaluating the hyper-posterior on the correct inputs...
#>  
#> The 'prior_mean_k' argument has not been specified. The hyper-prior  mean functions are thus set to be 0 everywhere.
#>  
#> Done!
#> 
```

### Display the resulting predictions

As before, a specific plotting function is provided. For MagmaClust, we
advise to use the heatmap representation in priority, as a mixture of
GPs may not be unimodal in general (and thus prevents the definition of
Credible Interval).

``` r
## Allocate individuals to their most probable cluster to colour them by clusters afterwards
data_train_with_clust = data_allocate_cluster(model_clust)

plot_magmaclust(pred = pred_clust,
                cluster = "all",
                data = magmaclust_pred,
                data_train = data_train_with_clust,
                col_clust = TRUE,
                prior_mean = model_clust$hyperpost$mean,
                y_grid = seq(0, 60, 0.5),
                heatmap = TRUE) 
```

<img src="man/figures/README-display_MagmaClust-1.png" width="100%" />

## Example: in 2-dimensions

Although unidimensional-input problems are easier to visualise, both
Magma and MagmaClust can also be applied with as many covariates as
desired in the model.

### Data generation

``` r
library(MagmaClustR)
## Dataset with 11 individuals, 10 reference input locations and a covariate
set.seed(2) 
data_dim2 <- simu_db(M = 11, N = 10, covariate = TRUE) 
## Split individuals into training and prediction sets, and define test points
dim2_train <- data_dim2 %>% subset(ID %in% 1:10)
dim2_pred <- data_dim2 %>% subset(ID == 11) %>% head(5)
dim2_test <- data_dim2 %>% subset(ID == 11) %>% tail(5)

data_dim2
#> # A tibble: 110 x 4
#>    ID    Output Input Covariate
#>    <chr>  <dbl> <dbl>     <dbl>
#>  1 1     -11.1   0.25     -2   
#>  2 1      -7.64  0.8       1.94
#>  3 1      -4.91  2         4.64
#>  4 1     -13.2   4.2      -3.7 
#>  5 1     -14.3   4.6      -4.24
#>  6 1     -13.0   6.2       0.68
#>  7 1     -14.1   6.75      0.55
#>  8 1     -20.3   7.95     -4.38
#>  9 1     -14.5   8.85      1.74
#> 10 1     -12.2   9.85      4.14
#> # ... with 100 more rows
```

### Training and prediction with Magma

``` r
model_dim2 <- train_magma(data = dim2_train)
#> The 'prior_mean' argument has not been specified. The hyper_prior mean function is thus set to be 0 everywhere.
#>  
#> The 'ini_hp_0' argument has not been specified. Random values of hyper-parameters for the mean process are used as initialisation.
#>  
#> Called from: train_magma(data = dim2_train)
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#238: if ("ID" %in% names(hp_0)) {
#>     hp_0 = hp_0[names(hp_0) != "ID"]
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#243: if (kern_i %>% is.function()) {
#>     if (ini_hp_i %>% is.null()) {
#>         stop("When using a custom kernel function the 'ini_hp_i' argument is ", 
#>             "mandatory, in order to provide the name of the hyper-parameters. ", 
#>             "You can use the function 'hp()' to easily generate a tibble of random", 
#>             " hyper-parameters with the desired format for initialisation.")
#>     }
#> } else {
#>     if (ini_hp_i %>% is.null()) {
#>         hp_i <- hp(kern_i, list_ID = list_ID, common_hp = common_hp, 
#>             noise = TRUE)
#>         cat("The 'ini_hp_i' argument has not been specified. Random values of", 
#>             "hyper-parameters for the individal processes are used as", 
#>             "initialisation.\n \n")
#>     }
#>     else if (!("ID" %in% names(ini_hp_i))) {
#>         hp_i <- tibble::tibble(ID = list_ID, dplyr::bind_rows(ini_hp_i))
#>     }
#>     else if (!(all(as.character(ini_hp_i$ID) %in% as.character(list_ID)) & 
#>         all(as.character(list_ID) %in% as.character(ini_hp_i$ID)))) {
#>         stop("The 'ID' column in 'ini_hp_i' is different from the 'ID' of the ", 
#>             "'data'.")
#>     }
#>     else {
#>         hp_i <- ini_hp_i
#>     }
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#253: if (ini_hp_i %>% is.null()) {
#>     hp_i <- hp(kern_i, list_ID = list_ID, common_hp = common_hp, 
#>         noise = TRUE)
#>     cat("The 'ini_hp_i' argument has not been specified. Random values of", 
#>         "hyper-parameters for the individal processes are used as", 
#>         "initialisation.\n \n")
#> } else if (!("ID" %in% names(ini_hp_i))) {
#>     hp_i <- tibble::tibble(ID = list_ID, dplyr::bind_rows(ini_hp_i))
#> } else if (!(all(as.character(ini_hp_i$ID) %in% as.character(list_ID)) & 
#>     all(as.character(list_ID) %in% as.character(ini_hp_i$ID)))) {
#>     stop("The 'ID' column in 'ini_hp_i' is different from the 'ID' of the ", 
#>         "'data'.")
#> } else {
#>     hp_i <- ini_hp_i
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#254: hp_i <- hp(kern_i, list_ID = list_ID, common_hp = common_hp, 
#>     noise = TRUE)
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#255: cat("The 'ini_hp_i' argument has not been specified. Random values of", 
#>     "hyper-parameters for the individal processes are used as", 
#>     "initialisation.\n \n")
#> The 'ini_hp_i' argument has not been specified. Random values of hyper-parameters for the individal processes are used as initialisation.
#>  
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#278: if (!("noise" %in% names(hp_i))) {
#>     if (common_hp) {
#>         hp_i <- hp_i %>% dplyr::mutate(hp(NULL, noise = T))
#>     }
#>     else {
#>         hp_i <- hp_i %>% dplyr::left_join(hp(NULL, list_ID = hp_i$ID, 
#>             noise = T), by = "ID")
#>     }
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#288: hp_i_ini <- hp_i
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#289: hp_0_ini <- hp_0
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#291: cv <- FALSE
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#292: logL_monitoring <- -Inf
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#293: seq_loglikelihood <- c()
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#296: for (i in 1:n_iter_max) {
#>     t_i_1 <- Sys.time()
#>     post <- e_step(db = data, m_0 = m_0, kern_0 = kern_0, kern_i = kern_i, 
#>         hp_0 = hp_0, hp_i = hp_i, pen_diag = pen_diag)
#>     if (fast_approx) {
#>         seq_loglikelihood <- logL_monitoring(hp_0 = hp_0, hp_i = hp_i, 
#>             db = data, m_0 = m_0, kern_0 = kern_0, kern_i = kern_i, 
#>             post_mean = post$mean, post_cov = post$cov, pen_diag = pen_diag)
#>         cv <- FALSE
#>         break
#>     }
#>     new_hp <- m_step(db = data, m_0 = m_0, kern_0 = kern_0, kern_i = kern_i, 
#>         old_hp_0 = hp_0, old_hp_i = hp_i, post_mean = post$mean, 
#>         post_cov = post$cov, common_hp = common_hp, pen_diag = pen_diag)
#>     new_hp_0 <- new_hp$hp_0
#>     new_hp_i <- new_hp$hp_i
#>     if (any(is.na(new_hp_0)) | any(is.na(new_hp_i))) {
#>         warning(paste0("The M-step encountered an error at iteration : ", 
#>             i))
#>         warning("Training has stopped and hyper-parameters values from the ", 
#>             "last valid iteration are returned.")
#>         break
#>     }
#>     new_logL_monitoring <- logL_monitoring(hp_0 = new_hp_0, hp_i = new_hp_i, 
#>         db = data, m_0 = m_0, kern_0 = kern_0, kern_i = kern_i, 
#>         post_mean = post$mean, post_cov = post$cov, pen_diag = pen_diag)
#>     diff_logL <- new_logL_monitoring - logL_monitoring
#>     if (diff_logL %>% is.nan()) {
#>         diff_logL <- -Inf
#>     }
#>     if (diff_logL < 0) {
#>         warning("The likelihood descreased. Possible numerical issues.")
#>     }
#>     hp_0 <- new_hp_0
#>     hp_i <- new_hp_i
#>     logL_monitoring <- new_logL_monitoring
#>     seq_loglikelihood <- c(seq_loglikelihood, logL_monitoring)
#>     eps <- diff_logL/abs(logL_monitoring)
#>     if (eps %>% is.nan()) {
#>         eps <- 1
#>     }
#>     t_i_2 <- Sys.time()
#>     paste0("EM algorithm, step ", i, ": ", difftime(t_i_2, t_i_1, 
#>         units = "secs") %>% round(2), " seconds \n \n") %>% cat()
#>     paste0("Value of the likelihood: ", logL_monitoring %>% round(5), 
#>         " --- Convergence ratio = ", eps %>% round(5), "\n \n") %>% 
#>         cat()
#>     if (abs(eps) < cv_threshold) {
#>         cat("The EM algorithm successfully converged, training is completed.", 
#>             "\n \n")
#>         cv <- TRUE
#>         break
#>     }
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#299: t_i_1 <- Sys.time()
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#302: post <- e_step(db = data, m_0 = m_0, kern_0 = kern_0, kern_i = kern_i, 
#>     hp_0 = hp_0, hp_i = hp_i, pen_diag = pen_diag)
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#313: if (fast_approx) {
#>     seq_loglikelihood <- logL_monitoring(hp_0 = hp_0, hp_i = hp_i, 
#>         db = data, m_0 = m_0, kern_0 = kern_0, kern_i = kern_i, 
#>         post_mean = post$mean, post_cov = post$cov, pen_diag = pen_diag)
#>     cv <- FALSE
#>     break
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#332: new_hp <- m_step(db = data, m_0 = m_0, kern_0 = kern_0, kern_i = kern_i, 
#>     old_hp_0 = hp_0, old_hp_i = hp_i, post_mean = post$mean, 
#>     post_cov = post$cov, common_hp = common_hp, pen_diag = pen_diag)
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#344: new_hp_0 <- new_hp$hp_0
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#345: new_hp_i <- new_hp$hp_i
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#348: if (any(is.na(new_hp_0)) | any(is.na(new_hp_i))) {
#>     warning(paste0("The M-step encountered an error at iteration : ", 
#>         i))
#>     warning("Training has stopped and hyper-parameters values from the ", 
#>         "last valid iteration are returned.")
#>     break
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#358: new_logL_monitoring <- logL_monitoring(hp_0 = new_hp_0, hp_i = new_hp_i, 
#>     db = data, m_0 = m_0, kern_0 = kern_0, kern_i = kern_i, post_mean = post$mean, 
#>     post_cov = post$cov, pen_diag = pen_diag)
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#370: diff_logL <- new_logL_monitoring - logL_monitoring
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#371: if (diff_logL %>% is.nan()) {
#>     diff_logL <- -Inf
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#375: if (diff_logL < 0) {
#>     warning("The likelihood descreased. Possible numerical issues.")
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#380: hp_0 <- new_hp_0
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#381: hp_i <- new_hp_i
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#382: logL_monitoring <- new_logL_monitoring
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#385: seq_loglikelihood <- c(seq_loglikelihood, logL_monitoring)
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#388: eps <- diff_logL/abs(logL_monitoring)
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#389: if (eps %>% is.nan()) {
#>     eps <- 1
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#394: t_i_2 <- Sys.time()
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#395: paste0("EM algorithm, step ", i, ": ", difftime(t_i_2, t_i_1, 
#>     units = "secs") %>% round(2), " seconds \n \n") %>% cat()
#> EM algorithm, step 1: 13.16 seconds 
#>  
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#402: paste0("Value of the likelihood: ", logL_monitoring %>% round(5), 
#>     " --- Convergence ratio = ", eps %>% round(5), "\n \n") %>% 
#>     cat()
#> Value of the likelihood: -242.84823 --- Convergence ratio = Inf
#>  
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#412: if (abs(eps) < cv_threshold) {
#>     cat("The EM algorithm successfully converged, training is completed.", 
#>         "\n \n")
#>     cv <- TRUE
#>     break
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#299: t_i_1 <- Sys.time()
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#302: post <- e_step(db = data, m_0 = m_0, kern_0 = kern_0, kern_i = kern_i, 
#>     hp_0 = hp_0, hp_i = hp_i, pen_diag = pen_diag)
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#313: if (fast_approx) {
#>     seq_loglikelihood <- logL_monitoring(hp_0 = hp_0, hp_i = hp_i, 
#>         db = data, m_0 = m_0, kern_0 = kern_0, kern_i = kern_i, 
#>         post_mean = post$mean, post_cov = post$cov, pen_diag = pen_diag)
#>     cv <- FALSE
#>     break
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#332: new_hp <- m_step(db = data, m_0 = m_0, kern_0 = kern_0, kern_i = kern_i, 
#>     old_hp_0 = hp_0, old_hp_i = hp_i, post_mean = post$mean, 
#>     post_cov = post$cov, common_hp = common_hp, pen_diag = pen_diag)
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#344: new_hp_0 <- new_hp$hp_0
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#345: new_hp_i <- new_hp$hp_i
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#348: if (any(is.na(new_hp_0)) | any(is.na(new_hp_i))) {
#>     warning(paste0("The M-step encountered an error at iteration : ", 
#>         i))
#>     warning("Training has stopped and hyper-parameters values from the ", 
#>         "last valid iteration are returned.")
#>     break
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#358: new_logL_monitoring <- logL_monitoring(hp_0 = new_hp_0, hp_i = new_hp_i, 
#>     db = data, m_0 = m_0, kern_0 = kern_0, kern_i = kern_i, post_mean = post$mean, 
#>     post_cov = post$cov, pen_diag = pen_diag)
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#370: diff_logL <- new_logL_monitoring - logL_monitoring
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#371: if (diff_logL %>% is.nan()) {
#>     diff_logL <- -Inf
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#375: if (diff_logL < 0) {
#>     warning("The likelihood descreased. Possible numerical issues.")
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#380: hp_0 <- new_hp_0
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#381: hp_i <- new_hp_i
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#382: logL_monitoring <- new_logL_monitoring
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#385: seq_loglikelihood <- c(seq_loglikelihood, logL_monitoring)
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#388: eps <- diff_logL/abs(logL_monitoring)
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#389: if (eps %>% is.nan()) {
#>     eps <- 1
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#394: t_i_2 <- Sys.time()
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#395: paste0("EM algorithm, step ", i, ": ", difftime(t_i_2, t_i_1, 
#>     units = "secs") %>% round(2), " seconds \n \n") %>% cat()
#> EM algorithm, step 2: 15.89 seconds 
#>  
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#402: paste0("Value of the likelihood: ", logL_monitoring %>% round(5), 
#>     " --- Convergence ratio = ", eps %>% round(5), "\n \n") %>% 
#>     cat()
#> Value of the likelihood: -231.94883 --- Convergence ratio = 0.04699
#>  
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#412: if (abs(eps) < cv_threshold) {
#>     cat("The EM algorithm successfully converged, training is completed.", 
#>         "\n \n")
#>     cv <- TRUE
#>     break
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#299: t_i_1 <- Sys.time()
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#302: post <- e_step(db = data, m_0 = m_0, kern_0 = kern_0, kern_i = kern_i, 
#>     hp_0 = hp_0, hp_i = hp_i, pen_diag = pen_diag)
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#313: if (fast_approx) {
#>     seq_loglikelihood <- logL_monitoring(hp_0 = hp_0, hp_i = hp_i, 
#>         db = data, m_0 = m_0, kern_0 = kern_0, kern_i = kern_i, 
#>         post_mean = post$mean, post_cov = post$cov, pen_diag = pen_diag)
#>     cv <- FALSE
#>     break
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#332: new_hp <- m_step(db = data, m_0 = m_0, kern_0 = kern_0, kern_i = kern_i, 
#>     old_hp_0 = hp_0, old_hp_i = hp_i, post_mean = post$mean, 
#>     post_cov = post$cov, common_hp = common_hp, pen_diag = pen_diag)
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#344: new_hp_0 <- new_hp$hp_0
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#345: new_hp_i <- new_hp$hp_i
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#348: if (any(is.na(new_hp_0)) | any(is.na(new_hp_i))) {
#>     warning(paste0("The M-step encountered an error at iteration : ", 
#>         i))
#>     warning("Training has stopped and hyper-parameters values from the ", 
#>         "last valid iteration are returned.")
#>     break
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#358: new_logL_monitoring <- logL_monitoring(hp_0 = new_hp_0, hp_i = new_hp_i, 
#>     db = data, m_0 = m_0, kern_0 = kern_0, kern_i = kern_i, post_mean = post$mean, 
#>     post_cov = post$cov, pen_diag = pen_diag)
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#370: diff_logL <- new_logL_monitoring - logL_monitoring
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#371: if (diff_logL %>% is.nan()) {
#>     diff_logL <- -Inf
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#375: if (diff_logL < 0) {
#>     warning("The likelihood descreased. Possible numerical issues.")
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#380: hp_0 <- new_hp_0
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#381: hp_i <- new_hp_i
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#382: logL_monitoring <- new_logL_monitoring
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#385: seq_loglikelihood <- c(seq_loglikelihood, logL_monitoring)
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#388: eps <- diff_logL/abs(logL_monitoring)
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#389: if (eps %>% is.nan()) {
#>     eps <- 1
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#394: t_i_2 <- Sys.time()
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#395: paste0("EM algorithm, step ", i, ": ", difftime(t_i_2, t_i_1, 
#>     units = "secs") %>% round(2), " seconds \n \n") %>% cat()
#> EM algorithm, step 3: 13.16 seconds 
#>  
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#402: paste0("Value of the likelihood: ", logL_monitoring %>% round(5), 
#>     " --- Convergence ratio = ", eps %>% round(5), "\n \n") %>% 
#>     cat()
#> Value of the likelihood: -231.6273 --- Convergence ratio = 0.00139
#>  
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#412: if (abs(eps) < cv_threshold) {
#>     cat("The EM algorithm successfully converged, training is completed.", 
#>         "\n \n")
#>     cv <- TRUE
#>     break
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#299: t_i_1 <- Sys.time()
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#302: post <- e_step(db = data, m_0 = m_0, kern_0 = kern_0, kern_i = kern_i, 
#>     hp_0 = hp_0, hp_i = hp_i, pen_diag = pen_diag)
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#313: if (fast_approx) {
#>     seq_loglikelihood <- logL_monitoring(hp_0 = hp_0, hp_i = hp_i, 
#>         db = data, m_0 = m_0, kern_0 = kern_0, kern_i = kern_i, 
#>         post_mean = post$mean, post_cov = post$cov, pen_diag = pen_diag)
#>     cv <- FALSE
#>     break
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#332: new_hp <- m_step(db = data, m_0 = m_0, kern_0 = kern_0, kern_i = kern_i, 
#>     old_hp_0 = hp_0, old_hp_i = hp_i, post_mean = post$mean, 
#>     post_cov = post$cov, common_hp = common_hp, pen_diag = pen_diag)
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#344: new_hp_0 <- new_hp$hp_0
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#345: new_hp_i <- new_hp$hp_i
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#348: if (any(is.na(new_hp_0)) | any(is.na(new_hp_i))) {
#>     warning(paste0("The M-step encountered an error at iteration : ", 
#>         i))
#>     warning("Training has stopped and hyper-parameters values from the ", 
#>         "last valid iteration are returned.")
#>     break
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#358: new_logL_monitoring <- logL_monitoring(hp_0 = new_hp_0, hp_i = new_hp_i, 
#>     db = data, m_0 = m_0, kern_0 = kern_0, kern_i = kern_i, post_mean = post$mean, 
#>     post_cov = post$cov, pen_diag = pen_diag)
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#370: diff_logL <- new_logL_monitoring - logL_monitoring
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#371: if (diff_logL %>% is.nan()) {
#>     diff_logL <- -Inf
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#375: if (diff_logL < 0) {
#>     warning("The likelihood descreased. Possible numerical issues.")
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#380: hp_0 <- new_hp_0
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#381: hp_i <- new_hp_i
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#382: logL_monitoring <- new_logL_monitoring
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#385: seq_loglikelihood <- c(seq_loglikelihood, logL_monitoring)
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#388: eps <- diff_logL/abs(logL_monitoring)
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#389: if (eps %>% is.nan()) {
#>     eps <- 1
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#394: t_i_2 <- Sys.time()
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#395: paste0("EM algorithm, step ", i, ": ", difftime(t_i_2, t_i_1, 
#>     units = "secs") %>% round(2), " seconds \n \n") %>% cat()
#> EM algorithm, step 4: 8.29 seconds 
#>  
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#402: paste0("Value of the likelihood: ", logL_monitoring %>% round(5), 
#>     " --- Convergence ratio = ", eps %>% round(5), "\n \n") %>% 
#>     cat()
#> Value of the likelihood: -231.61445 --- Convergence ratio = 6e-05
#>  
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#412: if (abs(eps) < cv_threshold) {
#>     cat("The EM algorithm successfully converged, training is completed.", 
#>         "\n \n")
#>     cv <- TRUE
#>     break
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#413: cat("The EM algorithm successfully converged, training is completed.", 
#>     "\n \n")
#> The EM algorithm successfully converged, training is completed. 
#>  
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#417: cv <- TRUE
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#418: break
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#422: if (!cv & (i == n_iter_max)) {
#>     warning("The EM algorithm has reached the maximum number of iterations ", 
#>         "before convergence, training might be sub-optimal \n \n")
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#430: if (!is.null(grid_inputs)) {
#>     cat("Start evaluating hyper-posterior distribution of the mean process", 
#>         "on the provided grid of inputs... \n \n")
#>     post <- hyperposterior(data = data, hp_0 = hp_0, hp_i = hp_i, 
#>         kern_0 = kern_0, kern_i = kern_i, prior_mean = prior_mean, 
#>         grid_inputs = grid_inputs, pen_diag = pen_diag)
#>     cat("Done!\n \n")
#> } else {
#>     post$pred <- tibble::tibble(Input = post$mean %>% dplyr::pull(.data$Input), 
#>         Mean = post$mean %>% dplyr::pull(.data$Output), Var = post$cov %>% 
#>             diag() %>% as.vector())
#> }
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#449: post$pred <- tibble::tibble(Input = post$mean %>% dplyr::pull(.data$Input), 
#>     Mean = post$mean %>% dplyr::pull(.data$Output), Var = post$cov %>% 
#>         diag() %>% as.vector())
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#457: fct_args <- list(data = data, prior_mean = prior_mean, ini_hp_0 = hp_0_ini, 
#>     ini_hp_i = hp_i_ini, kern_0 = kern_0, kern_i = kern_i, common_hp = common_hp, 
#>     grid_inputs = grid_inputs, pen_diag = pen_diag, n_iter_max = n_iter_max, 
#>     cv_threshold = cv_threshold)
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#471: t_2 <- Sys.time()
#> debug at C:/Users/user/Mon Drive/Travail/GitHub/MagmaClustR/R/training.R#473: list(hp_0 = hp_0, hp_i = hp_i, hyperpost = post, ini_args = fct_args, 
#>     seq_loglikelihood = seq_loglikelihood, converged = cv, training_time = difftime(t_2, 
#>         t_1, units = "secs")) %>% return()

pred_dim2  <- pred_magma(data = dim2_pred,
                         trained_model = model_dim2)
#> The hyper-posterior distribution of the mean process provided in 'hyperpost' argument isn't evaluated on the expected inputs.
#>  
#>  Start evaluating the hyper-posterior on the correct inputs...
#>  
#> The 'prior_mean' argument has not been specified. The hyper-prior mean function is thus set to be 0 everywhere.
#>  
#> Done!
#> 
```

<img src="man/figures/README-train_and_predict_Magma_in_2-D-1.png" width="100%" />
