
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MagmaClustR <img src="man/figures/logo.png" align="right" width="120" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/ArthurLeroy/MagmaClustR/workflows/R-CMD-check/badge.svg)](https://github.com/ArthurLeroy/MagmaClustR/actions)
[![codecov](https://codecov.io/gh/ArthurLeroy/MagmaClustR/branch/master/graph/badge.svg?token=KH7SOKOLKX)](https://app.codecov.io/gh/ArthurLeroy/MagmaClustR)
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
set.seed(28)
data_magma <- simu_db(M = 11, N = 10, common_input = FALSE)
## Split individuals into training and prediction sets, and define test points
magma_train <- data_magma %>% subset(ID %in% 1:10)
magma_pred <- data_magma %>% subset(ID == 11) %>% head(7)
magma_test <- data_magma %>% subset(ID == 11) %>% tail(3)

data_magma
#> # A tibble: 110 x 3
#>    ID    Output Input
#>    <chr>  <dbl> <dbl>
#>  1 1      -56.9  0   
#>  2 1      -57.4  0.7 
#>  3 1      -53.1  1   
#>  4 1      -62.8  1.5 
#>  5 1      -60.2  1.7 
#>  6 1      -81.4  5.2 
#>  7 1      -96.2  7   
#>  8 1     -102.   7.1 
#>  9 1      -96.4  8.05
#> 10 1      -86.4  9.65
#> # ... with 100 more rows
#> # i Use `print(n = ...)` to see more rows
```

As displayed above, any dataset processed in MagmaClustR should provide
columns named `ID`, `Input`, and `Output`. Any additional column would
be treated as a covariate (and thus define multi-dimensional inputs).

### Training and prediction with Magma

``` r
model <- train_magma(data = magma_train, common_hp = F)
#> The 'prior_mean' argument has not been specified. The hyper_prior mean function is thus set to be 0 everywhere.
#>  
#> The 'ini_hp_0' argument has not been specified. Random values of hyper-parameters for the mean process are used as initialisation.
#>  
#> The 'ini_hp_i' argument has not been specified. Random values of hyper-parameters for the individal processes are used as initialisation.
#>  
#> EM algorithm, step 1: 8.23 seconds 
#>  
#> Value of the likelihood: -434.32317 --- Convergence ratio = Inf
#>  
#> EM algorithm, step 2: 4.91 seconds 
#>  
#> Value of the likelihood: -401.89388 --- Convergence ratio = 0.08069
#>  
#> EM algorithm, step 3: 4.18 seconds 
#>  
#> Value of the likelihood: -400.50342 --- Convergence ratio = 0.00347
#>  
#> EM algorithm, step 4: 3.53 seconds 
#>  
#> Value of the likelihood: -400.21465 --- Convergence ratio = 0.00072
#>  
#> The EM algorithm successfully converged, training is completed. 
#> 

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
#> The 'hp' argument has not been specified. The 'train_gp()' function (with random initialisation) has been used to learn ML estimators for the hyper-parameters associated with the 'kern' argument.
#> 
```

<img src="man/figures/README-train_and_predict_Magma-1.png" width="80%" style="display: block; margin: auto;" />

Note that the `common_hp` and `grid_inputs` arguments are optional. They
respectively indicate that a specific set of hyper-parameters is trained
for each curve, and control the grid of values on which the prediction
is performed.

### Display the resulting predictions

Several other arguments are available in dedicated plotting functions to
offer extended options in the display of results. For instance, the GP
predictions can be represented as a heatmap of probabilities:

``` r
plot_gp(pred_gp = pred,
        data = magma_pred,
        data_train = magma_train,
        prior_mean = model$hyperpost$mean,
        heatmap = TRUE) 
```

<img src="man/figures/README-display_Magma-1.png" width="80%" style="display: block; margin: auto;" />

Additionally, it is also possible to create animated representations by
using functions that generate GIFs. For instance, below, the true
testing points have been represented as red dots and we can observe how
the prediction evolves as we add more data points to our prediction
dataset.

``` r

pred_gif  <- pred_gif(data = magma_pred,
                      trained_model = model,
                      grid_inputs = seq(0, 10, 0.01))
#>  => 1 => 2 => 3 => 4 => 5 => 6 => 7

plot_gif(pred_gp = pred_gif,
         data = magma_pred,
         data_train = magma_train,
         prior_mean = model$hyperpost$mean) + 
  ggplot2::geom_point(data = magma_test,
                      ggplot2::aes(x = Input, y = Output),
                      color = 'red', size = 2)
```

<img src="man/figures/README-gif_Magma-1.gif" width="80%" style="display: block; margin: auto;" />

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
#>  1 ID1-Clust1  -47.0  0.25
#>  2 ID1-Clust1  -39.5  0.8 
#>  3 ID1-Clust1  -29.7  2   
#>  4 ID1-Clust1  -38.7  4.2 
#>  5 ID1-Clust1  -43.4  4.6 
#>  6 ID1-Clust1  -53.6  6.2 
#>  7 ID1-Clust1  -57.5  6.75
#>  8 ID1-Clust1  -64.9  7.95
#>  9 ID1-Clust1  -53.3  8.85
#> 10 ID1-Clust1  -42.7  9.85
#> # ... with 110 more rows
#> # i Use `print(n = ...)` to see more rows
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
#> VEM algorithm, step 1: 44.73 seconds 
#>  
#> Value of the elbo: -454.05837 --- Convergence ratio = Inf
#>  
#> VEM algorithm, step 2: 14.69 seconds 
#>  
#> Value of the elbo: -409.05605 --- Convergence ratio = 0.11002
#>  
#> VEM algorithm, step 3: 14.16 seconds 
#>  
#> Value of the elbo: -408.61253 --- Convergence ratio = 0.00109
#>  
#> VEM algorithm, step 4: 14.46 seconds 
#>  
#> Value of the elbo: -408.20937 --- Convergence ratio = 0.00099
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
#> The mixture probability of the cluster K3 is 1. Therefore, the predictive distribution is Gaussian and the associated credible interval can be displayed.
```

<img src="man/figures/README-display_MagmaClust-1.png" width="80%" style="display: block; margin: auto;" />

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
#>  1 1      -47.0  0.25     -2   
#>  2 1      -39.5  0.8       1.94
#>  3 1      -29.7  2         4.64
#>  4 1      -38.7  4.2      -3.7 
#>  5 1      -43.4  4.6      -4.24
#>  6 1      -53.6  6.2       0.68
#>  7 1      -57.5  6.75      0.55
#>  8 1      -64.9  7.95     -4.38
#>  9 1      -53.3  8.85      1.74
#> 10 1      -42.7  9.85      4.14
#> # ... with 100 more rows
#> # i Use `print(n = ...)` to see more rows
```

### Training and prediction with Magma

``` r
model_dim2 <- train_magma(data = dim2_train)
#> The 'prior_mean' argument has not been specified. The hyper_prior mean function is thus set to be 0 everywhere.
#>  
#> The 'ini_hp_0' argument has not been specified. Random values of hyper-parameters for the mean process are used as initialisation.
#>  
#> The 'ini_hp_i' argument has not been specified. Random values of hyper-parameters for the individal processes are used as initialisation.
#>  
#> EM algorithm, step 1: 14.17 seconds 
#>  
#> Value of the likelihood: -547.76188 --- Convergence ratio = Inf
#>  
#> EM algorithm, step 2: 9.5 seconds 
#>  
#> Value of the likelihood: -494.5156 --- Convergence ratio = 0.10767
#>  
#> EM algorithm, step 3: 11.98 seconds 
#>  
#> Value of the likelihood: -466.4102 --- Convergence ratio = 0.06026
#>  
#> EM algorithm, step 4: 9.24 seconds 
#>  
#> Value of the likelihood: -449.17074 --- Convergence ratio = 0.03838
#>  
#> EM algorithm, step 5: 5 seconds 
#>  
#> Value of the likelihood: -435.574 --- Convergence ratio = 0.03122
#>  
#> EM algorithm, step 6: 5.4 seconds 
#>  
#> Value of the likelihood: -426.56255 --- Convergence ratio = 0.02113
#>  
#> EM algorithm, step 7: 5.17 seconds 
#>  
#> Value of the likelihood: -421.25636 --- Convergence ratio = 0.0126
#>  
#> EM algorithm, step 8: 5.15 seconds 
#>  
#> Value of the likelihood: -417.61643 --- Convergence ratio = 0.00872
#>  
#> EM algorithm, step 9: 5.11 seconds 
#>  
#> Value of the likelihood: -415.16894 --- Convergence ratio = 0.0059
#>  
#> EM algorithm, step 10: 4.99 seconds 
#>  
#> Value of the likelihood: -413.38881 --- Convergence ratio = 0.00431
#>  
#> EM algorithm, step 11: 4.8 seconds 
#>  
#> Value of the likelihood: -411.98102 --- Convergence ratio = 0.00342
#>  
#> EM algorithm, step 12: 4.25 seconds 
#>  
#> Value of the likelihood: -410.94035 --- Convergence ratio = 0.00253
#>  
#> EM algorithm, step 13: 4.25 seconds 
#>  
#> Value of the likelihood: -409.99426 --- Convergence ratio = 0.00231
#>  
#> EM algorithm, step 14: 4.25 seconds 
#>  
#> Value of the likelihood: -409.21728 --- Convergence ratio = 0.0019
#>  
#> EM algorithm, step 15: 4.35 seconds 
#>  
#> Value of the likelihood: -408.47281 --- Convergence ratio = 0.00182
#>  
#> EM algorithm, step 16: 4.33 seconds 
#>  
#> Value of the likelihood: -407.67549 --- Convergence ratio = 0.00196
#>  
#> EM algorithm, step 17: 4.42 seconds 
#>  
#> Value of the likelihood: -406.49068 --- Convergence ratio = 0.00291
#>  
#> EM algorithm, step 18: 4.67 seconds 
#>  
#> Value of the likelihood: -405.14487 --- Convergence ratio = 0.00332
#>  
#> EM algorithm, step 19: 4.22 seconds 
#>  
#> Value of the likelihood: -404.58873 --- Convergence ratio = 0.00137
#>  
#> EM algorithm, step 20: 4.31 seconds 
#>  
#> Value of the likelihood: -403.83303 --- Convergence ratio = 0.00187
#>  
#> EM algorithm, step 21: 4.51 seconds 
#>  
#> Value of the likelihood: -403.49586 --- Convergence ratio = 0.00084
#>  
#> The EM algorithm successfully converged, training is completed. 
#> 

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
#> Warning: Raster pixels are placed at uneven horizontal intervals and will be
#> shifted. Consider using geom_tile() instead.
#> Warning: Raster pixels are placed at uneven vertical intervals and will be
#> shifted. Consider using geom_tile() instead.
```

<img src="man/figures/README-train_and_predict_Magma_in_2-D-1.png" width="80%" style="display: block; margin: auto;" />
