
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
#>  1 1      -14.7  0   
#>  2 1      -15.7  0.7 
#>  3 1      -11.1  1   
#>  4 1      -19.3  1.5 
#>  5 1      -16.5  1.7 
#>  6 1      -33.5  5.2 
#>  7 1      -42.2  7   
#>  8 1      -48.9  7.1 
#>  9 1      -50.5  8.05
#> 10 1      -45.0  9.65
#> # ... with 100 more rows
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
#> EM algorithm, step 1: 1.25 seconds 
#>  
#> Value of the likelihood: -405.89184 --- Convergence ratio = Inf
#>  
#> EM algorithm, step 2: 0.87 seconds 
#>  
#> Value of the likelihood: -392.1544 --- Convergence ratio = 0.03503
#>  
#> EM algorithm, step 3: 0.86 seconds 
#>  
#> Value of the likelihood: -391.88908 --- Convergence ratio = 0.00068
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
#> VEM algorithm, step 1: 11.31 seconds 
#>  
#> Value of the elbo: -405.4094 --- Convergence ratio = Inf
#>  
#> VEM algorithm, step 2: 5.29 seconds 
#>  
#> Value of the elbo: -384.40578 --- Convergence ratio = 0.05464
#>  
#> VEM algorithm, step 3: 5.46 seconds 
#>  
#> Value of the elbo: -384.09759 --- Convergence ratio = 8e-04
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
#> The 'ini_hp_i' argument has not been specified. Random values of hyper-parameters for the individal processes are used as initialisation.
#>  
#> EM algorithm, step 1: 5.63 seconds 
#>  
#> Value of the likelihood: -243.32274 --- Convergence ratio = Inf
#>  
#> EM algorithm, step 2: 2.52 seconds 
#>  
#> Value of the likelihood: -232.15988 --- Convergence ratio = 0.04808
#>  
#> EM algorithm, step 3: 2.05 seconds 
#>  
#> Value of the likelihood: -231.72461 --- Convergence ratio = 0.00188
#>  
#> EM algorithm, step 4: 2.01 seconds 
#>  
#> Value of the likelihood: -231.68917 --- Convergence ratio = 0.00015
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
```

<img src="man/figures/README-train_and_predict_Magma_in_2-D-1.png" width="80%" style="display: block; margin: auto;" />
