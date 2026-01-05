# Plot MagmaClust predictions

Display MagmaClust predictions. According to the dimension of the
inputs, the graph may be a mean curve (dim inputs = 1) or a heatmap (dim
inputs = 2) of probabilities. Moreover, MagmaClust can provide credible
intervals only by visualising cluster-specific predictions (e.g. for the
most probable cluster). When visualising the full mixture-of-GPs
prediction, which can be multimodal, the user should choose between the
simple mean function or the full heatmap of probabilities (more
informative but slower).

## Usage

``` r
plot_magmaclust(
  pred_clust,
  cluster = "all",
  x_input = NULL,
  data = NULL,
  data_train = NULL,
  col_clust = FALSE,
  prior_mean = NULL,
  y_grid = NULL,
  heatmap = FALSE,
  samples = FALSE,
  nb_samples = 50,
  plot_mean = TRUE,
  alpha_samples = 0.3,
  prob_CI = 0.95,
  size_data = 3,
  size_data_train = 1,
  alpha_data_train = 0.5
)
```

## Arguments

- pred_clust:

  A list of predictions, typically coming from
  [`pred_magmaclust`](https://arthurleroy.github.io/MagmaClustR/reference/pred_magmaclust.md).
  Required elements: `pred`, `mixture`, `mixture_pred`.

- cluster:

  A character string, indicating which cluster to plot from. If 'all'
  (default) the mixture of GPs prediction is displayed as a mean curve
  (1-D inputs) or a mean heatmap (2-D inputs). Alternatively, if the
  name of one cluster is provided, the classic mean curve + credible
  interval is displayed (1-D inputs), or a heatmap with colour gradient
  for the mean and transparency gradient for the Credible Interval (2-D
  inputs).

- x_input:

  A vector of character strings, indicating which input should be
  displayed. If NULL (default) the 'Input' column is used for the
  x-axis. If providing a 2-dimensional vector, the corresponding columns
  are used for the x-axis and y-axis.

- data:

  (Optional) A tibble or data frame. Required columns: `Input` ,
  `Output`. Additional columns for covariates can be specified. This
  argument corresponds to the raw data on which the prediction has been
  performed.

- data_train:

  (Optional) A tibble or data frame, containing the training data of the
  MagmaClust model. The data set should have the same format as the
  `data` argument with an additional required column `ID` for
  identifying the different individuals/tasks. If provided, those data
  are displayed as backward colourful points (each colour corresponding
  to one individual or a cluster, see `col_clust` below).

- col_clust:

  A boolean indicating whether backward points are coloured according to
  the individuals or to their most probable cluster. If one wants to
  colour by clusters, a column `Cluster` shall be present in
  `data_train`. We advise to use
  [`data_allocate_cluster`](https://arthurleroy.github.io/MagmaClustR/reference/data_allocate_cluster.md)
  for automatically creating a well-formatted dataset from a trained
  MagmaClust model.

- prior_mean:

  (Optional) A list providing, for each cluster, a tibble containing
  prior mean parameters of the prediction. This argument typically comes
  as an outcome `hyperpost$mean`, available through the
  [`train_magmaclust`](https://arthurleroy.github.io/MagmaClustR/reference/train_magmaclust.md),
  [`pred_magmaclust`](https://arthurleroy.github.io/MagmaClustR/reference/pred_magmaclust.md)
  functions.

- y_grid:

  A vector, indicating the grid of values on the y-axis for which
  probabilities should be computed for heatmaps of 1-dimensional
  predictions. If NULL (default), a vector of length 50 is defined,
  ranging between the min and max 'Output' values contained in `pred`.

- heatmap:

  A logical value indicating whether the GP mixture should be
  represented as a heatmap of probabilities for 1-dimensional inputs. If
  FALSE (default), the mean curve (and associated Credible Interval if
  available) are displayed.

- samples:

  A logical value indicating whether the GP mixture should be
  represented as a collection of samples drawn from the posterior. If
  FALSE (default), the mean curve (and associated Credible Interval if
  available) are displayed.

- nb_samples:

  A number, indicating the number of samples to be drawn from the
  predictive posterior distribution. For two-dimensional graphs, only
  one sample can be displayed.

- plot_mean:

  A logical value, indicating whether the mean prediction should be
  displayed on the graph when `samples = TRUE`.

- alpha_samples:

  A number, controlling transparency of the sample curves.

- prob_CI:

  A number between 0 and 1 (default is 0.95), indicating the level of
  the Credible Interval associated with the posterior mean curve. If
  this this argument is set to 1, the Credible Interval is not
  displayed.

- size_data:

  A number, controlling the size of the `data` points.

- size_data_train:

  A number, controlling the size of the `data_train` points.

- alpha_data_train:

  A number, between 0 and 1, controlling transparency of the
  `data_train` points.

## Value

Visualisation of a MagmaClust prediction (optional: display data points,
training data points and the prior mean functions). For 1-D inputs, the
prediction is represented as a mean curve (and its associated 95%
Credible Interval for cluster-specific predictions), or as a heatmap of
probabilities if `heatmap` = TRUE. In the case of MagmaClust, the
heatmap representation should be preferred for clarity, although the
default display remains mean curve for quicker execution. For 2-D
inputs, the prediction is represented as a heatmap, where each couple of
inputs on the x-axis and y-axis are associated with a gradient of
colours for the posterior mean values, whereas the uncertainty is
indicated by the transparency (the narrower is the Credible Interval,
the more opaque is the associated colour, and vice versa). As for 1-D
inputs, Credible Interval information is only available for
cluster-specific predictions.

## Examples

``` r
TRUE
#> [1] TRUE
```
