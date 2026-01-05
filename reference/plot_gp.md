# Plot Magma or GP predictions

Display Magma or classic GP predictions. According to the dimension of
the inputs, the graph may be a mean curve + Credible Interval or a
heatmap of probabilities.

## Usage

``` r
plot_gp(
  pred_gp,
  x_input = NULL,
  data = NULL,
  data_train = NULL,
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

plot_magma(
  pred_gp,
  x_input = NULL,
  data = NULL,
  data_train = NULL,
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

- pred_gp:

  A tibble or data frame, typically coming from
  [`pred_magma`](https://arthurleroy.github.io/MagmaClustR/reference/pred_magma.md)
  or
  [`pred_gp`](https://arthurleroy.github.io/MagmaClustR/reference/pred_gp.md)
  functions. Required columns: 'Input', 'Mean', 'Var'. Additional
  covariate columns may be present in case of multi-dimensional inputs.

- x_input:

  A vector of character strings, indicating which input should be
  displayed. If NULL (default) the 'Input' column is used for the
  x-axis. If providing a 2-dimensional vector, the corresponding columns
  are used for the x-axis and y-axis.

- data:

  (Optional) A tibble or data frame. Required columns: 'Input',
  'Output'. Additional columns for covariates can be specified. This
  argument corresponds to the raw data on which the prediction has been
  performed.

- data_train:

  (Optional) A tibble or data frame, containing the training data of the
  Magma model. The data set should have the same format as the `data`
  argument with an additional required column 'ID' for identifying the
  different individuals/tasks. If provided, those data are displayed as
  backward colourful points (each colour corresponding to one
  individual/task).

- prior_mean:

  (Optional) A tibble or a data frame, containing the 'Input' and
  associated 'Output' prior mean parameter of the GP prediction.

- y_grid:

  A vector, indicating the grid of values on the y-axis for which
  probabilities should be computed for heatmaps of 1-dimensional
  predictions. If NULL (default), a vector of length 50 is defined,
  ranging between the min and max 'Output' values contained in
  `pred_gp`.

- heatmap:

  A logical value indicating whether the GP prediction should be
  represented as a heatmap of probabilities for 1-dimensional inputs. If
  FALSE (default), the mean curve and associated Credible Interval are
  displayed.

- samples:

  A logical value indicating whether the GP prediction should be
  represented as a collection of samples drawn from the posterior. If
  FALSE (default), the mean curve and associated Credible Interval are
  displayed.

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

Visualisation of a Magma or GP prediction (optional: display data
points, training data points and the prior mean function). For 1-D
inputs, the prediction is represented as a mean curve and its associated
95% Credible Interval, as a collection of samples drawn from the
posterior if `samples` = TRUE, or as a heatmap of probabilities if
`heatmap` = TRUE. For 2-D inputs, the prediction is represented as a
heatmap, where each couple of inputs on the x-axis and y-axis are
associated with a gradient of colours for the posterior mean values,
whereas the uncertainty is indicated by the transparency (the narrower
is the Credible Interval, the more opaque is the associated colour, and
vice versa)

## Examples

``` r
TRUE
#> [1] TRUE
```
