# Create a GIF of Magma or GP predictions

Create a GIF animation displaying how Magma or classic GP predictions
evolve and improve when the number of data points increase.

## Usage

``` r
plot_gif(
  pred_gp,
  x_input = NULL,
  data = NULL,
  data_train = NULL,
  prior_mean = NULL,
  y_grid = NULL,
  heatmap = FALSE,
  prob_CI = 0.95,
  size_data = 3,
  size_data_train = 1,
  alpha_data_train = 0.5,
  export_gif = FALSE,
  path = "gif_gp.gif",
  ...
)
```

## Arguments

- pred_gp:

  A tibble, typically coming from the
  [`pred_gif`](https://arthurleroy.github.io/MagmaClustR/reference/pred_gif.md)
  function. Required columns: 'Input', 'Mean', 'Var' and 'Index'.

- x_input:

  A vector of character strings, indicating which input should be
  displayed. If NULL(default) the 'Input' column is used for the x-axis.
  If providing a 2-dimensional vector, the corresponding columns are
  used for the x-axis and y-axis.

- data:

  (Optional) A tibble or data frame. Required columns: 'Input',
  'Output'. Additional columns for covariates can be specified. The
  'Input' column should define the variable that is used as reference
  for the observations (e.g. time for longitudinal data). The 'Output'
  column specifies the observed values (the response variable). The data
  frame can also provide as many covariates as desired, with no
  constraints on the column names. These covariates are additional
  inputs (explanatory variables) of the models that are also observed at
  each reference 'Input'.

- data_train:

  (Optional) A tibble or data frame, containing the training data of the
  Magma model. The data set should have the same format as the `data`
  argument with an additional column 'ID' for identifying the different
  individuals/tasks. If provided, those data are displayed as backward
  colourful points (each colour corresponding to one individual/task).

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
  FALSE (default), the mean curve and associated 95% CI are displayed.

- prob_CI:

  A number between 0 and 1 (default is 0.95), indicating the level of
  the Credible Interval associated with the posterior mean curve.

- size_data:

  A number, controlling the size of the `data` points.

- size_data_train:

  A number, controlling the size of the `data_train` points.

- alpha_data_train:

  A number, between 0 and 1, controlling transparency of the
  `data_train` points.

- export_gif:

  A logical value indicating whether the animation should be exported as
  a .gif file.

- path:

  A character string defining the path where the GIF file should be
  exported.

- ...:

  Any additional parameters that can be passed to the function
  [`transition_states`](https://gganimate.com/reference/transition_states.html)
  from the `gganimate` package.

## Value

Visualisation of a Magma or GP prediction (optional: display data
points, training data points and the prior mean function), where data
points are added sequentially for visualising changes in prediction as
information increases.

## Examples

``` r
TRUE
#> [1] TRUE
```
