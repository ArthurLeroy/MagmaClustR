#' Plot smoothed curves of raw data
#'
#' Display raw data under the Magma format as smoothed curves.
#'
#' @param data A data frame or tibble with format : ID, Input, Output.
#' @param cluster A boolean indicating whether data should be coloured by
#'   cluster. Requires a column named 'Cluster'.
#' @param legend A boolean indicating whether the legend should be displayed.
#'
#' @return Graph of smoothed curves of raw data.
#'
#' @examples
#' TRUE
plot_db <- function(data, cluster = F, legend = F) {
  ## Convert Cluster into factors for a better display
  data$ID = as.factor(data$ID)
  if (cluster) {
    ## Add a dummy column 'Cluster' if absent
    if(!('Cluster' %in% names(data))){data$Cluster = 1}
    ## Convert Cluster into factors for a better display
    data$Cluster = as.factor(data$Cluster)

    gg = ggplot2::ggplot(data) +
      ggplot2::geom_smooth(ggplot2::aes(
        x = .data$Input,
        y = .data$Output,
        group = .data$ID,
        color = .data$Cluster
      ),
      se = F
      ) +
      ggplot2::geom_point(ggplot2::aes(
        x = .data$Input,
        y = .data$Output,
        group = .data$ID,
        color = .data$Cluster
      )) +
      ggplot2::theme_classic()
  } else {
    gg = ggplot2::ggplot(data) +
      ggplot2::geom_smooth(ggplot2::aes(
        x =.data$Input,
        y = .data$Output,
        color = .data$ID
        ),
        se = F
        ) +
      ggplot2::geom_point(ggplot2::aes(
        x = .data$Input,
        y = .data$Output,
        color = .data$ID)) +
      ggplot2::theme_classic()
  }
  if(!legend)
  {
    gg = gg + ggplot2::guides(col = 'none')
  }
  return(gg)
}

#' Plot Magma or GP predictions
#'
#' Display Magma or classic GP predictions. According to the dimension of the
#' inputs, the graph may be a mean curve + Credible Interval or a heatmap of
#' probabilities.
#'
#' @param pred_gp A tibble or data frame, typically coming from
#'    \code{\link{pred_magma}} or \code{\link{pred_gp}} functions. Required
#'    columns: 'Input', 'Mean', 'Var'. Additional covariate columns may be
#'    present in case of multi-dimensional inputs.
#' @param x_input A vector of character strings, indicating which input should
#'    be displayed. If NULL(default) the 'Input' column is used for the x-axis.
#'    If providing a 2-dimensional vector, the corresponding columns are used
#'    for the x-axis and y-axis.
#' @param data (Optional) A tibble or data frame. Required columns: 'Input',
#'    'Output'. Additional columns for covariates can be specified.
#'    The 'Input' column should define the variable that is used as
#'    reference for the observations (e.g. time for longitudinal data). The
#'    'Output' column specifies the observed values (the response
#'    variable). The data frame can also provide as many covariates as desired,
#'    with no constraints on the column names. These covariates are additional
#'    inputs (explanatory variables) of the models that are also observed at
#'    each reference 'Input'.
#' @param data_train (Optional) A tibble or data frame, containing the training
#'    data of the Magma model. The data set should have the same format as the
#'    \code{data} argument with an additional column 'ID' for identifying the
#'    different individuals/tasks. If provided, those data are displayed as
#'    backward colourful points (each colour corresponding to one
#'    individual/task).
#' @param prior_mean (Optional) A tibble or a data frame, containing the 'Input'
#'    and associated 'Output' prior mean parameter of the GP prediction.
#' @param y_grid A vector, indicating the grid of values on the y-axis for which
#'    probabilities should be computed for heatmaps of 1-dimensional
#'    predictions. If NULL (default), a vector of length 50 is defined, ranging
#'    between the min and max 'Output' values contained in \code{pred_gp}.
#' @param heatmap A logical value indicating whether the GP prediction should be
#'    represented as a heatmap of probabilities for 1-dimensional inputs. If
#'    FALSE (default), the mean curve and associated 95%CI are displayed.
#' @param prob_CI A number between 0 and 1 (default is 0.95), indicating the
#'    level of the Credible Interval associated with the posterior mean curve.
#'
#' @return Visualisation of a Magma or GP prediction (optional: display data
#'    points, training data points and the prior mean function). For 1-D inputs,
#'    the prediction is represented as a mean curve and its associated 95%
#'    Credible Interval, or as a heatmap of probabilities if \code{heatmap} =
#'    TRUE. For 2-D inputs, the prediction is represented as a heatmap, where
#'    each couple of inputs on the x-axis and y-axis are associated
#'    with a gradient of colours for the posterior mean values, whereas the
#'    uncertainty is indicated by the transparency (the narrower is the 95%CI,
#'    the more opaque is the associated colour, and vice versa)
#' @export
#'
#' @examples
#' \dontrun{
#' ## 1-dimensional example
#' db <- simu_db(M = 1, covariate = FALSE)
#' pred_gp(db) %>%
#'   plot_gp(data = db)
#'
#' ## 2-dimensional example
#' db_2D <- simu_db(M = 1, covariate = TRUE)
#' pred_gp(db_2D) %>%
#'   plot_gp(data = db_2D)
#' }
plot_gp <- function(pred_gp,
                    x_input = NULL,
                    data = NULL,
                    data_train = NULL,
                    prior_mean = NULL,
                    y_grid = NULL,
                    heatmap = F,
                    prob_CI = 0.95) {
  ## Compute the quantile of the desired Credible Interval
  quant_ci <- stats::qnorm((1 + prob_CI) / 2)

  ## Check whether 'pred_gp' has a correct format
  if (pred_gp %>% is.data.frame()) {
    pred <- pred_gp
  } else if (is.list(pred_gp) &
    tryCatch(
      is.data.frame(pred_gp$pred),
      error = function(e) {
        FALSE
      }
    )) {
    pred <- pred_gp$pred
    ## Check whether the hyper-posterior distribution is provided and extract
    if (tryCatch(
      is.list(pred_gp$hyperpost),
      error = function(e) {
        FALSE
      }
    ) &
      tryCatch(
        is.data.frame(pred_gp$hyperpost$mean),
        error = function(e) {
          0
        }
      ) &
      is.null(prior_mean)) {
      prior_mean <- pred_gp$hyperpost$mean
    }
  } else {
    stop(
      "The 'pred_gp' argument should either be a list containing the 'pred' ",
      "element or a data frame. Please read ?plot_gp(), and use ",
      "pred_gp() or pred_magma() for making predictions under a correct format."
    )
  }

  ## Remove the 'Index' column if the prediction comes from 'pred_gif()'
  if (any("Index" %in% names(pred))) {
    index <- pred %>% dplyr::pull(.data$Index)
    pred <- pred %>% dplyr::select(-.data$Index)
  } else {
    index <- NULL
  }

  ## Rename the 'Output' column for enabling plot of the mean process in Magma
  if ("Output" %in% names(pred)) {
    pred <- pred %>% dplyr::rename("Mean" = .data$Output)
  }
  ## Get the inputs that should be used
  if (x_input %>% is.null()) {
    inputs <- pred %>% dplyr::select(-c(.data$Mean, .data$Var))
  } else {
    inputs <- pred[x_input]
  }

  ## Format the tibble for displaying the Credible Intervals
  pred <- pred %>%
    dplyr::mutate("CI_inf" = .data$Mean - quant_ci * sqrt(.data$Var)) %>%
    dplyr::mutate("CI_sup" = .data$Mean + quant_ci * sqrt(.data$Var)) %>%
    dplyr::mutate("CI_Width" = .data$CI_sup - .data$CI_inf)
  ## Display a heatmap if inputs are 2D
  if (ncol(inputs) == 2) {
    ## Add the 'Index' column if the prediction comes from 'pred_gif()'
    if (!is.null(index)) {
      pred <- pred %>% dplyr::mutate("Index" = index)
    }

    gg <- ggplot2::ggplot() +
      ggplot2::geom_raster(
        data = pred,
        ggplot2::aes_string(
          x = names(inputs)[1],
          y = names(inputs)[2],
          fill = "Mean",
          alpha = "CI_Width"
        ),
        interpolate = TRUE
      ) +
      ggplot2::scale_fill_distiller(palette = "RdPu", trans = "reverse") +
      ggplot2::scale_alpha_continuous(range = c(0.1, 1), trans = "reverse")

    if (!is.null(data)) {
      ## Round the 'Output' values to reduce size of labels on the graph
      data <-
        data %>% dplyr::mutate(Output = round(.data$Output, 1))

      gg <- gg + ggplot2::geom_label(
        data = data,
        ggplot2::aes_string(
          x = names(inputs)[1],
          y = names(inputs)[2],
          label = "Output",
          fill = "Output"
        ),
        size = 3
      )
    }
  } else {
    ## Check the dimension of the inputs a propose an adequate representation
    if (ncol(inputs) == 1) {
      if ((dplyr::n_distinct(inputs) != nrow(inputs)) & is.null(index)) {
        warning(
          "Some values on the x-axis appear multiple times, probably ",
          "resulting in an incorrect graphical representation. Please ",
          "consider recomputing predictions for more adequate inputs. "
        )
      }
    } else {
      warning(
        "Impossible to display inputs with dimensions greater than 2. The ",
        "graph then simply uses 'Input' as x_axis and 'Output' as y-axis. "
      )
      inputs <- inputs %>% dplyr::select(.data$Input)
    }

    ## Display a 'heatmap' if the argument is TRUE
    if (heatmap) {
      if (is.null(y_grid)) {
        y_grid <- seq(
          min(pred$Mean) - quant_ci * sqrt(max(pred$Var)),
          max(pred$Mean) + quant_ci * sqrt(max(pred$Var)),
          length.out = 500
        )
      }
      ## Define the columns needed to compute a prediction for all y-axis values
      col_to_nest <- c(names(inputs)[1], "Mean", "Var")

      ## Add the 'Index' column if the prediction comes from 'pred_gif()'
      if (!is.null(index)) {
        pred <- pred %>% dplyr::mutate("Index" = index)
        col_to_nest <- c(col_to_nest, "Index")
      }

      db_heat <- pred %>%
        tidyr::expand(tidyr::nesting_(col_to_nest),
          "Ygrid" = y_grid
        ) %>%
        dplyr::mutate("Proba" = 2 * (1 - stats::pnorm(abs((.data$Ygrid - .data$Mean) / sqrt(.data$Var)))))

      gg <- ggplot2::ggplot() +
        ggplot2::geom_raster(
          data = db_heat,
          ggplot2::aes_string(
            x = names(inputs)[1],
            y = "Ygrid",
            fill = "Proba"
          ),
          interpolate = TRUE
        ) +
        ggplot2::scale_fill_gradientn(
          colours = c(
            "white",
            "#FDE0DD",
            "#FCC5C0",
            "#FA9FB5",
            "#F768A1",
            "#DD3497",
            "#AE017E",
            "#7A0177"
          )
        ) +
        ggplot2::labs(fill = "Proba CI") +
        ggplot2::ylab("Output")
    } else {
      ## Display a classic curve otherwise
      ## Add the 'Index' column if the prediction comes from 'pred_gif()'
      if (!is.null(index)) {
        pred <- pred %>% dplyr::mutate("Index" = index)
      }

      gg <- ggplot2::ggplot() +
        ggplot2::geom_line(
          data = pred,
          ggplot2::aes_string(x = names(inputs)[1], y = "Mean"),
          color = "#DB15C1"
        ) +
        ggplot2::geom_ribbon(
          data = pred,
          ggplot2::aes_string(
            x = names(inputs)[1],
            ymin = "CI_inf",
            ymax = "CI_sup"
          ),
          alpha = 0.2,
          fill = "#FA9FB5"
        ) +
        ggplot2::ylab("Output")
    }

    ## Display the training data if provided
    if (!is.null(data_train)) {
      gg <- gg + ggplot2::geom_point(
        data = data_train,
        ggplot2::aes_string(
          x = names(inputs)[1],
          y = "Output",
          col = "ID"
        ),
        size = 0.5,
        alpha = 0.5
      ) +
        ggplot2::guides(color = FALSE)
    }
    ## Display the observed data if provided
    if (!is.null(data)) {
      gg <- gg + ggplot2::geom_point(
        data = data,
        ggplot2::aes_string(x = names(inputs)[1], y = "Output"),
        size = 2,
        shape = 20
      )
    }

    ## Display the (hyper-)prior mean process if provided
    if (!is.null(prior_mean)) {
      if (names(inputs)[1] %in% names(prior_mean)) {
        gg <- gg +
          ggplot2::geom_line(
            data = prior_mean,
            ggplot2::aes_string(x = names(inputs)[1], y = "Output"),
            linetype = "dashed"
          )
      } else {
        warning(
          "The ",
          names(inputs)[1],
          " column does not exist in the ",
          "'prior_mean' argument. The mean function cannot be displayed."
        )
      }
    }
  }

  (gg + ggplot2::theme_classic()) %>%
    return()
}

#' Display Realisation From Posterior GP
#'
#' A realisation of a posterior GP distribution is drawn and displayed.
#' According to the dimension of the inputs, the graph may be a curve or a
#' heatmap.
#'
#' @param pred_gp A tibble or data frame, typically coming from
#'    \code{\link{pred_magma}} or \code{\link{pred_gp}} functions. Required
#'    columns: 'Input', 'Mean', 'Var'. Additional covariate columns may be
#'    present in case of multi-dimensional inputs.
#' @param x_input A vector of character strings, indicating which input should
#'    be displayed. If NULL(default) the 'Input' column is used for the x-axis.
#'    If providing a 2-dimensional vector, the corresponding columns are used
#'    for the x-axis and y-axis.
#' @param data (Optional) A tibble or data frame, containing the data used in
#'    the GP prediction.
#' @param data_train (Optional) A tibble or data frame, containing the training
#'    data of the Magma model. The data set should have the same format as the
#'    \code{data} argument with an additional column 'ID' for identifying the
#'    different individuals/tasks. If provided, those data are displayed as
#'    backward colourful points (each colour corresponding to one
#'    individual/task).
#' @param prior_mean (Optional) A tibble or a data frame, containing the 'Input'
#'    and associated 'Output' prior mean parameter of the GP prediction.
#'
#' @return Draw and visualise from a posterior distribution from Magma or GP
#'    prediction (optional: display data points, training data points and the
#'    prior mean function).
#' @export
#'
#' @examples
#' \dontrun{
#' ## 1-dimensional example
#' db <- simu_db(M = 1, covariate = FALSE)
#' hp <- train_gp(db)
#' pred_gp(db, get_full_cov = TRUE, plot = TRUE) %>%
#'   sample_gp(data = db)
#'
#' ## 2-dimensional example
#' db_2D <- simu_db(M = 1, covariate = TRUE)
#' pred_gp(db_2D, get_full_cov = TRUE, plot = FALSE) %>%
#'   sample_gp(data = db_2D)
#' }
sample_gp <- function(pred_gp,
                      x_input = NULL,
                      data = NULL,
                      data_train = NULL,
                      prior_mean = NULL) {
  if (is.data.frame(pred_gp) | !is.list(pred_gp)) {
    stop(
      "The 'pred_gp' argument should be a list containing 'pred' and 'cov' ",
      "elements. Consider re-running the prediction function using the ",
      "argument 'get_full_cov' = TRUE."
    )
  }
  ## Get the inputs that should be used
  if (x_input %>% is.null()) {
    inputs <- pred_gp$pred %>% dplyr::select(-c(.data$Mean, .data$Var))
  } else {
    inputs <- pred_gp$pred[x_input]
  }

  ## Extract the predictions for further displaying
  input <- pred_gp$pred %>% dplyr::pull(.data$Input)
  mean <- pred_gp$pred %>% dplyr::pull(.data$Mean)
  cov <- pred_gp$cov

  sample <-
    tibble::tibble("Output" = mvtnorm::rmvnorm(1, mean, cov) %>% as.vector()) %>%
    dplyr::bind_cols(inputs)

  ## Display a heatmap if inputs are 2D
  if (ncol(inputs) == 2) {
    gg <- ggplot2::ggplot() +
      ggplot2::geom_raster(
        data = sample,
        ggplot2::aes_string(
          x = names(inputs)[1],
          y = names(inputs)[2],
          fill = "Output"
        ),
        interpolate = TRUE
      ) +
      ggplot2::scale_fill_distiller(palette = "RdPu", trans = "reverse")

    if (!is.null(data)) {
      ## Round the 'Output' values to reduce size of labels on the graph
      data <-
        data %>% dplyr::mutate(Output = round(.data$Output, 1))

      gg <- gg + ggplot2::geom_label(
        data = data,
        ggplot2::aes_string(
          x = names(inputs)[1],
          y = names(inputs)[2],
          label = "Output",
          fill = "Output"
        ),
        size = 3
      )
    }
  } else if (ncol(inputs) == 1) {
    if (dplyr::n_distinct(inputs) != nrow(inputs)) {
      warning(
        "Some values on the x-axis appear multiple times, probably resulting ",
        "in an incorrect graphical representation. Please consider ",
        "recomputing predictions for more adequate inputs. "
      )
    }
    gg <- ggplot2::ggplot() +
      ggplot2::geom_line(
        data = sample,
        ggplot2::aes_string(
          x = names(inputs)[1],
          y = "Output"
        ),
        color = "#DB15C1"
      )
    ## Display the training data if provided
    if (!is.null(data_train)) {
      gg <- gg + ggplot2::geom_point(
        data = data_train,
        ggplot2::aes_string(
          x = names(inputs)[1],
          y = "Output",
          col = "ID"
        ),
        size = 0.5,
        alpha = 0.5
      ) +
        ggplot2::guides(color = FALSE)
    }
    ## Display the observed data if provided
    if (!is.null(data)) {
      gg <- gg + ggplot2::geom_point(
        data = data,
        ggplot2::aes_string(
          x = names(inputs)[1],
          y = "Output"
        ),
        size = 2,
        shape = 18
      )
    }
  } else {
    warning(
      "Impossible to display inputs with dimensions greater than 2. The graph ",
      "then simply uses 'Input' as x_axis and 'Output' as y-axis. "
    )
    sample <- tibble::tibble(
      "Input" = input,
      "Output" = mvtnorm::rmvnorm(1, mean, cov) %>% as.vector()
    )
    gg <- ggplot2::ggplot() +
      ggplot2::geom_line(
        data = sample,
        ggplot2::aes(x = .data$Input, y = .data$Output),
        color = "#DB15C1"
      )

    ## Display the training data if provided
    if (!is.null(data_train)) {
      gg <- gg + ggplot2::geom_point(
        data = data_train,
        ggplot2::aes(
          x = .data$Input,
          y = .data$Output,
          col = .data$ID
        ),
        size = 0.5,
        alpha = 0.5
      ) +
        ggplot2::guides(color = FALSE)
    }
    ## Display the observed data if provided
    if (!is.null(data)) {
      gg <- gg +
        ggplot2::geom_point(
          data = data,
          ggplot2::aes(x = .data$Input, y = .data$Output),
          size = 2,
          shape = 18
        )
    }
    ## Display the (hyper-)prior mean process if provided
    if (!is.null(prior_mean)) {
      if (names(inputs)[1] %in% names(prior_mean)) {
        gg <- gg +
          ggplot2::geom_line(
            data = prior_mean,
            ggplot2::aes_string(x = names(inputs)[1], y = "Output"),
            linetype = "dashed"
          )
      } else {
        warning(
          "The ",
          names(inputs)[1],
          " column does not exist in the ",
          "'prior_mean' argument. The mean function cannot be displayed."
        )
      }
    }
  }

  (gg + ggplot2::theme_classic()) %>%
    return()
}

#' Create a GIF of Magma or GP predictions
#'
#' Create a GIF animation displaying how Magma or classic GP
#' predictions evolve and improve when the number of data points increase.
#'
#' @param pred_gp A tibble, typically coming from the \code{\link{pred_gif}}
#'    function. Required columns: 'Input', 'Mean', 'Var' and 'Index'.
#' @param x_input A vector of character strings, indicating which input should
#'    be displayed. If NULL(default) the 'Input' column is used for the x-axis.
#'    If providing a 2-dimensional vector, the corresponding columns are used
#'    for the x-axis and y-axis.
#' @param data (Optional) A tibble or data frame. Required columns: 'Input',
#'    'Output'. Additional columns for covariates can be specified.
#'    The 'Input' column should define the variable that is used as
#'    reference for the observations (e.g. time for longitudinal data). The
#'    'Output' column specifies the observed values (the response
#'    variable). The data frame can also provide as many covariates as desired,
#'    with no constraints on the column names. These covariates are additional
#'    inputs (explanatory variables) of the models that are also observed at
#'    each reference 'Input'.
#' @param data_train (Optional) A tibble or data frame, containing the training
#'    data of the Magma model. The data set should have the same format as the
#'    \code{data} argument with an additional column 'ID' for identifying the
#'    different individuals/tasks. If provided, those data are displayed as
#'    backward colourful points (each colour corresponding to one
#'    individual/task).
#' @param prior_mean (Optional) A tibble or a data frame, containing the 'Input'
#'    and associated 'Output' prior mean parameter of the GP prediction.
#' @param y_grid A vector, indicating the grid of values on the y-axis for which
#'    probabilities should be computed for heatmaps of 1-dimensional
#'    predictions. If NULL (default), a vector of length 50 is defined, ranging
#'    between the min and max 'Output' values contained in \code{pred_gp}.
#' @param heatmap A logical value indicating whether the GP prediction should be
#'    represented as a heatmap of probabilities for 1-dimensional inputs. If
#'    FALSE (default), the mean curve and associated 95% CI are displayed.
#' @param prob_CI A number between 0 and 1 (default is 0.95), indicating the
#'    level of the Credible Interval associated with the posterior mean curve.
#' @param export_gif A logical value indicating whether the animation should
#'    be exported as a .gif file.
#' @param path A character string defining the path where the GIF file should be
#'    exported.
#' @param ... Any additional parameters that can be passed to the function
#'    \code{\link[gganimate]{transition_states}} from the \code{gganimate}
#'    package.
#'
#' @return Visualisation of a Magma or GP prediction (optional: display data
#'    points, training data points and the prior mean function), where data
#'    points are added sequentially for visualising changes in prediction as
#'    information increases.
#' @export
#'
#' @examples
#' \dontrun{
#' ## 2-dimensional example
#' db <- simu_db(M = 1, covariate = FALSE)
#' hp <- train_gp(db)
#' pred_gif(db, hp = hp) %>%
#'   plot_gif(data = db)
#'
#' ## 2-dimensional example
#' db_2D <- simu_db(M = 1, covariate = TRUE)
#' hp_2D <- train_gp(db_2D)
#' pred_gif(db_2D, hp = hp_2D) %>%
#'   plot_gif(data = db_2D)
#' }
plot_gif <- function(pred_gp,
                     x_input = NULL,
                     data = NULL,
                     data_train = NULL,
                     prior_mean = NULL,
                     y_grid = NULL,
                     heatmap = F,
                     prob_CI = 0.95,
                     export_gif = FALSE,
                     path = "gif_gp.gif",
                     ...) {
  ## If 'heatmap' is TRUE, a grid of values on the y-axis is define
  if (heatmap) {
    if (is.null(y_grid)) {
      y_grid <- seq(
        min(pred_gp$Mean) - stats::qnorm((1 + prob_CI) / 2) * sqrt(max(pred_gp$Var)),
        max(pred_gp$Mean) + stats::qnorm((1 + prob_CI) / 2) * sqrt(max(pred_gp$Var)),
        length.out = 50
      )
    }
  }
  ## If 'data' is provided, format for a correct displaying in the GIF
  if (!is.null(data)) {
    data_anim <- tibble::tibble()
    for (j in 1:nrow(data)) {
      ## Extract the sample of the 'j' first data points
      data_anim <- data %>%
        dplyr::slice(1:j) %>%
        dplyr::mutate(Index = j) %>%
        dplyr::bind_rows(data_anim)
    }
  } else {
    data_anim <- NULL
  }
  ## Visualise the GP predictions and create the animation
  gg <- plot_gp(
    pred_gp = pred_gp,
    x_input = x_input,
    data = data_anim,
    data_train = data_train,
    prior_mean = prior_mean,
    y_grid = y_grid,
    heatmap = heatmap,
    prob_CI = prob_CI
  ) +
    gganimate::transition_states(.data$Index, ...)

  if (export_gif) {
    gganimate::animate(gg, renderer = gganimate::gifski_renderer(path))
  }

  gg %>% return()
}
