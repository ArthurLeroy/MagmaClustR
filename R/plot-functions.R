#' Plot smoothed curves of raw data
#'
#' Display raw data under the Magma format as smoothed curves.
#'
#' @param data A data frame or tibble with format : ID, Input, Output for
#'  single output configurations; Task_ID, Input_ID, Input, Output_ID, Output
#'  for multi-output configurations.
#' @param cluster A boolean indicating whether data should be coloured by
#'   cluster. Requires a column named 'Cluster'.
#' @param legend A boolean indicating whether the legend should be displayed.
#'
#' @return Graph of smoothed curves of raw data.
#'
#' @examples
#' TRUE
plot_db <- function(data, cluster = FALSE, legend = FALSE) {
  if(all(c("ID", "Input", "Output") %in% names(data))){
    # Single output case
    ## Convert Cluster into factors for a better display
    data$ID <- as.factor(data$ID)
    if (cluster) {
      ## Add a dummy column 'Cluster' if absent
      if (!("Cluster" %in% names(data))) {
        data$Cluster <- 1
      }
      ## Convert Cluster into factors for a better display
      data$Cluster <- as.factor(data$Cluster)

      gg <- ggplot2::ggplot(data) +
        ggplot2::geom_smooth(ggplot2::aes(
          x = .data$Input,
          y = .data$Output,
          group = .data$ID,
          color = .data$Cluster
        ),
        se = F,
        span = 0.5
        ) +
        ggplot2::geom_point(ggplot2::aes(
          x = .data$Input,
          y = .data$Output,
          group = .data$ID,
          color = .data$Cluster
        )) +
        ggplot2::theme_classic()
    } else {
      gg <- ggplot2::ggplot(data) +
        ggplot2::geom_smooth(ggplot2::aes(
          x = .data$Input,
          y = .data$Output,
          color = .data$ID
        ),
        se = F,
        span = 0.5
        ) +
        ggplot2::geom_point(ggplot2::aes(
          x = .data$Input,
          y = .data$Output,
          color = .data$ID
        )) +
        ggplot2::theme_classic()
    }
    if (!legend) {
      gg <- gg + ggplot2::guides(color = "none")
    }
  } else {
    # Multi-output case
    data$Task_ID <- as.factor(data$Task_ID)
    if (cluster) {
      ## Add a dummy column 'Cluster' if absent
      if (!("Cluster" %in% names(data))) {
        data$Cluster <- 1
      }
      ## Convert Cluster into factors for a better display
      data$Cluster <- as.factor(data$Cluster)

      gg <- ggplot2::ggplot(data) +
        ggplot2::geom_point(
          data = data,
          ggplot2::aes(x = .data$Input,
                       y = .data$Output,
                       group = .data$Task_ID,
                       color = .data$Cluster),
          size = 1.5,
          alpha = 0.6
        ) +
        ggplot2::geom_smooth(
          data = data,
          ggplot2::aes(x = .data$Input,
                       y = .data$Output,
                       group = .data$Task_ID,
                       color = .data$Cluster),
          se = F,
          linewidth = 0.3,
          alpha = 0.4,
          span = 0.5
        ) +
        ggplot2::facet_wrap(~.data$Output_ID,
                            labeller = ggplot2::as_labeller(function(x) paste("Output", x)),
                            scales = "free_y") +
        ggplot2::scale_color_brewer(palette = "Set1") +
        ggplot2::theme_classic() +
        ggplot2::theme(
          strip.background = ggplot2::element_blank(),
          strip.text = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.2))
        ) +
        ggplot2::labs(
          y = "Output Value",
          x = "Input",
          color = "Cluster"
        )
    } else {
      gg <- ggplot2::ggplot(data) +
        ggplot2::geom_smooth(ggplot2::aes(
          x = .data$Input,
          y = .data$Output,
          color = .data$Task_ID
        ),
        se = F,
        span = 0.5
        ) +
        ggplot2::geom_point(ggplot2::aes(
          x = .data$Input,
          y = .data$Output,
          color = .data$Task_ID
        )) +
        ggplot2::facet_wrap(~.data$Output_ID,
                            labeller = ggplot2::as_labeller(function(x) paste("Output", x)),
                            scales = "free_y") +
        ggplot2::theme_classic() +
        ggplot2::theme(
          strip.background = ggplot2::element_blank(),
          strip.text = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.2))
        ) +
        ggplot2::labs(
          y = "Output Value",
          x = "Input"
        )
    }
    if (!legend) {
      gg <- gg + ggplot2::guides(color = "none")
    }
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
#'    be displayed. If NULL (default) the 'Input' column is used for the x-axis.
#'    If providing a 2-dimensional vector, the corresponding columns are used
#'    for the x-axis and y-axis.
#' @param data (Optional) A tibble or data frame. Required columns: 'Input',
#'    'Output'. Additional columns for covariates can be specified. This
#'    argument corresponds to the raw data on which the prediction has been
#'    performed.
#' @param data_train (Optional) A tibble or data frame, containing the training
#'    data of the Magma model. The data set should have the same format as the
#'    \code{data} argument with an additional required column 'ID' for
#'    identifying the different individuals/tasks. If provided, those data are
#'    displayed as backward colourful points (each colour corresponding to one
#'    individual/task).
#' @param prior_mean (Optional) A tibble or a data frame, containing the 'Input'
#'    and associated 'Output' prior mean parameter of the GP prediction.
#' @param y_grid A vector, indicating the grid of values on the y-axis for which
#'    probabilities should be computed for heatmaps of 1-dimensional
#'    predictions. If NULL (default), a vector of length 50 is defined, ranging
#'    between the min and max 'Output' values contained in \code{pred_gp}.
#' @param heatmap A logical value indicating whether the GP prediction should be
#'    represented as a heatmap of probabilities for 1-dimensional inputs. If
#'    FALSE (default), the mean curve and associated Credible Interval are
#'    displayed.
#' @param samples A logical value indicating whether the GP prediction should be
#'    represented as a collection of samples drawn from the posterior. If
#'    FALSE (default), the mean curve and associated Credible Interval are
#'    displayed.
#' @param nb_samples A number, indicating the number of samples to be drawn from
#'    the predictive posterior distribution. For two-dimensional graphs, only
#'    one sample can be displayed.
#' @param plot_mean A logical value, indicating whether the mean prediction
#'    should be displayed on the graph when \code{samples = TRUE}.
#' @param alpha_samples A number, controlling transparency of the sample curves.
#' @param prob_CI A number between 0 and 1 (default is 0.95), indicating the
#'    level of the Credible Interval associated with the posterior mean curve.
#'    If this this argument is set to 1, the Credible Interval is not displayed.
#' @param size_data A number, controlling the size of the \code{data} points.
#' @param size_data_train A number, controlling the size of the
#'    \code{data_train} points.
#' @param alpha_data_train A number, between 0 and 1, controlling transparency
#'    of the \code{data_train} points.
#'
#' @return Visualisation of a Magma or GP prediction (optional: display data
#'    points, training data points and the prior mean function). For 1-D
#'    inputs, the prediction is represented as a mean curve and its associated
#'    95%  Credible Interval, as a collection of samples drawn from the
#'    posterior if \code{samples} = TRUE, or as a heatmap of probabilities if
#'    \code{heatmap} = TRUE. For 2-D inputs, the prediction is represented as a
#'    heatmap, where each couple of inputs on the x-axis and y-axis are
#'    associated with a gradient of colours for the posterior mean values,
#'    whereas the uncertainty is indicated by the transparency (the narrower is
#'    the Credible Interval, the more opaque is the associated colour, and vice
#'    versa)
#'
#' @export
#'
#' @examples
#' TRUE
plot_gp <- function(
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
) {
  if (prob_CI < 0 | prob_CI > 1) {
    stop("The 'prob_CI' argument should be a number between 0 and 1.")
  }
  ## Compute the quantile of the desired Credible Interval
  quant_ci <- stats::qnorm((1 + prob_CI) / 2)

  ## Check whether 'pred_gp' has a correct format
  if (pred_gp %>% is.data.frame()) {
    pred <- pred_gp
  } else if (
    is.list(pred_gp) &
      tryCatch(
        is.data.frame(pred_gp$pred),
        error = function(e) {
          FALSE
        }
      )
  ) {
    pred <- pred_gp$pred
    ## Check whether the hyper-posterior distribution is provided and extract
    if (
      tryCatch(
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
        is.null(prior_mean)
    ) {
      prior_mean <- pred_gp$hyperpost$mean
    }
  } else {
    stop(
      "The 'pred_gp' argument should either be a list containing the 'pred' ",
      "element or a data frame. Please read ?plot_gp(), and use ",
      "pred_gp() or pred_magma() for making predictions under a correct format."
    )
  }

  ## Remove 'ID' column if present
  if ("ID" %in% names(pred)) {
    pred <- pred %>% dplyr::select(-.data$ID)
  }

  ## Remove 'Reference' column if present
  if ("Reference" %in% names(pred)) {
    pred <- pred %>% dplyr::select(-.data$Reference)
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
    if (all(x_input %in% names(pred_gp))) {
      inputs <- pred[x_input]
    } else {
      stop("The names in the 'x_input' argument don't exist in 'pred_gp'.")
    }
  }

  ## Format the tibble for displaying the Credible Intervals
  pred <- pred %>%
    dplyr::mutate("CI_inf" = .data$Mean - quant_ci * sqrt(.data$Var)) %>%
    dplyr::mutate("CI_sup" = .data$Mean + quant_ci * sqrt(.data$Var)) %>%
    dplyr::mutate("CI_Width" = .data$CI_sup - .data$CI_inf)

  ## If no CI (i.e. prob_CI = 1), then no transparency (i.e. alpha = 1)
  if (prob_CI == 0) {
    pred$CI_width <- 1
  }

  ## Display a heatmap if inputs are 2D
  if (ncol(inputs) == 2) {
    if (samples) {
      ## Display samples from the posterior
      gg <- plot_samples(
        pred = pred_gp,
        x_input = x_input,
        nb_samples = nb_samples,
        plot_mean = plot_mean,
        alpha_samples = alpha_samples
      )
    } else {
      ## Add the 'Index' column if the prediction comes from 'pred_gif()'
      if (!is.null(index)) {
        pred <- pred %>% dplyr::mutate("Index" = index)
      }

      gg <- ggplot2::ggplot() +
        ggplot2::geom_raster(
          data = pred,
          ggplot2::aes(
            x = .data[[names(inputs)[1]]],
            y = .data[[names(inputs)[2]]],
            fill = .data$Mean,
            alpha = .data$CI_Width
          ),
          interpolate = TRUE
        ) +
        ggplot2::scale_fill_distiller(palette = "RdPu", trans = "reverse") +
        ggplot2::scale_alpha_continuous(range = c(0.1, 1), trans = "reverse")
    }

    if (!is.null(data)) {
      ## Round the 'Output' values to reduce size of labels on the graph
      data <- data %>% dplyr::mutate(Output = round(.data$Output, 1))

      gg <- gg +
        ggplot2::geom_label(
          data = data,
          ggplot2::aes(
            x = .data[[names(inputs)[1]]],
            y = .data[[names(inputs)[2]]],
            label = .data$Output,
            fill = .data$Output
          ),
          size = 3
        )
    }

    ## If some day I want to add the feature of displaying data_train in 2D

    # if (!is.null(data_train)) {
    #   ## Round the 'Output' values to reduce size of labels on the graph
    #   data_train <- data_train %>% dplyr::mutate(Output = round(.data$Output, 1))
    #
    #   gg <- gg + ggplot2::geom_label(
    #     data = data_train,
    #     ggplot2::aes(
    #       x = .data[[names(inputs)[1]]],
    #       y = .data[[names(inputs)[2]]],
    #       label = .data$Output,
    #       colour = .data$ID
    #     ),
    #     size = 3
    #   ) + ggplot2::guides(colour = 'none')
    # }
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
        tidyr::expand(
          tidyr::nesting(!!!rlang::syms(col_to_nest)),
          "Ygrid" = y_grid
        ) %>%
        dplyr::mutate(
          "Proba" = 2 *
            (1 -
              stats::pnorm(abs((.data$Ygrid - .data$Mean) / sqrt(.data$Var))))
        )

      gg <- ggplot2::ggplot() +
        ggplot2::geom_raster(
          data = db_heat,
          ggplot2::aes(
            x = .data[[names(inputs)[1]]],
            y = .data$Ygrid,
            fill = .data$Proba
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
    } else if (samples) {
      ## Display samples from the posterior
      gg <- plot_samples(
        pred = pred_gp,
        x_input = x_input,
        nb_samples = nb_samples,
        plot_mean = plot_mean,
        alpha_samples = alpha_samples
      )
    } else {
      ## Display a classic curve otherwise
      ## Add the 'Index' column if the prediction comes from 'pred_gif()'
      if (!is.null(index)) {
        pred <- pred %>% dplyr::mutate("Index" = index)
      }

      gg <- ggplot2::ggplot() +
        ggplot2::geom_line(
          data = pred,
          ggplot2::aes(x = .data[[names(inputs)[1]]], y = .data$Mean),
          color = "#DB15C1"
        ) +
        ggplot2::geom_ribbon(
          data = pred,
          ggplot2::aes(
            x = .data[[names(inputs)[1]]],
            ymin = .data$CI_inf,
            ymax = .data$CI_sup
          ),
          alpha = 0.2,
          fill = "#FA9FB5"
        ) +
        ggplot2::ylab("Output")
    }

    ## Display the training data if provided
    if (!is.null(data_train)) {
      gg <- gg +
        ggplot2::geom_point(
          data = data_train,
          ggplot2::aes(
            x = .data[[names(inputs)[1]]],
            y = .data$Output,
            col = .data$ID
          ),
          size = size_data_train,
          alpha = alpha_data_train
        ) +
        ggplot2::guides(color = "none")
    }
    ## Display the observed data if provided
    if (!is.null(data)) {
      gg <- gg +
        ggplot2::geom_point(
          data = data,
          ggplot2::aes(x = .data[[names(inputs)[1]]], y = .data$Output),
          size = size_data,
          shape = 20
        )
    }

    ## Display the (hyper-)prior mean process if provided
    if (!is.null(prior_mean)) {
      if (names(inputs)[1] %in% names(prior_mean)) {
        gg <- gg +
          ggplot2::geom_line(
            data = prior_mean,
            ggplot2::aes(x = .data[[names(inputs)[1]]], y = .data$Output),
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

#' @rdname plot_gp
#' @export
plot_magma <- plot_gp

#' Display realisations from a (mixture of) GP prediction
#'
#' Display samples drawn from the posterior of a GP, Magma or
#' MagmaClust prediction. According to the dimension of the inputs, the graph
#' may represent curves or a heatmap.
#'
#' @param pred A list, typically coming from \code{\link{pred_gp}},
#'    \code{\link{pred_magma}} or \code{\link{pred_magmaclust}} functions, using
#'    the argument 'get_full_cov = TRUE'. Required elements: \code{pred},
#'    \code{cov}. This argument is needed if \code{samples} is missing.
#' @param samples A tibble or data frame, containing the samples generated from
#'    a GP, Magma, or MagmaClust prediction. Required columns: \code{Input},
#'    \code{Sample}, \code{Output}.  This argument is needed if \code{pred}
#'    is missing.
#' @param nb_samples A number, indicating the number of samples to be drawn from
#'    the predictive posterior distribution. For two-dimensional graphs, only
#'    one sample can be displayed.
#' @param x_input A vector of character strings, indicating which 'column'
#'    should be displayed in the case of multidimensional inputs. If
#'    NULL(default) the Input' column is used for the x-axis. If providing a
#'    2-dimensional vector, the corresponding columns are used for the x-axis
#'    and the y-axis.
#' @param plot_mean A logical value, indicating whether the mean prediction
#'    should be displayed on the graph.
#' @param alpha_samples A number, controlling transparency of the sample curves.
#'
#' @return Graph of samples drawn from a posterior distribution of a GP,
#'    Magma, or MagmaClust prediction.
#' @export
#'
#' @examples
#' TRUE
plot_samples <- function(
  pred = NULL,
  samples = NULL,
  nb_samples = 50,
  x_input = NULL,
  plot_mean = TRUE,
  alpha_samples = 0.3
) {
  ## Check whether 'samples' or 'pred' exist
  if (is.null(samples) & is.null(pred)) {
    stop("Either 'sample' or 'pred' is needed as an argument.")
  }

  ## If provided, check format of 'pred' and extract the mixture prediction
  if (!is.null(pred)) {
    ## Check 'pred' format
    if (!(is.list(pred) & ('cov' %in% names(pred)))) {
      stop(
        "The 'pred' argument should be a list containing 'pred' and 'cov' ",
        "elements. Consider re-running the prediction function using the ",
        "argument 'get_full_cov' = TRUE."
      )
    }

    ## Check whether 'pred' is a GP/Magma or a MagmaClust prediction
    if ('mixture_pred' %in% names(pred)) {
      ## If 'samples' is not provided, draw new samples
      if (is.null(samples)) {
        samples = sample_magmaclust(pred_clust = pred, nb_samples = nb_samples)
      }

      mean_pred = pred$mixture_pred
    } else {
      ## If 'samples' is not provided, draw new samples
      if (is.null(samples)) {
        samples = sample_gp(pred_gp = pred, nb_samples = nb_samples)
      }

      mean_pred = pred$pred
    }
  }

  ## Get the inputs that should be used
  if (x_input %>% is.null()) {
    inputs <- samples %>% dplyr::select(-c(.data$Sample, .data$Output))
  } else {
    inputs <- samples[x_input]
  }

  ## Display a heatmap if inputs are 2D
  if (ncol(inputs) == 2) {
    ## Extract only one sample when displaying in 2D
    samples <- samples %>%
      dplyr::filter(.data$Sample == unique(samples$Sample)[1]) %>%
      dplyr::select(-.data$Sample)

    gg <- ggplot2::ggplot() +
      ggplot2::geom_raster(
        data = samples,
        ggplot2::aes(
          x = .data[[names(inputs)[1]]],
          y = .data[[names(inputs)[2]]],
          fill = .data$Output
        ),
        interpolate = TRUE
      ) +
      ggplot2::scale_fill_distiller(palette = "RdPu", trans = "reverse")
  } else if (ncol(inputs) == 1) {
    gg <- ggplot2::ggplot() +
      ggplot2::geom_line(
        data = samples,
        ggplot2::aes(
          x = .data[[names(inputs)[1]]],
          y = .data$Output,
          group = .data$Sample
        ),
        color = "#FA9FB5",
        alpha = alpha_samples
      ) +
      ggplot2::guides(group = "none")

    if (plot_mean) {
      if (is.null(pred)) {
        warning("The 'pred' argument is needed to display the mean prediction.")
      } else {
        gg = gg +
          ggplot2::geom_line(
            data = mean_pred,
            ggplot2::aes(x = .data[[names(inputs)[1]]], y = .data$Mean),
            color = "#DB15C1"
          )
      }
    }
  } else {
    stop(
      "Impossible to display inputs with dimensions greater than 2. Please ",
      "provide two elements or less in the 'x_axis' argument."
    )
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
#' @param size_data A number, controlling the size of the \code{data} points.
#' @param size_data_train A number, controlling the size of the
#'    \code{data_train} points.
#' @param alpha_data_train A number, between 0 and 1, controlling transparency
#'    of the \code{data_train} points.
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
#' TRUE
plot_gif <- function(
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
) {
  ## If 'heatmap' is TRUE, a grid of values on the y-axis is define
  if (heatmap) {
    if (is.null(y_grid)) {
      y_grid <- seq(
        min(pred_gp$Mean) -
          stats::qnorm((1 + prob_CI) / 2) * sqrt(max(pred_gp$Var)),
        max(pred_gp$Mean) +
          stats::qnorm((1 + prob_CI) / 2) * sqrt(max(pred_gp$Var)),
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
    prob_CI = prob_CI,
    size_data = size_data,
    size_data_train = size_data_train,
    alpha_data_train = alpha_data_train
  ) +
    gganimate::transition_states(.data$Index, ...)

  if (export_gif) {
    gganimate::animate(gg, renderer = gganimate::gifski_renderer(path))
  }

  gg %>%
    invisible() %>%
    return()
}

#' Plot MagmaClust predictions
#'
#' Display MagmaClust predictions. According to the dimension of the
#' inputs, the graph may be a mean curve (dim inputs = 1) or a heatmap
#' (dim inputs = 2) of probabilities. Moreover, MagmaClust can provide credible
#' intervals only by visualising cluster-specific predictions (e.g. for the most
#'  probable cluster). When visualising the full mixture-of-GPs prediction,
#'  which can be multimodal, the user should choose between the simple mean
#'  function or the full heatmap of probabilities (more informative but slower).
#'
#' @param pred_clust A list of predictions, typically coming from
#'    \code{\link{pred_magmaclust}}. Required elements: \code{pred},
#'    \code{mixture}, \code{mixture_pred}.
#' @param cluster A character string, indicating which cluster to plot from.
#'    If 'all' (default) the mixture of GPs prediction is displayed as a mean
#'    curve (1-D inputs) or a mean heatmap (2-D inputs). Alternatively, if the
#'    name of one cluster is provided, the classic mean curve + credible
#'    interval is displayed (1-D inputs), or a heatmap with colour gradient for
#'    the mean and transparency gradient for the Credible Interval (2-D inputs).
#' @param x_input A vector of character strings, indicating which input should
#'    be displayed. If NULL (default) the 'Input' column is used for the x-axis.
#'    If providing a 2-dimensional vector, the corresponding columns are used
#'    for the x-axis and y-axis.
#' @param data (Optional) A tibble or data frame. Required columns: \code{Input}
#'    , \code{Output}. Additional columns for covariates can be specified. This
#'    argument corresponds to the raw data on which the prediction has been
#'    performed.
#' @param data_train (Optional) A tibble or data frame, containing the training
#'    data of the MagmaClust model. The data set should have the same format as
#'    the \code{data} argument with an additional required column \code{ID} for
#'    identifying the different individuals/tasks. If provided, those data are
#'    displayed as backward colourful points (each colour corresponding to one
#'    individual or a cluster, see \code{col_clust} below).
#' @param data_test (Optional) A tibble or data frame. Required columns:
#'    \code{Input}, \code{Output}. These are data points withheld from the
#'    prediction and overlaid on the plot as a distinct layer, allowing visual
#'    evaluation of prediction quality at unobserved locations.
#' @param col_clust A boolean indicating whether backward points are coloured
#'    according to the individuals or to their most probable cluster. If one
#'    wants to colour by clusters, a column \code{Cluster} shall be present in
#'    \code{data_train}. We advise to use \code{\link{data_allocate_cluster}}
#'    for automatically creating a well-formatted dataset from a trained
#'    MagmaClust model.
#' @param prior_mean (Optional) A list providing, for each cluster, a
#'    tibble containing prior mean parameters of the prediction. This argument
#'    typically comes as an outcome \code{hyperpost$mean}, available through
#'    the \code{\link{train_magmaclust}}, \code{\link{pred_magmaclust}}
#'    functions.
#' @param y_grid A vector, indicating the grid of values on the y-axis for which
#'    probabilities should be computed for heatmaps of 1-dimensional
#'    predictions. If NULL (default), a vector of length 50 is defined, ranging
#'    between the min and max 'Output' values contained in \code{pred}.
#' @param heatmap A logical value indicating whether the GP mixture should be
#'    represented as a heatmap of probabilities for 1-dimensional inputs. If
#'    FALSE (default), the mean curve (and associated Credible Interval if
#'    available) are displayed.
#' @param samples A logical value indicating whether the GP mixture should be
#'    represented as a collection of samples drawn from the posterior. If
#'    FALSE (default), the mean curve (and associated Credible Interval if
#'    available) are displayed.
#' @param nb_samples A number, indicating the number of samples to be drawn from
#'    the predictive posterior distribution. For two-dimensional graphs, only
#'    one sample can be displayed.
#' @param plot_mean A logical value, indicating whether the mean prediction
#'    should be displayed on the graph when \code{samples = TRUE}.
#' @param alpha_samples A number, controlling transparency of the sample curves.
#' @param prob_CI A number between 0 and 1 (default is 0.95), indicating the
#'    level of the Credible Interval associated with the posterior mean curve.
#'    If this this argument is set to 1, the Credible Interval is not displayed.
#' @param size_data A number, controlling the size of the \code{data} points.
#' @param size_data_train A number, controlling the size of the
#'    \code{data_train} points.
#' @param alpha_data_train A number, between 0 and 1, controlling transparency
#'    of the \code{data_train} points.
#'
#' @return Visualisation of a MagmaClust prediction (optional: display data
#'    points, training data points and the prior mean functions). For 1-D
#'    inputs, the prediction is represented as a mean curve (and its associated
#'    95% Credible Interval for cluster-specific predictions), or as a heatmap
#'    of probabilities if \code{heatmap} = TRUE. In the case of MagmaClust,
#'    the heatmap representation should be preferred for clarity, although the
#'    default display remains mean curve for quicker execution. For 2-D inputs,
#'    the prediction is represented as a heatmap, where each couple of inputs on
#'    the x-axis and y-axis are associated with a gradient of colours for the
#'    posterior mean values, whereas the uncertainty is indicated by the
#'    transparency (the narrower is the Credible Interval, the more opaque is
#'    the associated colour, and vice versa). As for 1-D inputs, Credible
#'    Interval information is only available for cluster-specific predictions.
#'
#' @export
#'
#' @examples
#' TRUE
plot_magmaclust <- function(
  pred_clust,
  cluster = "all",
  x_input = NULL,
  data = NULL,
  data_train = NULL,
  data_test = NULL,
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
) {
  ## Check prob_CI format
  if (prob_CI < 0 | prob_CI > 1) {
    stop("The 'prob_CI' argument should be a number between 0 and 1.")
  }

  ## Check format for prediction
  if (!is.list(pred_clust)) {
    stop("Wrong format for 'pred_clust', please read ?plot_magmaclust().")
  }
  ## Check presence of 'pred'
  if (!("pred" %in% names(pred_clust))) {
    stop("Wrong format for 'pred_clust', please read ?plot_magmaclust().")
  }
  ## Check presence of 'mixture'
  if (!("mixture" %in% names(pred_clust))) {
    stop("Wrong format for 'pred_clust', please read ?plot_magmaclust().")
  }
  ## Check presence of 'mixture_pred'
  if (!("mixture_pred" %in% names(pred_clust))) {
    stop("Wrong format for 'pred_clust', please read ?plot_magmaclust().")
  }

  pred <- pred_clust$pred
  mixture <- pred_clust$mixture
  mixture_pred <- pred_clust$mixture_pred
  ID_k <- names(pred)

  ## Checker whether we can provide Credible Interval
  if (cluster == "all") {
    ## Check whether one cluster's proba is 1 (or really close)
    max_clust <- proba_max_cluster(mixture)
    if (round(max_clust$Proba[1], 3) == 1) {
      cluster <- max_clust$Cluster
      cat(
        "The mixture probability of the cluster",
        cluster,
        "is 1. Therefore,",
        "the predictive distribution is Gaussian and the associated",
        "credible interval can be displayed. \n\n"
      )
      ## Create dummy variable for indicating the type of prediction
      all_clust <- FALSE

      ## Compute the quantile of the desired Credible Interval
      quant_ci <- stats::qnorm((1 + prob_CI) / 2)
    } else {
      all_clust <- TRUE
    }
  } else {
    all_clust <- FALSE
    ## Compute the quantile of the desired Credible Interval
    quant_ci <- stats::qnorm((1 + prob_CI) / 2)
  }

  ## Select the appropriate tibble for displaying predictions
  if (all_clust) {
    pred_gp <- mixture_pred
  } else {
    ## Check the name provided in 'cluster'
    if (!(cluster %in% ID_k)) {
      stop(
        "The cluster's name provided in 'cluster' does not exist in) ",
        "'pred_clust'."
      )
    }
    ## Remove the 'Proba' column if selecting cluster-specific prediction
    pred_gp <- pred[[cluster]] %>% dplyr::select(-.data$Proba)

    ## Get the 'Proba' value to display in the Title
    proba <- pred[[cluster]] %>%
      dplyr::pull(.data$Proba) %>%
      unique()
  }

  ## Get the inputs that should be used
  if (x_input %>% is.null()) {
    inputs <- pred_gp %>% dplyr::select(-c(.data$ID, .data$Mean, .data$Var))
  } else {
    if (all(x_input %in% names(pred_gp))) {
      inputs <- pred_gp[x_input]
    } else {
      stop("The names in the 'x_input' argument don't exist in 'pred_clust'.")
    }
  }

  if (samples) {
    ## Display samples drawn from the posterior mixture of GPs
    gg <- plot_samples(
      pred = pred_clust,
      nb_samples = nb_samples,
      x_input = x_input,
      plot_mean = plot_mean,
      alpha_samples = alpha_samples
    )

    ## Display the observed data if provided
    if (!is.null(data)) {
      ## Check dimension of the inputs
      if (ncol(inputs) == 2) {
        ## Display labels if 2-D
        gg <- gg +
          ggplot2::geom_label(
            data = data,
            ggplot2::aes(
              x = .data[[names(inputs)[1]]],
              y = .data[[names(inputs)[2]]],
              label = .data$Output,
              fill = .data$Output
            ),
            size = 3
          )
      } else if (ncol(inputs) == 1) {
        ## Display points if 1-D
        gg <- gg +
          ggplot2::geom_point(
            data = data,
            ggplot2::aes(x = .data[[names(inputs)[1]]], y = .data$Output),
            size = size_data,
            shape = 20
          )
      }
    }
    ## Define the adequate title
    gtitle <- paste0("Samples from a mixture of GPs")
  } else if (all_clust) {
    ## GP visualisation without Credible Interval
    gg <- plot_gp(
      pred_gp = pred_gp,
      x_input = x_input,
      data = data,
      y_grid = y_grid,
      heatmap = heatmap,
      prob_CI = 0,
      size_data = size_data
    )

    ## Define the adequate title
    gtitle <- paste0("Mixture of GP predictions")
  } else {
    ## Classic GP visualisation for cluster-specific predictions
    gg <- plot_gp(
      pred_gp = pred_gp,
      x_input = x_input,
      data = data,
      y_grid = y_grid,
      heatmap = heatmap,
      prob_CI = prob_CI,
      size_data = size_data
    )

    ## Define the adequate title
    gtitle <- paste0("Cluster ", cluster, " -- Proba = ", mixture[[cluster]])
  }

  ## Display training data if available
  if (!is.null(data_train)) {
    ## Change colours of background points depending on 'col_clust'
    if (col_clust) {
      ## Check whether 'data_train' provides a 'Cluster' column
      if (!("Cluster" %in% names(data_train))) {
        cat(
          "The 'data_train' argument does not provide a 'Cluster' column.",
          "Therefore, training data remain coloured by individual. \n \n"
        )
      } else {
        ## Certify that 'Cluster' is discrete
        data_train$Cluster <- as.factor(data_train$Cluster)

        ## Colour training data plot by cluster
        gg <- gg +
          ggplot2::geom_point(
            data = data_train,
            ggplot2::aes(
              x = .data[[names(inputs)[1]]],
              y = .data$Output,
              col = .data$Cluster
            ),
            size = size_data_train,
            alpha = alpha_data_train
          )
      }
    } else {
      gg <- gg +
        ggplot2::geom_point(
          data = data_train,
          ggplot2::aes(
            x = .data[[names(inputs)[1]]],
            y = .data$Output,
            fill = .data$ID
          ),
          shape = 21,
          size = size_data_train,
          alpha = alpha_data_train
        ) +
        ggplot2::guides(fill = "none")
    }
  }

  ## Add testing point if provided

  if (!is.null(data_test)) {
    gg <- gg +
      ggplot2::geom_point(
        data = data_test,
        ggplot2::aes(x = .data[[names(inputs)[1]]], y = .data$Output),
        colour = "red",
        size = size_data,
        shape = 18
      )
  }

  ## Display the prior mean process if provided
  if (!is.null(prior_mean)) {
    ## Extract 'mean' if the user provides the full 'hyperpost'
    ## Bind the tibbles of hyper-posterior mean processes
    mean_k <- prior_mean %>% dplyr::bind_rows(.id = "Cluster")

    ## Display the mean functions if available
    if (names(inputs)[1] %in% names(mean_k)) {
      gg <- gg +
        ggplot2::geom_line(
          data = mean_k,
          ggplot2::aes(
            x = .data[[names(inputs)[1]]],
            y = .data$Output,
            col = .data$Cluster
          ),
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

  ## Change scale colour palette
  gg <- gg + ggplot2::scale_color_brewer(palette = "Set1")

  (gg + ggplot2::ggtitle(gtitle)) %>%
    return()
}

#' Plot time series grouped by cluster
#'
#' After training a MagmaClust model, this function assigns each individual
#' to its most likely cluster and plots all time series in that cluster
#' together — one panel per cluster, arranged in a grid.
#'
#' @param trained_model The object returned by \code{\link{train_magmaclust}}.
#' @param panels A boolean, indicating whether a single plot with all clusters
#'    should be displayed (FALSE), or a panel of plots for each cluster (TRUE).
#' @param grid_inputs The grid of input (reference Input and covariates) values
#'    on which the hyperposterior should be re-evaluated. If \code{NULL}
#'    (default), the hyperposterior evaluated on the data grid (from training)
#'    is returned. Ideally, this argument should be a
#'    tibble or a data frame, providing the same columns as \code{data}, except
#'    'Output'. Nonetheless, in cases where \code{data} provides only one
#'    'Input' column, the \code{grid_inputs} argument can be a vector.
#' @param remove_empty A boolean, indicating whether empty clusters should be
#'    removed from the plot.
#' @param displayed_clusters A vector of characters, indicating which subset of
#'    clusters should be displayed. When \code{NULL}, they are all plotted.
#' @param ncol A number, indicating the number of columns in the panel grid.
#'    If NULL (default), the value is chosen automatically: 1 for a single
#'    cluster, 2 for 2–4 clusters, and 3 for 5 or more.
#' @param prob_CI A number between 0 and 1 (default is 0.95), indicating the
#'    level of the Credible Interval associated with the posterior mean curve.
#'    If this this argument is set to 1, the Credible Interval is not displayed.
#' @param x_input A character string, specifying the Input used in the x-axis.
#'     As this visual representation is only available in 1-dimension, when
#'     multiple Inputs exist, the marginal according to the chosen Input will be
#'     displayed. Default is \code{"Input"}.
#' @param y_lab A character string, specifying the y-axis label. Default is
#'    \code{"Output"}.
#' @param alpha_point A number, between 0 and 1, controlling transparency of
#'    the individual data points.
#' @param size_point A number, controlling the size of the data points.
#' @param show_prob A logical value, indicating whether the mean assignment
#'    probability should be appended to each panel title (e.g.
#'    \emph{K1 (n = 5, mean prob = 0.91)}). Default is \code{TRUE}.
#' @param show_legend A logical value, indicating whether a colour legend
#'    mapping line colours to individual IDs should be displayed. Set to
#'    \code{FALSE} for large clusters where the legend would be unreadable.
#'    Default is \code{TRUE}.
#'
#' @return A visual representation of clusters resulting from MagmaClust
#'         training. If \code{panels = FALSE}, a single ggplot object is
#'         returned, and a (cowplot) panel for each cluster
#'         is displayed when \code{panels = TRUE}.
#'
#' @seealso \code{\link{train_magmaclust}}, \code{\link{plot_magmaclust}},
#'    \code{\link{data_allocate_cluster}}
#'
#' @export
#'
#' @examples
#' TRUE
plot_clusters <- function(
    trained_model,
    panels = FALSE,
    grid_inputs = NULL,
    remove_empty = TRUE,
    displayed_clusters = NULL,
    x_input = "Input",
    y_lab = "Output",
    ncol = NULL,
    prob_CI = 0.95,
    alpha_point = 0.5,
    size_point = 1.5,
    show_prob = TRUE,
    show_legend = FALSE) {

  ## Extract cluster assignments from the trained model
  mixture <- trained_model$hyperpost$mixture

  if(grid_inputs %>% is.null()){
    hyperpost = trained_model$hyperpost
  } else {
    ## Recompute hyperposterior mean processes on a fine grid
    hyperpost = quiet(hyperposterior_clust(trained_model = trained_model,
                                           grid_inputs = grid_inputs))
  }


  ## Extract the dataset used to train the model
  data = trained_model$ini_args$data

  ## Get the names of clusters
  k_cols <- names(mixture)

  ## Assign task to the most probable cluster
  data_clustered <- trained_model %>%
    data_allocate_cluster()

  assignments = data_clustered %>%
    dplyr::select(.data$ID, .data$Cluster, .data$Proba)


  ## Extract the list of all clusters' names
  all_clust = trained_model$hp_k$ID %>% unique()

  ## If users want to display a subset of clusters
  if(displayed_clusters %>% is.null()){

    ## Remove empty cluster
    if(remove_empty){

      ## Identify non-empty clusters and filter to requested subset
      non_empty <- assignments %>%
        dplyr::distinct(.data$Cluster) %>%
        dplyr::pull(.data$Cluster)

      ## Identify empty clusters
      empty = all_clust[ which(!(all_clust %in% non_empty) )]

      if( !(length(empty) == 0) )
      {
        cat("Clusters", empty, "are empty, and were removed from the plot.")
      }

      clusters_to_plot <- non_empty

    } else{
      clusters_to_plot <- all_clust
    }

  } else {
    clusters_to_plot <- displayed_clusters
  }

  if(panels == FALSE){
    mean_k <- hyperpost$mean %>%
      dplyr::bind_rows(.id = "Cluster") %>%
      dplyr::filter(.data$Cluster %in% clusters_to_plot)

    data_k = data_clustered %>%
      dplyr::filter(.data$Cluster %in% clusters_to_plot)

    ## Display the mean functions if available
    gg = ggplot2::ggplot() +
      ggplot2::geom_line(
        data = mean_k,
        ggplot2::aes(
          x = .data[[x_input]],
          y = .data$Output,
          col = .data$Cluster
        ),
        linetype = "dashed"
      ) +
      ggplot2::geom_point(data = data_k, ggplot2::aes(
        x = .data[[x_input]],
        y = .data$Output,
        col = .data$Cluster), size = size_point, alpha = alpha_point) +
      ggplot2::theme_classic()

    return(gg)
  }

  ## Build one ggplot panel per cluster
  loop_panel = function(k) {

    quant_ci <- stats::qnorm((1 + prob_CI) / 2)

    hyperpost_k = hyperpost$pred[[k]] %>%
      dplyr::mutate("CI_inf" = .data$Mean - quant_ci * sqrt(.data$Var)) %>%
      dplyr::mutate("CI_sup" = .data$Mean + quant_ci * sqrt(.data$Var) )

    k_data <- data_clustered %>%
      dplyr::filter(.data$Cluster == k)

    n_ids <- dplyr::n_distinct(k_data$ID)

    mean_prob <- assignments %>%
      dplyr::filter(.data$Cluster == k) %>%
      dplyr::pull(.data$Proba) %>%
      mean() %>%
      round(2)

    ## Define panel title with or without assignment probability
    title <- if(show_prob){
      paste0(k, "  (n = ", n_ids, ", mean probability = ", mean_prob, ")")
    } else {
      paste0(k, "  (n = ", n_ids, ")")
    }

    gg <- ggplot2::ggplot() +
      ## Plot data points
      ggplot2::geom_point(data = k_data, ggplot2::aes(
        x = .data[[x_input]],
        y = .data$Output,
        col = .data$ID), size = size_point, alpha = alpha_point) +
      ggplot2::labs(title = title, x = x_input, y = y_lab) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        legend.position = "right",
        legend.justification = "top") +
      ## Add mean process
      ggplot2::geom_line(
        data = hyperpost_k,
        ggplot2::aes(
          x = .data[[x_input]],
          y = .data$Mean),
        linetype = "dashed") +
      ggplot2::geom_ribbon(
        data = hyperpost_k,
        ggplot2::aes(
          x = .data[[x_input]],
          ymin = .data$CI_inf,
          ymax = .data$CI_sup
        ),
        alpha = 0.2,
        fill = "#FA9FB5"
      )


    ## Remove legend if requested
    if (!show_legend){gg <- gg + ggplot2::theme(legend.position = "none")}

    return(gg)
  }

  panels <- lapply(clusters_to_plot,loop_panel)

  ## Return a single panel directly, or arrange multiple panels in a grid
  if (length(panels) == 1) return(panels[[1]])

  cowplot::plot_grid(plotlist = panels, ncol = ncol, align = "hv") %>%
    return()
}
