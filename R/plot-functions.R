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
plot_db <- function(data, cluster = FALSE, legend = FALSE) {
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
    gg <- ggplot2::ggplot(data) +
      ggplot2::geom_smooth(ggplot2::aes(
        x = .data$Input,
        y = .data$Output,
        color = .data$ID
      ),
      se = F
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
#'    'Input_ID', 'Output', 'Output_ID'. This argument corresponds to the raw
#'    data on which the prediction has been performed.
#' @param data_train (Optional) A tibble or data frame, containing the training
#'    data of the Magma model. The data set should have the same format as the
#'    \code{data} argument with an additional required column 'Task_ID' for
#'    identifying the different individuals/tasks. If provided, those data are
#'    displayed as backward colorful points (each color corresponding to one
#'    individual/task).
#' @param prior_mean (Optional) A tibble or a data frame, containing the 'Input'
#'    and associated 'Output'and 'Output_ID' prior mean parameter of the GP prediction.
#' @param y_grid A vector, indicating the grid of values on the y-axis for which
#'    probabilities should be computed for heatmaps of 1-dimensional
#'    predictions. If NULL (default), a vector of length 50 is defined, ranging
#'    between the min and max 'Output' values for each Output_ID contained in
#'    \code{pred_gp}.
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
#' @return A list of visualisation of a Magma (single or multi-output) or GP
#'    prediction (single or multi-output) (optional: display data points,
#'    training data points and the prior mean function). Each element of the
#'    list corresponds to a prediction associated to one Output_ID; therefore,
#'    the list is as long as there are different Output_ID in the data.
#'    For 1-D inputs, the prediction is represented as a mean curve and its
#'    associated 95%  Credible Interval, as a collection of samples drawn from
#'    the posterior if \code{samples} = TRUE, or as a heatmap of probabilities
#'    if \code{heatmap} = TRUE. For 2-D inputs, the prediction is represented as
#'    a heatmap, where each couple of inputs on the x-axis and y-axis are
#'    associated with a gradient of colours for the posterior mean values,
#'    whereas the uncertainty is indicated by the transparency (the narrower is
#'    the Credible Interval, the more opaque is the associated colour, and vice
#'    versa).
#'
#' @export
#'
#' @examples
#' TRUE
plot_gp <- function(pred_gp,
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
                    alpha_data_train = 0.5) {

  if (prob_CI < 0 | prob_CI > 1) {
    stop("The 'prob_CI' argument should be a number between 0 and 1.")
  }
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
      "element or a data frame."
    )
  }

  ## Remove 'ID' column if present
  if ("Task_ID" %in% names(pred)) {
    pred <- pred %>% dplyr::select(-Task_ID)
  }

  ## Remove 'Reference' column if present
  if ("Reference" %in% names(pred)) {
    pred <- pred %>% dplyr::select(-Reference)
  }

  ## Remove the 'Index' column if the prediction comes from 'pred_gif()'
  if (any("Index" %in% names(pred))) {
    index <- pred %>% dplyr::pull(Index)
    pred <- pred %>% dplyr::select(-Index)
  } else {
    index <- NULL
  }

  ## Rename the 'Output' column for enabling plot of the mean process in Magma
  if ("Output" %in% names(pred)) {
    pred <- pred %>% dplyr::rename("Mean" = Output)
  }

  all_samples <- NULL
  if (samples) {
    if (!is.list(pred_gp) || is.null(pred_gp$cov)) {
      stop("If 'samples = TRUE', 'pred_gp' should be an object from ",
           "pred_gp() with 'get_full_cov = TRUE'.")
    }

    # Call sample_gp()
    all_samples <- tryCatch(
      sample_gp(pred_gp = pred_gp, nb_samples = nb_samples),
      error = function(e) {
        stop(paste("Error when calling sample_gp():", e$message))
      }
    )
    if (!"Output_ID" %in% names(all_samples)) {
      stop("sample_gp() does not contain an 'Output_ID' column.")
    }
  }

  if (!"Output_ID" %in% names(pred)) {
    stop("'Output_ID' not found in 'pred'.")
  }

  unique_outputs <- as.character(unique(pred$Output_ID))

  ## Loop on each output
  plot_list <- purrr::map(unique_outputs, function(current_output_id) {
    ## Subset pred only on the current Output_ID
    pred_subset <- pred %>% dplyr::filter(Output_ID == current_output_id)

    # Subset data only on the current Output_ID
    data_subset <- NULL
    if (!is.null(data)) {
      data_subset <- if ("Output_ID" %in% names(data)) {
        data %>% dplyr::filter(Output_ID == current_output_id)
      } else {
        data
      }
    }

    # Subset data_train only on the current Output_ID
    data_train_subset <- NULL
    if (!is.null(data_train)) {
      data_train_subset <- if ("Output_ID" %in% names(data_train)) {
        data_train %>% dplyr::filter(Output_ID == current_output_id)
      } else {
        data_train
      }
    }

    # Subset prior_mean only on the current Output_ID
    prior_mean_subset <- NULL
    if (!is.null(prior_mean)) {
      prior_mean_subset <- if ("Output_ID" %in% names(prior_mean)) {
        prior_mean %>% dplyr::filter(Output_ID == current_output_id)
      } else {
        prior_mean
      }
    }

    # Subset samples only on the current Output_ID
    samples_subset <- NULL
    if (samples && !is.null(all_samples)) {
      samples_subset <- all_samples %>% dplyr::filter(Output_ID == current_output_id)
    }

    # Extract inputs corresponding to the current Output_ID
    if (x_input %>% is.null()) {
      inputs <- pred_subset %>% dplyr::select(-c(Mean, Var, Output_ID))
    } else {
      if (all(x_input %in% names(pred_subset))) {
        inputs <- pred_subset[x_input]
      } else {
        stop("The names in the 'x_input' argument don't exist in 'pred_gp'.")
      }
    }

    # Compute IC for the current Output_ID
    pred_subset <- pred_subset %>%
      dplyr::mutate("CI_inf" = Mean - quant_ci * sqrt(Var)) %>%
      dplyr::mutate("CI_sup" = Mean + quant_ci * sqrt(Var)) %>%
      dplyr::mutate("CI_Width" = CI_sup - CI_inf)

    if (prob_CI == 0) {
      pred_subset$CI_width <- 1
    }

    ## Display a heatmap if inputs are 2D
    if (ncol(inputs) == 2) {

      if (samples){
        samples_2d <- samples_subset %>%
          dplyr::filter(.data$Sample == unique(samples_subset$Sample)[1]) %>%
          dplyr::select(- .data$Sample)

        gg <- ggplot2::ggplot() +
          ggplot2::geom_raster(
            data = samples_2d,
            ggplot2::aes_string(
              x = names(inputs)[1],
              y = names(inputs)[2],
              fill = "Output"
            ),
            interpolate = TRUE
          ) +
          ggplot2::scale_fill_distiller(palette = "RdPu", trans = "reverse")

      } else {
        if (!is.null(index)) {
          pred_subset <- pred_subset %>%
            dplyr::mutate("Index" = index[pred$Output_ID == current_output_id])
        }

        gg <- ggplot2::ggplot() +
          ggplot2::geom_raster(
            data = pred_subset,
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
      }

      # Add data points
      if (!is.null(data_subset)) {
        data_subset <- data_subset %>%
          tidyr::pivot_wider(
            names_from = Input_ID,
            values_from = Input,
            names_prefix = "Input_"
          ) %>%
          # Keep 6 significant digits for Inputs to avoid numerical issues
          dplyr::mutate(across(starts_with("Input_"), ~ round(.x, 6))) %>%
          rowwise() %>%
          dplyr::mutate(
            Reference = paste(
              # Create output's prefix
              paste0("o", Output_ID),
              # Create the reference for each Output_ID
              paste(c_across(starts_with("Input_")), collapse = ":"),
              # Join output's prefix and reference
              sep = ";"
            )
          ) %>%
          dplyr::ungroup()

        gg <- gg + ggplot2::geom_label(
          data = data_subset,
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
      if (ncol(inputs) == 1) {
        if ((dplyr::n_distinct(inputs) != nrow(inputs)) & is.null(index)) {
          warning("Some values on the x-axis appear multiple times... ")
        }
      } else {
        warning("Impossible to display inputs with dimensions greater than 2... ")
        inputs <- inputs %>% dplyr::select(Input)
      }

      ## Display a 'heatmap' if the argument is TRUE
      if (heatmap) {
        if (is.null(y_grid)) {
          y_grid <- seq(
            min(pred_subset$Mean) - quant_ci * sqrt(max(pred_subset$Var)),
            max(pred_subset$Mean) + quant_ci * sqrt(max(pred_subset$Var)),
            length.out = 500
          )
        }
        col_to_nest <- c(names(inputs)[1], "Mean", "Var")

        if (!is.null(index)) {
          pred_subset <- pred_subset %>%
            dplyr::mutate("Index" = index[pred$Output_ID == current_output_id])
          col_to_nest <- c(col_to_nest, "Index")
        }

        db_heat <- pred_subset %>%
          tidyr::expand(tidyr::nesting(!!!rlang::syms(col_to_nest)),
                        "Ygrid" = y_grid
          ) %>%
          dplyr::mutate("Proba" = 2 *
                          (1 - stats::pnorm(abs((data$Ygrid - data$Mean) / sqrt(data$Var)))))

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
            colours = c("white", "#FDE0DD", "#FCC5C0", "#FA9FB5",
                        "#F768A1", "#DD3497", "#AE017E", "#7A0177")
          ) +
          ggplot2::labs(fill = "Proba CI") +
          ggplot2::ylab("Output")

      } else if (samples){
        gg <- ggplot2::ggplot() +
          ggplot2::geom_line(
            data = samples_subset,
            ggplot2::aes_string(
              x = names(inputs)[1],
              y = "Output",
              group = "Sample"
            ),
            color = "#FA9FB5",
            alpha = alpha_samples
          ) +
          ggplot2::guides(group = "none")

        if(plot_mean){
          if(is.null(pred_subset)){
            warning(paste("For Output_ID '", current_output_id,
                          "': 'pred' is necessary to display the mean."))
          } else {
            gg = gg + ggplot2::geom_line(
              data = pred_subset,
              ggplot2::aes_string(x = names(inputs)[1], y = "Mean"),
              color = "#DB15C1"
            )
          }
        }

      } else {
        if (!is.null(index)) {
          pred_subset <- pred_subset %>%
            dplyr::mutate("Index" = index[pred$Output_ID == current_output_id])
        }

        gg <- ggplot2::ggplot() +
          ggplot2::geom_line(
            data = pred_subset,
            ggplot2::aes_string(x = names(inputs)[1], y = "Mean"),
            color = "#DB15C1"
          ) +
          ggplot2::geom_ribbon(
            data = pred_subset,
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

      ## Add training data points
      if (!is.null(data_train_subset)) {
        data_train_subset <- data_train_subset %>%
          tidyr::pivot_wider(
            names_from = Input_ID,
            values_from = Input,
            names_prefix = "Input_"
          ) %>%
          # Keep 6 significant digits for Inputs to avoid numerical issues
          dplyr::mutate(across(starts_with("Input_"), ~ round(.x, 6))) %>%
          rowwise() %>%
          dplyr::mutate(
            Reference = paste(
              # Create output's prefix
              paste0("o", Output_ID),
              # Create the reference for each Output_ID
              paste(c_across(starts_with("Input_")), collapse = ":"),
              # Join output's prefix and reference
              sep = ";"
            )
          ) %>%
          dplyr::ungroup()
        gg <- gg + ggplot2::geom_point(
          data = data_train_subset,
          ggplot2::aes_string(
            x = names(inputs)[1],
            y = "Output",
            col = "Task_ID"
          ),
          size = size_data_train,
          alpha = alpha_data_train
        ) + ggplot2::guides(color = "none")
      }
      if (!is.null(data_subset)) {
        data_subset <- data_subset %>%
          tidyr::pivot_wider(
            names_from = Input_ID,
            values_from = Input,
            names_prefix = "Input_"
          ) %>%
          # Keep 6 significant digits for Inputs to avoid numerical issues
          dplyr::mutate(across(starts_with("Input_"), ~ round(.x, 6))) %>%
          rowwise() %>%
          dplyr::mutate(
            Reference = paste(
              # Create output's prefix
              paste0("o", Output_ID),
              # Create the reference for each Output_ID
              paste(c_across(starts_with("Input_")), collapse = ":"),
              # Join output's prefix and reference
              sep = ";"
            )
          ) %>%
          dplyr::ungroup()

        gg <- gg + ggplot2::geom_point(
          data = data_subset,
          ggplot2::aes_string(x = names(inputs)[1], y = "Output"),
          size = size_data,
          shape = 20
        )
      }

      ## Add prior mean
      if (!is.null(prior_mean_subset)) {
        if (names(inputs)[1] %in% names(prior_mean_subset)) {
          gg <- gg +
            ggplot2::geom_line(
              data = prior_mean_subset,
              ggplot2::aes_string(x = names(inputs)[1], y = "Output"),
              linetype = "dashed"
            )
        } else {
          warning(
            "The ", names(inputs)[1], " column does not exist in 'prior_mean'."
          )
        }
      }
    }

    gg <- gg + ggplot2::theme_classic() +
      ggplot2::labs(title = paste("Prediction for Output_ID:", current_output_id))

    return(gg)

  }) %>%
    stats::setNames(unique_outputs)

  return(plot_list)
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
#'    \code{Sample}, \code{Output_ID}, \code{Output}.  This argument is needed if \code{pred}
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
#' @return List of graph of samples drawn from a posterior distribution of a GP,
#'    Magma, or MagmaClust prediction.
#' @export
#'
#' @examples
#' TRUE
plot_samples <- function(pred = NULL,
                         samples = NULL,
                         nb_samples = 50,
                         x_input = NULL,
                         plot_mean = TRUE,
                         alpha_samples = 0.3
                         ) {
  ## Check whether 'samples' or 'pred' exist
  if(is.null(samples) & is.null(pred) ){
      stop("Either 'sample' or 'pred' is needed as an argument.")
  }

  ## If provided, check format of 'pred' and extract the mixture prediction
  if(!is.null(pred)){

    ## Check 'pred' format
    if (!(is.list(pred) & ('cov' %in% names(pred))) ) {
      stop(
        "The 'pred' argument should be a list containing 'pred' and 'cov' ",
        "elements. Consider re-running the prediction function using the ",
        "argument 'get_full_cov' = TRUE."
      )
    }

    ## Check whether 'pred' is a GP/Magma or a MagmaClust prediction
    if ('mixture_pred' %in% names(pred)){

      ## If 'samples' is not provided, draw new samples
      if(is.null(samples)){
        samples = sample_magmaclust(pred_clust = pred, nb_samples = nb_samples)
      }

      mean_pred = pred$mixture_pred

    } else {

      ## If 'samples' is not provided, draw new samples
      if(is.null(samples)){
        samples = sample_gp(pred_gp = pred, nb_samples = nb_samples)
      }

      mean_pred = pred$pred
    }
  }

  if (!"Output_ID" %in% names(samples)) {
    stop("'samples' tibble should contain 'Output_ID' column.")
  }

  ## Get the unique ID of Outputs
  unique_outputs <- unique(samples$Output_ID)

  ## Create a plot for each Output_ID
  plot_list <- purrr::map(unique_outputs, function(current_output_id) {

    # Subset samples only on the current Output_ID
    samples_subset <- samples %>%
      dplyr::filter(Output_ID == current_output_id)

    mean_pred_subset <- NULL
    if (plot_mean && !is.null(pred) && exists("mean_pred") && !is.null(mean_pred)) {
      mean_pred_subset <- mean_pred %>%
        dplyr::filter(Output_ID == current_output_id)
    }

    # Determine inputs on which we want to plot (subset relative to the current Output_ID)
    if (x_input %>% is.null()) {
      inputs <- samples_subset %>% dplyr::select(-c(Sample, Output, Output_ID))
    } else {
      inputs <- samples_subset[x_input]
    }

    ## Display a heatmap if inputs are 2D
    if (ncol(inputs) == 2) {
      ## Extract only one sample when displaying in 2D
      samples_2d <- samples_subset %>%
        dplyr::filter(Sample == unique(samples_subset$Sample)[1]) %>%
        dplyr::select(- Sample)

      gg <- ggplot2::ggplot() +
        ggplot2::geom_raster(
          data = samples_2d,
          ggplot2::aes_string(
            x = names(inputs)[1],
            y = names(inputs)[2],
            fill = "Output"
          ),
          interpolate = TRUE
        ) +
        ggplot2::scale_fill_distiller(palette = "RdPu", trans = "reverse")

    } else if (ncol(inputs) == 1) {
      ## Plot 1D
      gg <- ggplot2::ggplot() +
        ggplot2::geom_line(
          data = samples_subset,
          ggplot2::aes_string(
            x = names(inputs)[1],
            y = "Output",
            group = "Sample"
          ),
          color = "#FA9FB5",
          alpha = alpha_samples
        ) +
        ggplot2::guides(group = "none")

      if(plot_mean){
        if(is.null(pred) || is.null(mean_pred_subset)){
          warning(paste("For Output_ID '", current_output_id,
                        "': 'pred' is necessary to compute the mean"))
        } else {
          gg = gg + ggplot2::geom_line(
            data = mean_pred_subset,
            ggplot2::aes_string(x = names(inputs)[1], y = "Mean"),
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

    gg <- gg +
      ggplot2::theme_classic() +
      ggplot2::labs(title = paste("Prediction for Output_ID:", current_output_id))

    # Return plot for the current Output_ID
    return(gg)

  }) %>%
    stats::setNames(unique_outputs) # Name plots in the list

  ## Return list of plots
  return(plot_list)
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
plot_gif <- function(pred_gp,
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
                     ...) {
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

  gg %>% return()
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
plot_magmaclust <- function(pred_clust,
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
        "The mixture probability of the cluster", cluster, "is 1. Therefore,",
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

  if(samples) {
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
      } else if (ncol(inputs) == 1){
        ## Display points if 1-D
        gg <- gg + ggplot2::geom_point(
          data = data,
          ggplot2::aes_string(x = names(inputs)[1], y = "Output"),
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
            ggplot2::aes_string(
              x = names(inputs)[1],
              y = "Output", col = "Cluster"
            ),
            size = size_data_train,
            alpha = alpha_data_train
          )
      }
    } else {
      gg <- gg +
        ggplot2::geom_point(
          data = data_train,
          ggplot2::aes_string(
            x = names(inputs)[1],
            y = "Output", fill = "ID"
          ),
          shape = 21,
          size = size_data_train,
          alpha = alpha_data_train
        ) + ggplot2::guides(fill = "none")
    }
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
          ggplot2::aes_string(
            x = names(inputs)[1],
            y = "Output",
            col = "Cluster"
          ),
          linetype = "dashed"
        )
    } else {
      warning(
        "The ", names(inputs)[1], " column does not exist in the ",
        "'prior_mean' argument. The mean function cannot be displayed."
      )
    }
  }

  ## Change scale colour palette
  gg <- gg + ggplot2::scale_color_brewer(palette = "Set1")

  (gg + ggplot2::ggtitle(gtitle)) %>%
    return()
}
