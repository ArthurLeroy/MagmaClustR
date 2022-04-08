#' Plot MagmaClust or GP predictions
#'
#' @param pred A tibble or data frame, typically coming from
#' \code{\link{pred_magmaclust}} function.
#'    Required columns: 'Input', 'Mean', 'Var'. Additional covariate columns may be
#'    present in case of multi-dimensional inputs.
#' @param cluster A string indicating the cluster to plot from or
#' 'all' for the full GPs mixture.
#' @param data tibble of observational data, columns required : 'Input', 'Output'
#' @param data_train tibble of training dataset, columns required : 'Input', 'Output'
#' @param col_clust A boolean indicating whether we color according to
#' clusters or individuals.
#' @param prob_CI A number between 0 and 1 (default is 0.95), indicating the
#'    level of the Credible Interval associated with the posterior mean curve.
#' @param x_input A vector of character strings, indicating which input should
#'    be displayed. If NULL(default) the 'Input' column is used for the x-axis.
#'    If providing a 2-dimensional vector, the corresponding columns are used
#'    for the x-axis and y-axis.
#' @param y_grid A vector, indicating the grid of values on the y-axis for which
#'    probabilities should be computed for heatmaps of 1-dimensional
#'    predictions.
#' @param heatmap A logical value indicating whether the GP prediction should be
#'    represented as a heatmap of probabilities for 1-dimensional inputs. If
#'    FALSE (default), the mean curve and associated 95%CI are displayed.
#' @param prior_mean (Optional) A tibble or a data frame, containing the 'Input'
#'    and associated 'Output' prior mean parameter of the GP prediction.
#'
#' @return Plot of the predicted curve of the GP with the 0.95 confidence interval
#' @export
#'
#' @examples
#'TRUE
#'
plot_magmaclust = function(pred,
                         x_input = NULL,
                         cluster = 'all',
                         data = NULL,
                         data_train = NULL,
                         col_clust = FALSE,
                         y_grid = NULL,
                         prior_mean = NULL,
                         heatmap = FALSE,
                         prob_CI = 0.95)
{
  ## Compute the quantile of the desired Credible Interval
  quant_ci <- stats::qnorm((1 + prob_CI) / 2)

  ## Check whether 'pred' has a correct format
  if (pred %>% is.data.frame()) {
    pred_gp <- pred
  } else if (is.list(pred) &
      tryCatch(is.data.frame(pred$Prediction), error = function(e) FALSE)) {
    pred_gp <- pred_gp$Prediction
  } else {
    stop(
      "The 'pred_gp' argument should either be a list containing the 'pred' ",
      "element or a data frame. Please read ?plot_magmaclust(), and use ",
      "pred_magmaclust() for making predictions under a correct format."
    )
  }


  if(cluster == 'all'){
    mean_all = 0
    var_all = 0
  for(k in names(pred_gp))
  {
    mean_all = mean_all + pred_gp[[k]]$hp_k_mixture * pred_gp[[k]]$Mean
    var_all = var_all + pred_gp[[k]]$hp_k_mixture * pred_gp[[k]]$Var
  }
  #browser()
  pred = tibble::tibble(pred_gp[[1]]%>%
                          dplyr::select(-c(.data$Mean,.data$Var,.data$hp_k_mixture)),
                        'Mean' = mean_all,
                        'Var' = var_all)
  }
  else{pred = pred_gp[[cluster]]}

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
  }
  else {
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
      data <- data %>% dplyr::mutate(Output = round(.data$Output, 1))

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
      col_to_nest = c(names(inputs)[1], "Mean", "Var")

      ## Add the 'Index' column if the prediction comes from 'pred_gif()'
      if (!is.null(index)) {
        pred <- pred %>% dplyr::mutate("Index" = index)
        col_to_nest = c(col_to_nest, 'Index')
      }

      db_heat <- pred %>%
        tidyr::expand(
          tidyr::nesting_(col_to_nest),
          "Ygrid" = y_grid
        ) %>%
        dplyr::mutate(
          "Proba" = 2 * (1 - stats::pnorm(
            abs((.data$Ygrid - .data$Mean) / sqrt(.data$Var)))
          )
        )

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
            "white", "#FDE0DD", "#FCC5C0", "#FA9FB5", "#F768A1",
            "#DD3497", "#AE017E", "#7A0177"
          )
        ) +
        ggplot2::labs(fill = "Proba CI") +
        ggplot2::ylab("Output")

    } else { ## Display a classic curve otherwise
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
          fill = '#FA9FB5'
        ) +
        ggplot2::ylab("Output")
    }

    ## Display the training data if provided
    if (!is.null(data_train)) {
      gg <- gg + ggplot2::geom_point(
        data = data_train,
        ggplot2::aes_string(x = names(inputs)[1], y = "Output", col = "ID"),
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
      colori <- 1
      for(k in prior_mean %>% names()){
        if (names(inputs)[1] %in% names(prior_mean[[k]])) {
          gg <- gg +
            ggplot2::geom_line(
              data = prior_mean[[k]],
              ggplot2::aes_string(x = names(inputs)[1], y = "Output"),
              linetype = "dashed",
              color = colori
              )
          } else {
            warning(
              "The ", names(inputs)[1], " column does not exist in the ",
              "'prior_mean' argument. The mean function cannot be displayed."
        )
          }
        colori <- colori + 1
      }

    }
  }

  (gg +  ggplot2::theme_classic()) %>%
    return()
}

