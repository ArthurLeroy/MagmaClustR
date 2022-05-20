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
#'    individual or a cluster, see \code'{col_clust} below).
#' @param col_clust A boolean indicating whether backward points are coloured
#'    according to the individuals or to their predicted cluster. If one wants
#'    to colour by clusters, a column \code{Cluster} shall be present in
#'    \code{data_train} and its value should refer to the most probable cluster
#'    for each individual.
#' @param prior_mean (Optional) A list providing, for each cluster, a
#'    tibble containing prior mean parameters of the prediction. This argument
#'    typically comes as an outcome \code{hyperpost$mean}, available through
#'    the \code{\link{train_magmaclust}}, \code{\link{pred_magmaclust}}
#'    functions.
#' @param y_grid A vector, indicating the grid of values on the y-axis for which
#'    probabilities should be computed for heatmaps of 1-dimensional
#'    predictions. If NULL (default), a vector of length 50 is defined, ranging
#'    between the min and max 'Output' values contained in \code{pred}.
#' @param heatmap A logical value indicating whether the GP prediction should be
#'    represented as a heatmap of probabilities for 1-dimensional inputs. If
#'    FALSE (default), the mean curve (and associated Credible Interval if
#'    available) are displayed.
#' @param prob_CI A number between 0 and 1 (default is 0.95), indicating the
#'    level of the Credible Interval associated with the posterior mean curve.
#'    If this this argument is set to 1, the Credible Interval is not displayed.
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
#'TRUE
#'
plot_magmaclust = function(pred_clust,
                         cluster = 'all',
                         x_input = NULL,
                         data = NULL,
                         data_train = NULL,
                         col_clust = FALSE,
                         prior_mean = NULL,
                         y_grid = NULL,
                         heatmap = FALSE,
                         prob_CI = 0.95)
{
  ## Check prob_CI format
  if(prob_CI < 0 | prob_CI > 1)
  {
    stop("The 'prob_CI' argument should be a number between 0 and 1.")
  }

  ## Check format for prediction
  if(!is.list(pred_clust)){
    stop("Wrong format for 'pred_clust', please read ?plot_magmaclust().")
  }
  ## Check presence of 'pred'
  if( !('pred' %in% names(pred_clust)) ){
    stop("Wrong format for 'pred_clust', please read ?plot_magmaclust().")
  }
  ## Check presence of 'mixture'
  if( !('mixture' %in% names(pred_clust)) ){
    stop("Wrong format for 'pred_clust', please read ?plot_magmaclust().")
  }
  ## Check presence of 'pred_mixture'
  if( !('mixture_pred' %in% names(pred_clust)) ){
    stop("Wrong format for 'pred_clust', please read ?plot_magmaclust().")
  }

  pred = pred_clust$pred
  mixture = pred_clust$mixture
  mixture_pred = pred_clust$mixture_pred
  ID_k = names(pred)

  ## Checker whether we can provide Credible Interval
  if(cluster == 'all'){
    ## Check whether one cluster's proba is 1 (or really close)
    max_clust = proba_max_cluster(mixture)
    if( round(max_clust$Proba, 3) == 1){

        cluster = max_clust$Cluster
        cat(
          "The mixture probability of the cluster", cluster,"is 1. Therefore,",
          "the predictive distribution is Gaussian and the associated",
          "credible interval can be displayed. \n\n"
        )
      ## Create dummy variable for indicating the type of prediction
      all_clust = FALSE

      ## Compute the quantile of the desired Credible Interval
      quant_ci <- stats::qnorm((1 + prob_CI) / 2)
    } else {
      all_clust = TRUE
    }

  } else {
    all_clust = FALSE
    ## Compute the quantile of the desired Credible Interval
    quant_ci <- stats::qnorm((1 + prob_CI) / 2)
  }

  ## Select the appropriate tibble for displaying predictions
  if( all_clust ){
    pred_gp = mixture_pred

  } else {
    ## Check the name provided in 'cluster'
    if( !(cluster %in% ID_k) ){
      stop("The cluster's name provided in 'cluster' does not exist in) " ,
           "'pred_clust'.")
    }
    ## Remove the 'Proba' column if selecting cluster-specific prediction
    pred_gp <- pred[[cluster]] %>% dplyr::select(- .data$Proba)

    ## Get the 'Proba' value to display in the Title
    proba = pred[[cluster]] %>%
      dplyr::pull(.data$Proba) %>%
      unique()
  }

  ## Get the inputs that should be used
  if (x_input %>% is.null()) {
    inputs <- pred_gp %>% dplyr::select(-c(.data$ID, .data$Mean, .data$Var))
  } else {
    if( all(x_input %in% names(pred_gp)) ){
      inputs <- pred_gp[x_input]
    } else {
      stop("The names in the 'x_input' argument don't exist in 'pred_clust'.")
    }
  }

  if(all_clust){
    ## GP visualisation without Credible Interval
    gg <- plot_gp(
      pred_gp = pred_gp,
      x_input = x_input,
      data = data,
      data_train = data_train,
      y_grid = y_grid,
      heatmap = heatmap,
      prob_CI = 0)


    gtitle = paste0("Mixture of GP predictions")
  } else {
    if(col_clust)
    {
      ## Check whether 'data_train' provides a 'Cluster' column
      if(!('Cluster' %in% names(data_train))){
        cat("The 'data_train' argument does not provide a 'Cluster' column.",
            "Therefore, training data remain coloured by individual. \n \n")
      } else{
        ## Certify that 'Cluster' is discrete
        data_train$Cluster = as.factor(data_train$Cluster)

        ## Classic GP visualisation with modified training data plot
        gg <- plot_gp(
          pred_gp = pred_gp,
          x_input = x_input,
          data = data,
          y_grid = y_grid,
          heatmap = heatmap,
          prob_CI = prob_CI
          ) +
          ggplot2::geom_point(
            data = data_train,
            ggplot2::aes_string(x = names(inputs)[1],
                                y = "Output", col = "Cluster"),
            size = 0.5,
            alpha = 0.5
          ) +
          ggplot2::scale_color_brewer(palette="Set1")
      }
    } else{
      ## Classic GP visualisation for cluster-specific predictions
      gg <- plot_gp(
        pred_gp = pred_gp,
        x_input = x_input,
        data = data,
        data_train = data_train,
        y_grid = y_grid,
        heatmap = heatmap,
        prob_CI = prob_CI)
    }

    gtitle = paste0("Cluster ", cluster, " -- Proba = ", mixture[[cluster]])
  }

  ## Display the prior mean process if provided
  if (!is.null(prior_mean)) {
    ## Bind the tibbles of hyper-posterior mean processes
    mean_k = prior_mean %>% dplyr::bind_rows(.id = 'Cluster')

    ## Display the mean functions if available
    if (names(inputs)[1] %in% names(mean_k)) {
      gg <- gg +
        ggplot2::geom_line(
          data = mean_k,
          ggplot2::aes_string(
            x = names(inputs)[1],
            y = 'Output',
            col = 'Cluster'),
          linetype = 'dashed'
        ) +
        ggplot2::scale_color_brewer(palette="Set1")
    } else {
      warning(
        "The ", names(inputs)[1], " column does not exist in the ",
        "'prior_mean' argument. The mean function cannot be displayed."
      )
    }
  }


  (gg + ggplot2::ggtitle(gtitle)) %>%
  return()
}
#
#   ## Get the inputs that should be used
#   if (x_input %>% is.null()) {
#     inputs <- pred_gp %>%
#       dplyr::select(-c(.data$ID, .data$Mean, .data$Var))
#   }
#   else {
#     inputs <- pred_gp[x_input]
#   }
#
#   if(!all_clust){
#     ## Format the tibble for displaying the Credible Intervals
#     pred_gp <- pred_gp %>%
#       dplyr::mutate("CI_inf" = .data$Mean - quant_ci * sqrt(.data$Var)) %>%
#       dplyr::mutate("CI_sup" = .data$Mean + quant_ci * sqrt(.data$Var)) %>%
#       dplyr::mutate("CI_Width" = .data$CI_sup - .data$CI_inf)
#   }
#
#
#   ## Display a heatmap if inputs are 2D
#   if (ncol(inputs) == 2) {
#
#     gg <- ggplot2::ggplot() +
#       ggplot2::geom_raster(
#         data = pred_gp,
#         ggplot2::aes_string(
#           x = names(inputs)[1],
#           y = names(inputs)[2],
#           fill = "Mean",
#           alpha = "CI_Width"
#         ),
#         interpolate = TRUE
#       ) +
#       ggplot2::scale_fill_distiller(palette = "RdPu", trans = "reverse") +
#       ggplot2::scale_alpha_continuous(range = c(0.1, 1), trans = "reverse")
#
#     if (!is.null(data)) {
#       ## Round the 'Output' values to reduce size of labels on the graph
#       data <- data %>% dplyr::mutate(Output = round(.data$Output, 1))
#
#       gg <- gg + ggplot2::geom_label(
#         data = data,
#         ggplot2::aes_string(
#           x = names(inputs)[1],
#           y = names(inputs)[2],
#           label = "Output",
#           fill = "Output"
#         ),
#         size = 3
#       )
#     }
#   } else {
#     ## Check the dimension of the inputs a propose an adequate representation
#     if (ncol(inputs) == 1) {
#       if ((dplyr::n_distinct(inputs) != nrow(inputs)) ) {
#         warning(
#           "Some values on the x-axis appear multiple times, probably ",
#           "resulting in an incorrect graphical representation. Please ",
#           "consider recomputing predictions for more adequate inputs. "
#         )
#       }
#     } else {
#       warning(
#         "Impossible to display inputs with dimensions greater than 2. The ",
#         "graph then simply uses 'Input' as x_axis and 'Output' as y-axis. "
#       )
#       inputs <- inputs %>% dplyr::select(.data$Input)
#     }
#
#     ## Display a 'heatmap' if the argument is TRUE
#     if (heatmap) {
#       if (is.null(y_grid)) {
#         y_grid <- seq(
#           min(pred$Mean) - quant_ci * sqrt(max(pred$Var)),
#           max(pred$Mean) + quant_ci * sqrt(max(pred$Var)),
#           length.out = 500
#         )
#       }
#       ## Define the columns needed to compute a prediction for all y-axis values
#       col_to_nest = c(names(inputs)[1], "Mean", "Var")
#
#       db_heat <- pred %>%
#         tidyr::expand(
#           tidyr::nesting_(col_to_nest),
#           "Ygrid" = y_grid
#         ) %>%
#         dplyr::mutate(
#           "Proba" = 2 * (1 - stats::pnorm(
#             abs((.data$Ygrid - .data$Mean) / sqrt(.data$Var)))
#           )
#         )
#
#       gg <- ggplot2::ggplot() +
#         ggplot2::geom_raster(
#           data = db_heat,
#           ggplot2::aes_string(
#             x = names(inputs)[1],
#             y = "Ygrid",
#             fill = "Proba"
#           ),
#           interpolate = TRUE
#         ) +
#         ggplot2::scale_fill_gradientn(
#           colours = c(
#             "white", "#FDE0DD", "#FCC5C0", "#FA9FB5", "#F768A1",
#             "#DD3497", "#AE017E", "#7A0177"
#           )
#         ) +
#         ggplot2::labs(fill = "Proba CI") +
#         ggplot2::ylab("Output")
#
#     } else { ## Display a classic curve otherwise
#
#       gg <- ggplot2::ggplot() +
#         ggplot2::geom_line(
#           data = pred,
#           ggplot2::aes_string(x = names(inputs)[1], y = "Mean"),
#           color = "#DB15C1"
#         ) +
#         ggplot2::geom_ribbon(
#           data = pred,
#           ggplot2::aes_string(
#             x = names(inputs)[1],
#             ymin = "CI_inf",
#             ymax = "CI_sup"
#           ),
#           alpha = 0.2,
#           fill = '#FA9FB5'
#         ) +
#         ggplot2::ylab("Output")
#     }
#
#     ## Display the training data if provided
#     if (!is.null(data_train)) {
#       gg <- gg + ggplot2::geom_point(
#         data = data_train,
#         ggplot2::aes_string(x = names(inputs)[1], y = "Output", col = "ID"),
#         size = 0.5,
#         alpha = 0.5
#       ) +
#         ggplot2::guides(color = FALSE)
#     }
#     ## Display the observed data if provided
#     if (!is.null(data)) {
#       gg <- gg + ggplot2::geom_point(
#         data = data,
#         ggplot2::aes_string(x = names(inputs)[1], y = "Output"),
#         size = 2,
#         shape = 20
#       )
#     }
#   ## Display the (hyper-)prior mean process if provided
#     if (!is.null(prior_mean)) {
#       colori <- 1
#       for(k in prior_mean %>% names()){
#         if (names(inputs)[1] %in% names(prior_mean[[k]])) {
#           gg <- gg +
#             ggplot2::geom_line(
#               data = prior_mean[[k]],
#               ggplot2::aes_string(x = names(inputs)[1], y = "Output"),
#               linetype = "dashed",
#               color = colori
#               )
#           } else {
#             warning(
#               "The ", names(inputs)[1], " column does not exist in the ",
#               "'prior_mean' argument. The mean function cannot be displayed."
#         )
#           }
#         colori <- colori + 1
#       }
#
#     }
#   }
#
#   (gg +  ggplot2::theme_classic()) %>%
#     return()
# }
#
