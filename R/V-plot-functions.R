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
#' @param heatmap A logical value indicating whether the GP prediction should be
#'    represented as a heatmap of probabilities for 1-dimensional inputs. If
#'    FALSE (default), the mean curve (and associated Credible Interval if
#'    available) are displayed.
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
                         prob_CI = 0.95,
                         size_data = 3,
                         size_data_train = 1,
                         alpha_data_train = 0.5)
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
  ## Check presence of 'mixture_pred'
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
      y_grid = y_grid,
      heatmap = heatmap,
      prob_CI = 0,
      size_data = size_data,
      size_data_train = size_data_train,
      alpha_data_train = alpha_data_train)

    ## Define the adequate title
    gtitle = paste0("Mixture of GP predictions")
  } else {
    ## Classic GP visualisation for cluster-specific predictions
    gg <- plot_gp(
      pred_gp = pred_gp,
      x_input = x_input,
      data = data,
      y_grid = y_grid,
      heatmap = heatmap,
      prob_CI = prob_CI,
      size_data = size_data,
      size_data_train = size_data_train,
      alpha_data_train = alpha_data_train)

    ## Define the adequate title
    gtitle = paste0("Cluster ", cluster, " -- Proba = ", mixture[[cluster]])
  }

  ## Display training data if available
  if(!is.null(data_train)){
    ## Change colours of background points depending on 'col_clust'
    if(col_clust)
    {
      ## Check whether 'data_train' provides a 'Cluster' column
      if(!('Cluster' %in% names(data_train))) {
        cat("The 'data_train' argument does not provide a 'Cluster' column.",
            "Therefore, training data remain coloured by individual. \n \n")
      } else{
        ## Certify that 'Cluster' is discrete
        data_train$Cluster = as.factor(data_train$Cluster)

        ## Colour training data plot by cluster
        gg <- gg +
          ggplot2::geom_point(
            data = data_train,
            ggplot2::aes_string(x = names(inputs)[1],
                                y = "Output", col = "Cluster"),
            size = size_data_train,
            alpha = alpha_data_train
          )
      }
    } else {
      gg <- gg +
        ggplot2::geom_point(
          data = data_train,
          ggplot2::aes_string(x = names(inputs)[1],
                              y = "Output", fill = "ID"),
          shape = 21,
          size = size_data_train,
          alpha = alpha_data_train
        ) + ggplot2::guides(fill = 'none')
    }
  }

  ## Display the prior mean process if provided
  if (!is.null(prior_mean)) {

    ## Extract 'mean' if the user provides the full 'hyperpost'
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
        )
    } else {
      warning(
        "The ", names(inputs)[1], " column does not exist in the ",
        "'prior_mean' argument. The mean function cannot be displayed."
      )
    }
  }

  ## Change scale colour palette
  gg <- gg + ggplot2::scale_color_brewer(palette="Set1")

  (gg + ggplot2::ggtitle(gtitle)) %>%
  return()
}
