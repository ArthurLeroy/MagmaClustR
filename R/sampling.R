#' Draw samples from a posterior GP/Magma distribution
#'
#' @param pred_gp A list, typically coming from
#'    \code{\link{pred_magma}} or \code{\link{pred_gp}} functions, with argument
#'    'get_full_cov = TRUE'. Required elements: \code{pred}, \code{cov}.
#' @param nb_samples A number, indicating the number of samples to be drawn from
#'    the predictive posterior distribution. For two-dimensional graphs, only
#'    one sample can be displayed.
#'
#' @return A tibble or data frame, containing the samples generated from
#'    a GP prediction. Format: \code{Input}, \code{Sample}, \code{Output_ID},
#'    \code{Output}.
#' @export
#'
#' @examples
#' TRUE
sample_gp = function(
    pred_gp,
    nb_samples = 50){
  ## Extract the GP prediction
  pred <- pred_gp$pred

  ## Remove 'ID' if present
  if ("Task_ID" %in% names(pred)) {
    pred = pred %>% dplyr::select(-Task_ID)
  }

  ## Extract the full covariance matrix
  cov <- pred_gp$cov

  if (!"Output_ID" %in% names(pred)) {
    stop("'pred_gp' should contain an 'Output_ID' column.")
  }

  ## Get the unique ID of Outputs
  unique_outputs <- unique(pred$Output_ID)

  ## Sample on each output and combine results
  all_samples <- purrr::map_dfr(unique_outputs, function(current_output_id) {
    # Subset pred only on the current Output_ID
    output_indices <- which(pred$Output_ID == current_output_id)

    if (length(output_indices) == 0) {
      # Do nothing if the current Output_ID does not have any prediction point
      return(NULL)
    }

    # Filter the pred dataframe on the current Output_ID indices
    pred_subset <- pred[output_indices, ]

    # Subset mean on the current Output_ID
    mean_subset <- pred_subset %>% dplyr::pull(Mean)

    # Extract corresponding inputs
    inputs_subset <- pred_subset %>% dplyr::select(-c(Mean, Var, Output_ID))

    # Subset cov on the current Output_ID thanks to inputs_subset
    cov_subset <- cov[output_indices, output_indices, drop = FALSE]

    # Draw samples for the current Output_ID
    drawn_samples <- mvtnorm::rmvnorm(n = nb_samples,
                                      mean = mean_subset,
                                      sigma = cov_subset,
                                      checkSymmetry = FALSE)

    # Organize tibble
    tidied_samples <- drawn_samples %>%
      t() %>%
      magrittr::set_colnames(paste0("s_", 1:nb_samples)) %>%
      tibble::as_tibble() %>%
      dplyr::bind_cols(inputs_subset) %>%
      tidyr::pivot_longer(
        cols = dplyr::starts_with("s_"),
        names_to = "Sample",
        names_prefix = "s_",
        values_to = "Output"
      ) %>%
      dplyr::mutate(Output_ID = current_output_id)

    return(tidied_samples)
  })

  return(all_samples)
}

#' @rdname sample_gp
#' @export
sample_magma <- sample_gp


#' Draw samples from a MagmaClust posterior distribution
#'
#' @param pred_clust A list, typically coming from
#'    \code{\link{pred_magmaclust}}, with argument get_full_cov = TRUE'.
#'    Required elements: \code{pred}, \code{cov}, \code{mixture}.
#' @param nb_samples A number, indicating the number of samples to be drawn from
#'    the predictive posterior distribution. For two-dimensional graphs, only
#'    one sample can be displayed.
#'
#' @return A tibble or data frame, containing the samples generated from
#'    a GP prediction. Format: \code{Cluster}, \code{Proba}, \code{Input},
#'    \code{Sample}, \code{Output}.
#' @export
#'
#' @examples
#' TRUE
sample_magmaclust = function(
    pred_clust,
    nb_samples = 50){

  floop = function(k){

  ## Extract the GP prediction
  pred <- pred_clust$pred[[k]]

  ## Remove 'ID' if present
  if ("ID" %in% names(pred)) {
    pred = pred %>% dplyr::select(- .data$ID)
  }

  ## Extract parameters and inputs from the prediction
  inputs <- pred %>% dplyr::select(-c(.data$Mean, .data$Var))
  mean <- pred %>% dplyr::pull(.data$Mean)
  cov <- pred_clust$cov[[k]]
  weight <- pred_clust$mixture[[k]]

  #Draw samples and format the tibble
  mvtnorm::rmvnorm(nb_samples, mean, cov, checkSymmetry = FALSE) %>%
    t() %>%
    magrittr::set_colnames(1:nb_samples) %>%
    tibble::as_tibble() %>%
    dplyr::bind_cols(inputs) %>%
    tidyr::pivot_longer(- names(inputs) ,
                        names_to= "Sample",
                        values_to = "Output") %>%
    dplyr::mutate('Proba' = weight,
                  'Cluster' = k, .before = 1) %>%
    return()
  }
  names(pred_clust$pred) %>%
    lapply(floop) %>%
    dplyr::bind_rows() %>%
    dplyr::group_by(.data$Sample, .data$Input) %>%
    dplyr::summarise('Output' = sum(.data$Proba * .data$Output),
                     .groups = "drop") %>%
    return()
}
