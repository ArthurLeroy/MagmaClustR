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
#'    a GP prediction. Format: \code{Input}, \code{Sample}, \code{Output}.
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
  if ("ID" %in% names(pred)) {
    pred = pred %>% dplyr::select(- .data$ID)
  }

  ## Extract parameters and inputs from the prediction
  inputs <- pred %>% dplyr::select(-c(.data$Mean, .data$Var))
  mean <- pred %>% dplyr::pull(.data$Mean)
  cov <- pred_gp$cov

  #Draw samples and format the tibble
  mvtnorm::rmvnorm(nb_samples, mean, cov) %>%
    t() %>%
    tibble::as_tibble(.name_repair = 'unique') %>%
    dplyr::bind_cols(inputs) %>%
    tidyr::pivot_longer(- names(inputs) ,
                        names_to= "Sample",
                        values_to = "Output") %>%
    return()
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
  mvtnorm::rmvnorm(nb_samples, mean, cov) %>%
    t() %>%
    tibble::as_tibble(.name_repair = 'unique') %>%
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
    dplyr::summarise('Output' = sum(.data$Proba * .data$Output)) %>%
    dplyr::ungroup() %>%
    return()
}
