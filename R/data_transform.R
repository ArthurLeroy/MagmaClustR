laplace_matching = function(data, likelihood, eps = 1e-1) {

  if (likelihood == 'Bernoulli') {
    ## Compute alpha et beta parameters of the Beta pseudo-likelihood
    ## Then transform into Gaussian observations via Laplace-Matching
    pseudo_data = data %>%
      dplyr::mutate(
        alpha = .data$Output == 1,
        beta = .data$Output == 0
      ) %>%
      dplyr::mutate(Output = log((.data$alpha + eps) / (.data$beta + eps))) %>%
      dplyr::mutate(
        Var_output = (.data$alpha + .data$beta + 2 * eps) /
          ((.data$alpha + eps) * (.data$beta + eps))
      ) %>%
      dplyr::select(c(.data$ID, .data$Input, .data$Output, .data$Var_output))
  }

  return(pseudo_data)
}

revert_laplace_matching = function(sample, likelihood = 'Bernoulli') {
  if (likelihood == 'Beta') {
    ## Apply the logistic transform to revert back to [0,1]
    revert_pred = sample %>%
      dplyr::mutate(Output = 1 / (1 + exp(-.data$Output)))
  }

  return(revert_pred)
}
