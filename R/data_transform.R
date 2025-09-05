# sparse_batch_data(data) = function(data,
#                                    size_grid,
#                                    grid_inputs = NULL){
#   data %>%
#     dplyr::mutate(Batch_input = grid_inputs[match_closest(db$Input, grid_inputs)]) %>%
#     return()
# }

laplace_matching = function(data, likelihood = 'Beta', eps = 1e-1){

  ## Initialise a dummy Batch column if inputs are not batched
  if( !('Batch_input' %in% names(data)) ){
    data = data %>% dplyr::mutate(Batch_input = .data$Input)
  }

  if(likelihood == 'Beta'){

    ## Compute alpha et beta parameters of the Beta pseudo-likelihood
    ## Then transform into Gaussian observations via Laplace-Matching
    pseudo_data = data %>%
      dplyr::group_by(.data$ID, .data$Batch_input) %>%
      dplyr::mutate(
             alpha = sum(.data$Output == T),
             beta = sum(.data$Output == F)) %>%
      dplyr::mutate(Output = log((.data$alpha + eps) / (.data$beta + eps))) %>%
      dplyr::mutate(Var_output=
                      (.data$alpha+.data$beta+2*eps) /
                      ((.data$alpha+eps)*(.data$beta+eps))) %>%
      dplyr::select(c(.data$ID, .data$Batch_input, .data$Output)) %>%
      dplyr::rename(Input = .data$Batch_input) %>%
      dplyr::distinct() %>%
      dplyr::ungroup()
  }

  return(pseudo_data)
}

revert_laplace_matching = function(sample, likelihood = 'Beta'){

  if(likelihood == 'Beta'){
    ## Apply the logistic transform to revert back to [0,1]
    revert_pred = sample %>%
      dplyr::mutate(Output = 1/(1+exp(-.data$Output)))
  }

  return(revert_pred)
}
