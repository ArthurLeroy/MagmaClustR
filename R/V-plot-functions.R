#' Plot smoothed raw data
#'
#' @param db database with format : ID, Input, Output
#' @param cluster boolean indicating whether you want cluster or not
#' @param legend boolean indicating whether want legend or not
#'
#' @return TRUE
#' @export
#'
#' @examples
V_plot_db = function(db, cluster = F, legend = F)
{
  # if(cluster)
  # {
  #   ggplot2::ggplot(db) + ggplot2::geom_smooth(ggplot2::aes(Input, Output, group = ID, color = Cluster)) +
  #     ggplot2::geom_point(ggplot2::aes(Input, Output, group = ID, color = Cluster)) + ggplot2::guides(col = legend)
  # }
  # else
  # {
  #   ggplot2::ggplot(db) + ggplot2::geom_smooth(ggplot2::aes(Input, Output, color = ID)) +
  #     ggplot2::geom_point(ggplot2::aes(Input, Output, color = ID)) + ggplot2::guides(col = legend)
  # }
}

#' plot_animate_clust
#'
#' @param pred_gp pred_gp
#' @param ygrid ygrid
#' @param data data
#' @param data_train data_train
#' @param mean_k mean_k
#' @param file file
#'
#' @return
#' TRUE
#' @export
#'
#' @examples
plot_animate_clust = function(pred_gp, ygrid, data = NULL, data_train = NULL, mean_k = NULL, file = "gganim.gif")
{ ## pred_gp : tibble coming out of the pred_gp_animate() function, columns required : 'Timestamp', 'Mean', 'Var'
  ## data : tibble of observational data, columns required : 'Timestamp', 'Output' (Optional)
  ####
  ## return : plot the animated curves of the GP with the 0.95 confidence interval (optional display raw data)

  # gg = plot_heat_clust(pred_gp, ygrid, data, data_train, mean_k, col_clust = F, animate = T) +
  #   transition_states(Nb_data, transition_length = 2, state_length = 1)
  # animate(gg, renderer = gifski_renderer(file)) %>% return()
}
