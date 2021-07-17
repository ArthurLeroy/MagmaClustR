#' Plot smoothed raw data
#'
#' @param db database with format : ID, Input, Output
#'
#' @return Visualize smoothed raw data.
#'
#' @examples
#' TRUE
plot_db = function(db)
{
  ggplot2::ggplot(db) +
    ggplot2::geom_smooth(ggplot2::aes(.data$Input,
                                      .data$Output,
                                      color = .data$ID), se = F) +
    ggplot2::geom_point(ggplot2::aes(.data$Input,
                                     .data$Output,
                                     color = .data$ID)) +
    ggplot2::theme_classic()
}

#' Plot your Gaussian Process
#'
#' @param pred_gp tibble coming out of the pred_gp() function, columns required : 'input', 'Mean', 'Var'
#' @param data tibble of observational data, columns required : 'input', 'Output' (Optional)
#' @param data_train tibble of training data, columns required : 'ID', input', 'Output' (Optional)
#' @param mean tibble of estimated value for the mean process, columns required : 'input', 'Mean' (Optional)
#' @param mean_CI boolean indicating whether the mean curve should be displayed along with CI95 band or not (Optional)
#'
#' @return plot the predicted curve of the GP with the 0.95 confidence interval (optional display raw data)
#' @export
#'
#' @examples
#' TRUE
plot_gp = function(pred_gp, data = NULL, data_train = NULL, mean = NULL, mean_CI = F)
{
  gg = ggplot2::ggplot() +
    ggplot2::geom_line(data = pred_gp, ggplot2::aes(x = pred_gp$input, y = pred_gp$Mean), color = 'blue') +
    ggplot2::geom_ribbon(data = pred_gp, ggplot2::aes(x = pred_gp$input, ymin = pred_gp$Mean - 1.96 * sqrt(pred_gp$Var),
                                                      ymax = pred_gp$Mean +  1.96 * sqrt(pred_gp$Var)), alpha = 0.2) + ggplot2::ylab('Output')

  ## Display the raw data and/or mean (with or without its CI) if provided
  if(!is.null(data_train)){gg = gg + ggplot2::geom_point(data = data_train, ggplot2::aes(x = data_train $input, y = data_train $Output, col = data_train$ID ),
                                                         size = 0.5, shape = 4)}
  if(!is.null(data)){gg = gg + ggplot2::geom_point(data = data, ggplot2::aes(x = data$input, y = data$Output), size = 2, shape = 18)}
  if(!is.null(mean)){gg = gg + ggplot2::geom_line(data = mean, ggplot2::aes(x = mean$input, y = mean$Mean), linetype = 'dashed')}
  if(mean_CI){gg = gg + ggplot2::geom_ribbon(data = mean, ggplot2::aes(x = mean$input, ymin = mean$Mean - 1.96 * sqrt(mean$Var),
                                                                       ymax = mean$Mean +  1.96 * sqrt(mean$Var)), alpha = 0.4)}

  return(gg)
}

#' Title
#'
#' @param pred_gp tibble coming out of the pred_gp() function, columns required : 'input', 'Mean', 'Var'
#' @param data tibble of observational data for the new individual, columns required : 'input', 'Output' (optional)
#' @param data_train tibble of training data for the model, , columns required : 'input', 'Output' (optional)
#' @param mean (optional)
#' @param ygrid grid on y-axis to plot the heatmap on (optional, default = range of data +- 2 sd)
#' @param interactive boolean indicating whether the output plot should be interactive (plotly)
#' @param CI boolean. If True, display the heatmap of Credible Intervals. If False, display the heatmap of likelihood
#'
#' @return A heatmap illustrating the GP predition results (Optionnal diplay of raw data and mean process)
#' @export
#'
#' @examples
#' TRUE
plot_heat =  function(pred_gp, data = NULL, data_train = NULL, mean = NULL, ygrid = NULL, interactive = F, CI = T)
{
  if(is.null(ygrid)){ygrid = seq(min(pred_gp$Mean) - 2 * sqrt(max(pred_gp$Var)) ,
                                 max(pred_gp$Mean) + 2 * sqrt(max(pred_gp$Var)), 0.1)}

  ## mutate from plotly or dplyr ?

  if(CI)
  {
    db_heat = tidyr::expand(pred_gp, tidyr::nesting(pred_gp$input, pred_gp$Mean, pred_gp$Var), 'Ygrid' = ygrid) %>%
      plotly::mutate('Proba' = 2 * stats::pnorm(abs((.data$Ygrid - pred_gp$Mean)/ sqrt(pred_gp$Var))) - 1)
    gg = ggplot2::ggplot(db_heat) + ggplot2::geom_tile(ggplot2::aes(pred_gp$input, .data$Ygrid, fill = .data$Proba)) + ggplot2::scale_fill_distiller(palette = "RdPu") +
      ggplot2::theme_minimal() + ggplot2::ylab("Output") + ggplot2::labs(fill = "Proba CI")
  }
  else
  {
    db_heat = tidyr::expand(pred_gp, tidyr::nesting(pred_gp$input, pred_gp$Mean, pred_gp$Var), 'Ygrid' = ygrid) %>%
      plotly::mutate('Proba' = stats::dnorm(.data$Ygrid, mean = pred_gp$Mean, sd = sqrt(pred_gp$Var)) )
    gg = ggplot2::ggplot(db_heat) + ggplot2::geom_tile(ggplot2::aes(pred_gp$input, .data$Ygrid, fill = .data$Proba)) +
      ggplot2::scale_fill_distiller(palette = "RdPu", trans = "reverse") +
      ggplot2::theme_minimal() + ggplot2::ylab("Output") + ggplot2::labs(fill = "Likelihood")
  }
  ## Display data/training data/mean process if provided
  if(!is.null(data_train)){gg = gg + ggplot2::geom_point(data = data_train, ggplot2::aes(x = data_train$input, y = data_train$Output, col = data_train$ID), shape = 4)}
  if(!is.null(data)){gg = gg + ggplot2::geom_point(data = data, ggplot2::aes(x = data$input, y = data$Output), size = 2, shape = 18)}
  if(!is.null(mean)){gg = gg + ggplot2::geom_line(data = mean, ggplot2::aes(x = mean$input, y = mean$Mean), linetype = 'dashed')}

  ## Turn into an interactive plotly (html plot)
  if(interactive){gg = plotly::ggplotly(gg)}

  return(gg)
}

#' Plot the animated (GIF) curves of the GP
#'
#' @param pred_gp tibble coming out of the pred_gp() function, columns required : 'input', 'Mean', 'Var'
#' @param data tibble of observational data, columns required : 'input', 'Output' (Optional)
#' @param data_train tibble of training data, columns required : 'ID', input', 'Output' (Optional)
#' @param mean tibble of estimated value for the mean process, columns required : 'input', 'Mean' (Optional)
#' @param mean_CI boolean indicating whether the mean curve should be displayed along with CI95 band or not (Optional)
#' @param file The animation file
#'
#' @return plot the animated curves of the GP (GIF) with the 0.95 confidence interval (optional display raw data)
#' @export
#'
#' @examples
#' TRUE
plot_animate = function(pred_gp, data = NULL, data_train = NULL, mean = NULL, mean_CI = F, file = "gganim.gif")
{
  gg = plot_gp(pred_gp, data, data_train, mean, mean_CI)  +
    gganimate::transition_states(.data$Nb_data, transition_length = 2, state_length = 1)
  gganimate::animate(gg, renderer = gganimate::gifski_renderer(file)) %>% return()
}
