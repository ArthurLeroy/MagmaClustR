#' Plot smoothed raw data
#'
#' @param db database with format : ID, Input, Output
#' @param cluster boolean indicating whether you want cluster or not
#' @param legend boolean indicating whether want legend or not
#'
#' @return
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
