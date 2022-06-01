#' MagmaClustR : Clustering and Prediction using Multi-Task Gaussian Processes
#'
#' The \strong{MagmaClustR} package implements two main algorithms, called
#'  \emph{Magma} and \emph{MagmaClust}, using a multi-task GPs model to perform
#'   predictions for supervised learning problems. Theses approaches leverage
#'   the learning of cluster-specific mean processes, which are common across
#'   similar tasks, to provide enhanced prediction performances (even far from
#'   data) at a linear computational cost (in the number of tasks).
#'   \emph{MagmaClust} is a generalisation of \emph{Magma} where the tasks are
#'   simultaneously clustered into groups, each being associated to a specific
#'   mean process. User-oriented functions in the package are decomposed into
#'   training, prediction and plotting functions. Some basic features of
#'   standard GPs are also implemented.
#'
#' @section Details:
#' For a quick introduction to \pkg{MagmaClustR}, please refer to the README at
#' \url{https://github.com/ArthurLeroy/MagmaClustR}
#'
#' @section Author(s):
#' Arthur Leroy, Pierre Pathe and Pierre Latouche \cr
#' Maintainer: Arthur Leroy - \email{arthur.leroy.pro@@gmail.com}
#'
#' @section References:
#' Arthur Leroy, Pierre Latouche, Benjamin Guedj, and Servane Gey. \cr
#' MAGMA: Inference and Prediction with Multi-Task Gaussian Processes.
#' *Machine Learning*, 2022,
#' \url{https://link.springer.com/article/10.1007/s10994-022-06172-1}
#'
#' Arthur Leroy, Pierre Latouche, Benjamin Guedj, and Servane Gey. \cr
#' Cluster-Specific Predictions with Multi-Task Gaussian Processes.
#' *PREPRINT*, Nov. 2020, \url{https://arxiv.org/abs/2011.07866}
#'
#' @section Examples:
#'
#' ### Simulate a dataset, train and predict with Magma \cr
#' set.seed(42) \cr
#' data_magma <- simu_db(M = 11, N = 10, K = 1) \cr
#' magma_train <- data_magma %>% subset(ID %in% 1:10) \cr
#' magma_test <- data_magma %>% subset(ID == 11) %>% head(5) \cr
#'
#' magma_model <- train_magma(data = magma_train) \cr
#' magma_pred  <- pred_magma(data = magma_test, trained_model = magma_model,
#' grid_inputs = seq(0, 10, 0.01)) \cr
#'
#' ### Simulate a dataset, train and predict with MagmaClust \cr
#' set.seed(42) \cr
#' data_magmaclust <- simu_db(M = 4, N = 10, K = 3) \cr
#' list_ID = unique(data_magmaclust$ID) \cr
#' magmaclust_train <- data_magmaclust %>% subset(ID %in% list_ID\[1:11\]) \cr
#' magmaclust_test <- data_magmaclust %>% subset(ID == list_ID\[12\]) %>%
#'  head(5)\cr
#'
#' magmaclust_model <- train_magmaclust(data = magmaclust_train) \cr
#' magmaclust_pred  <- pred_magmaclust(data = magmaclust_test, \cr
#'   trained_model = magmaclust_model, grid_inputs = seq(0, 10, 0.01)) \cr
#'
#' @docType package
#' @name MagmaClustR
NULL
#> NULL
