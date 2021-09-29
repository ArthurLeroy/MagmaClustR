#' foo: A package for computating the notorious bar statistic
#'
#' The foo package provides three categories of important functions:
#' foo, bar and baz.
#'
#' @section Details:
#' The foo functions ...
#'
#'@section Author(s):
#' The foo functions ...
#'
#' @section References:
#' Arthur Leroy, Pierre Latouche, Benjamin Guedj, and Servane Gey.
#' MAGMA: Inference and Prediction with Multi-Task Gaussian Processes.
#' PREPRINT arXiv:2007.10731 [cs, stat], July 2020
#'
#' Arthur Leroy, Pierre Latouche, Benjamin Guedj, and Servane Gey.
#' Cluster-Specific Predictions with Multi-Task Gaussian Processes.
#' PREPRINT arXiv:2011.07866 [cs, LG], Nov. 2020
#'
#' @section Examples:
#' #Simulation of datasets to be trained and predict to. \cr
#' data_train <- simu_db(covariate = FALSE, common_input = FALSE) \cr
#' data_pred <- simu_db(M=1, covariate = FALSE) \cr
#'
#' #Predictive distribution in Magma \cr
#' train <- train_magma(data_train) \cr
#' pred_magma(data_pred, trained_model = train)
#'
#' #Predictive distribution in Magma with cluster (MagmaClust) \cr
#' training_test = train_magma_VEM(data_train) \cr
#' pred_magma_clust(data_pred, trained_magmaclust = training_test)
#'
#' @docType package
#' @name MagmaClustR
NULL
#> NULL
