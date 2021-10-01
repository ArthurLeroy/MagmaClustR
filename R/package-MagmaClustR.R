#' MagmaClustR : A package for Modeling, predicting and clustering
#'
#' The MagmaClustR package provides two categories of important functions:
#' training and prediction. Each part divide in whether or not we want to
#' do a clustering.
#'
#' @section Details:
#' For a quick introduction to \pkg{MagmaClustR} see the vignette
#' \href{../doc/Introduction_MagmaClustR.html}{Introduction of MagmaClustR}. \cr
#' Or run the vignette with the code : \cr
#' \code{vignette("Introduction_MagmaClustR", package = "MagmaClustR")} \cr
#'
#' For a deepening to the functions using by \pkg{MagmaClustR} see the vignette
#' \href{../doc/Details.html}{Introduction of MagmaClustR}. \cr
#' Or run the vignette with the code : \cr
#' \code{vignette("Details", package = "MagmaClustR")} \cr
#'
#'@section Author(s):
#' Arthur Leroy, Pierre Latouche and Pierre Pathe \cr
#' Maintainer: Arthur Leroy \email{arthur.leroy.pro@@gmail.com}
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
#' set.seed(5) \cr
#'
#' #Simulation of datasets to be trained and predict to. \cr
#' data_train <- simu_db(covariate = FALSE) \cr
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
