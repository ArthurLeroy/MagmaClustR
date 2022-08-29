#' Simulate a batch a data
#'
#' Simulate a batch a output data, corresponding to one individual, coming from
#' a GP with a the Squared Exponential kernel as covariance structure, and
#' specified hyper-parameters and input.
#'
#' @param ID An identification code, whether numeric or character.
#' @param input A vector of numbers. The input variable that is used as
#'  'reference' for input and outputs.
#' @param covariate A vector of numbers. An additional input variable, observed
#'  along with each reference input.
#' @param mean A vector of numbers. Prior mean values of the GP.
#' @param v A number. The variance hyper-parameter of the SE kernel.
#' @param l A number. The lengthscale hyper-parameter of the SE kernel.
#' @param sigma A number. The noise hyper-parameter.
#'
#' @return A tibble containing a batch of output data along with input and
#'  additional information for a simulated individual.
#'
#' @importFrom rlang .data
#'
#' @keywords internal
#'
#' @examples
#' TRUE
simu_indiv_se <- function(ID, input, covariate, mean, v, l, sigma) {
  db <- tibble::tibble(
    "ID" = ID,
    "Output" = mvtnorm::rmvnorm(
      1,
      mean + covariate,
      kern_to_cov(
        input,
        "SE",
        tibble::tibble(se_variance = v, se_lengthscale = l)
      ) +
        diag(sigma, length(input), length(input))
    ) %>% as.vector(),
    "Input" = input,
    "Covariate" = covariate,
    "se_variance" = v,
    "se_lengthscale" = l,
    "noise" = sigma
  )
  return(db)
}


#' Draw a number
#'
#' Draw uniformly a number within a specified interval
#'
#' @param int An interval of values we want to draw uniformly in.
#'
#' @return A 2-decimals-rounded random number
#'
#' @keywords internal
#'
#' @examples
#' TRUE
draw <- function(int) {
  stats::runif(1, int[1], int[2]) %>%
    round(2) %>%
    return()
}

#' Simulate a dataset tailored for MagmaClustR
#'
#' Simulate a complete training dataset, which may be representative of various
#' applications. Several flexible arguments allow adjustment of the number of
#' individuals, of observed inputs, and the values of many parameters
#' controlling the data generation.
#'
#' @param M An integer. The number of individual per cluster.
#' @param N An integer. The number of observations per individual.
#' @param K An integer. The number of underlying clusters.
#' @param covariate A logical value indicating whether the dataset should
#'    include an additional input covariate named 'Covariate'.
#' @param grid A vector of numbers defining a grid of observations
#'    (i.e. the reference inputs).
#' @param common_input A logical value indicating whether the reference inputs
#'    are common to all individual.
#' @param common_hp  A logical value indicating whether the hyper-parameters are
#'   common to all individual. If TRUE and K>1, the hyper-parameters remain
#'   different between the clusters.
#' @param add_hp A logical value indicating whether the values of
#'    hyper-parameters should be added as columns in the dataset.
#' @param add_clust A logical value indicating whether the name of the
#'    clusters should be added as a column in the dataset.
#' @param int_mu_v A vector of 2 numbers, defining an interval of admissible
#'    values for the variance hyper-parameter of the mean process' kernel.
#' @param int_mu_l A vector of 2 numbers, defining an interval of admissible
#'    values for the lengthscale hyper-parameter of the mean process' kernel.
#' @param int_i_v A vector of 2 numbers, defining an interval of admissible
#'    values for the variance hyper-parameter of the individual process' kernel.
#' @param int_i_l A vector of 2 numbers, defining an interval of admissible
#'    values for the lengthscale hyper-parameter of the individual process'
#'    kernel.
#' @param int_i_sigma A vector of 2 numbers, defining an interval of admissible
#'    values for the noise hyper-parameter.
#' @param m0_intercept A vector of 2 numbers, defining an interval of admissible
#'    values for the intercept of m0.
#' @param m0_slope A vector of 2 numbers, defining an interval of admissible
#'    values for the slope of m0.
#' @param int_covariate A vector of 2 numbers, defining an interval of
#'    admissible values for the covariate inputs.
#'
#' @return A full dataset of simulated training data.
#' @export
#'
#' @examples
#' ## Generate a dataset with 3 clusters of 4 individuals, observed at 10 inputs
#' data = simu_db(M = 4, N = 10, K = 3)
#'
#' ## Generate a 2-D dataset with an additional input 'Covariate'
#' data = simu_db(covariate = TRUE)
#'
#' ## Generate a dataset where input locations are different among individuals
#' data = simu_db(common_input = FALSE)
#'
#' ## Generate a dataset with an additional column indicating the true clusters
#' data = simu_db(K = 3, add_clust = TRUE)
simu_db <- function(M = 10,
                    N = 10,
                    K = 1,
                    covariate = FALSE,
                    grid = seq(0, 10, 0.05),
                    common_input = TRUE,
                    common_hp = TRUE,
                    add_hp = FALSE,
                    add_clust = FALSE,
                    int_mu_v = c(4, 5),
                    int_mu_l = c(0, 1),
                    int_i_v = c(1, 2),
                    int_i_l = c(0, 1),
                    int_i_sigma = c(0, 0.2),
                    m0_slope = c(-5, 5),
                    m0_intercept = c(-50, 50),
                    int_covariate = c(-5, 5)) {
  if (common_input) {
    t_i <- sample(grid, N, replace = F) %>% sort()
  }

  floop_k <- function(k) {
    m_0 <- draw(m0_intercept) + draw(m0_slope) * grid
    mu_v <- draw(int_mu_v)
    mu_l <- draw(int_mu_l)

    db_0 <- simu_indiv_se(
      ID = "0",
      input = grid,
      covariate = 0,
      mean = m_0,
      v = mu_v,
      l = mu_l,
      sigma = 0
    )
    if (common_hp) {
      i_v <- draw(int_i_v)
      i_l <- draw(int_i_l)
      i_sigma <- draw(int_i_sigma)
    }

    floop_i <- function(i, k) {
      if (!common_input) {
        t_i <- sample(grid, N, replace = F) %>% sort()
      }
      if (!common_hp) {
        i_v <- draw(int_i_v)
        i_l <- draw(int_i_l)
        i_sigma <- draw(int_i_sigma)
      }
      mean_i <- db_0 %>%
        dplyr::filter(.data$Input %in% t_i) %>%
        dplyr::pull(.data$Output)

      covariate_i <- stats::runif(N, int_covariate[1], int_covariate[2]) %>%
        round(2)

      if (K > 1) {
        ID <- paste0("ID", i, "-Clust", k)
      } else {
        ID <- as.character(i)
      }

      simu_indiv_se(
        ID = ID,
        input = t_i,
        covariate = covariate_i,
        mean = mean_i,
        v = i_v,
        l = i_l,
        sigma = i_sigma
      ) %>%
        return()
    }
    db <- sapply(seq_len(M), floop_i, k = k, simplify = F, USE.NAMES = T) %>%
      dplyr::bind_rows()
    if (!add_hp) {
      db <- db %>% dplyr::select(-c("se_variance", "se_lengthscale", "noise"))
    }

    if (!covariate) {
      db <- db %>% dplyr::select(-.data$Covariate)
    }

    if (add_clust) {
      db <- db %>% dplyr::mutate("Cluster" = k)
    }
    return(db)
  }

  lapply(seq_len(K), floop_k) %>%
    dplyr::bind_rows() %>%
    return()
}

#' Run a k-means algoithm to initialise clusters' allocation
#'
#' @param data A tibble containing common Input and associated Output values
#'   to cluster.
#' @param k A number of clusters assumed for running the kmeans algorithm.
#' @param nstart A number, indicating how many re-starts of kmeans are set.
#' @param summary A boolean, indicating whether we want an outcome summary
#'
#' @return A tibble containing the initial clustering obtained through kmeans.
#'
#' @keywords internal
#'
#' @examples
#' TRUE
ini_kmeans <- function(data, k, nstart = 50, summary = FALSE) {
  # if (!identical(
  #   unique(data$Input),
  #   data %>%
  #     dplyr::filter(.data$ID == unique(data$ID)[[1]]) %>%
  #     dplyr::pull(.data$Input)
  # )) {
    floop <- function(i) {
      obs_i <- data %>%
        dplyr::filter(.data$ID == i) %>%
        dplyr::pull(.data$Output)
      tibble::tibble(
        "ID" = i,
        "Input" = seq_len(3),
        "Output" = c(min(obs_i), mean(obs_i), max(obs_i))
      ) %>%
        return()
    }
    db_regular <- unique(data$ID) %>%
      lapply(floop) %>%
      dplyr::bind_rows() %>%
      dplyr::select(c(.data$ID, .data$Input, .data$Output))
  # } else {
  # db_regular <- data %>% dplyr::select(c(.data$ID, .data$Input, .data$Output))
  # }

  res <- db_regular %>%
    tidyr::spread(key = .data$Input, value = .data$Output) %>%
    dplyr::select(-.data$ID) %>%
    stats::kmeans(centers = k, nstart = nstart)

  if (summary) {
    res %>% print()
  }

  broom::augment(
    res,
    db_regular %>% tidyr::spread(key = .data$Input, value = .data$Output)
  ) %>%
    dplyr::select(c(.data$ID, .data$.cluster)) %>%
    dplyr::rename(Cluster_ini = .data$.cluster) %>%
    dplyr::mutate(Cluster_ini = paste0("K", .data$Cluster_ini)) %>%
    return()
}


#' Mixture initialisation with kmeans
#'
#' Provide an initial kmeans allocation of the individuals/tasks in a dataset
#' into a definite number of clusters, and return the associated mixture
#' probabilities.
#'
#' @param data A tibble or data frame. Required columns: \code{ID}, \code{Input}
#'    , \code{Output}.
#' @param k A number, indicating the number of clusters.
#' @param name_clust A vector of characters. Each element should correspond to
#'    the name of one cluster.
#' @param nstart A number of restart used in the underlying kmeans algorithm
#'
#' @return A tibble indicating for each \code{ID} in which cluster it belongs
#'    after a kmeans initialisation.
#'
#' @keywords internal
#'
#' @examples
#' TRUE
ini_mixture <- function(data, k, name_clust = NULL, nstart = 50) {
  db_ini <- ini_kmeans(data, k, nstart) %>%
    dplyr::mutate(value = 1) %>%
    tidyr::spread(key = .data$Cluster_ini, value = .data$value, fill = 0)

  if (!is.null(name_clust)) {
    names(db_ini) <- c("ID", name_clust)
  }

  return(db_ini)
}
