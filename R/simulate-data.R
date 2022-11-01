#' Simulate a batch of data
#'
#' Simulate a batch of output data, corresponding to one individual, coming from
#' a GP with a the Squared Exponential kernel as covariance structure, and
#' specified hyper-parameters and input.
#'
#' @param ID An identification code, whether numeric or character.
#' @param input A vector of numbers. The input variable that is used as
#'  'reference' for input and outputs.
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
simu_indiv_se <- function(ID, input, mean, v, l, sigma) {
  if(is.data.frame(input)){
    db <- tibble::tibble(
      "ID" = ID,
      input,
      "Output" = mvtnorm::rmvnorm(
        1,
        mean,
        kern_to_cov(
          input,
          "SE",
          tibble::tibble(se_variance = v, se_lengthscale = l)
        ) +
          diag(sigma, length(input$Reference), length(input$Reference))
      ) %>% as.vector(),
      "se_variance" = v,
      "se_lengthscale" = l,
      "noise" = sigma
    )
  }else{
    db <- tibble::tibble(
      "ID" = ID,
      "Input" = input,
      "Output" = mvtnorm::rmvnorm(
        1,
        mean,
        kern_to_cov(
          input,
          "SE",
          tibble::tibble(se_variance = v, se_lengthscale = l)
        ) +
          diag(sigma, length(input), length(input))
      ) %>% as.vector(),
      "se_variance" = v,
      "se_lengthscale" = l,
      "noise" = sigma
    )
  }
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
#' @param grid_cov A vector of numbers defining a grid of observations
#'    (i.e. the covariate reference inputs).
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
#' @param lambda_int A vector of 2 numbers, defining an interval of admissible
#'    values for the lambda parameter of the 2D exponential.
#' @param m_int A vector of 2 numbers, defining an interval of admissible
#'    values for the mean of the 2D exponential.
#' @param lengthscale_int A vector of 2 numbers, defining an interval of
#'    admissible values for the lengthscale parameter of the 2D exponential.
#' @param m0_intercept A vector of 2 numbers, defining an interval of admissible
#'    values for the intercept of m0.
#' @param m0_slope A vector of 2 numbers, defining an interval of admissible
#'    values for the slope of m0.
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
                    grid_cov = seq(0, 10, 0.5),
                    common_input = TRUE,
                    common_hp = TRUE,
                    add_hp = FALSE,
                    add_clust = FALSE,
                    int_mu_v = c(4, 5),
                    int_mu_l = c(0, 1),
                    int_i_v = c(1, 2),
                    int_i_l = c(0, 1),
                    int_i_sigma = c(0, 0.2),
                    lambda_int = c(30,40),
                    m_int = c(0,10),
                    lengthscale_int = c(30, 40),
                    m0_slope = c(-5, 5),
                    m0_intercept = c(-50, 50)) {
  if(covariate){
    exponential_mean <- function(x, y, lambda, m1, m2, lengthscale){
      return (lambda * base::exp(-((x-m1)^2 + (y-m2)^2)/lengthscale))
    }
    t_0_tibble <- tidyr::expand_grid(Input = grid, Covariate = grid_cov) %>%
        purrr::modify(signif) %>%
        tidyr::unite("Reference", sep = ":", remove = FALSE) %>%
        dplyr::arrange(.data$Reference)

    t_0 <- t_0_tibble %>% dplyr::pull(.data$Reference)

    if (common_input) {
        t_i_input <- sample(grid, N, replace = F)
        t_i_covariate <- sample(grid_cov, N, replace = F)
        t_i_tibble <- tibble::tibble(Input = t_i_input,
                                     Covariate = t_i_covariate) %>%
          purrr::modify(signif) %>%
          tidyr::unite("Reference", sep = ":", remove = FALSE) %>%
          dplyr::arrange(.data$Reference)

      t_i <- t_i_tibble %>% dplyr::pull(.data$Reference)

    }

    floop_k <- function(k) {

      m_0 <- tibble::tibble(
        Input = t_0_tibble$Input,
        Covariate = t_0_tibble$Covariate) %>%
        purrr::modify(signif) %>%
        tidyr::unite("Reference", sep = ":", remove = FALSE) %>%
        dplyr::mutate(Output =  exponential_mean(.data$Input,
                                          .data$Covariate,
                                          lambda = draw(lambda_int),
                                          m1 = draw(m_int),
                                          m2 = draw(m_int),
                                          lengthscale = draw(lengthscale_int)))

        if (common_hp) {
          i_v <- draw(int_i_v)
          i_l <- draw(int_i_l)
          i_sigma <- draw(int_i_sigma)
        }

        floop_i <- function(i, k) {
          if (!common_input) {
            t_i_input <- sample(grid, N, replace = F)
            t_i_covariate <- sample(grid_cov, N, replace = F)
            t_i_tibble <- tibble::tibble(Input = t_i_input,
                                         Covariate = t_i_covariate) %>%
              purrr::modify(signif) %>%
              tidyr::unite("Reference", sep = ":", remove = FALSE) %>%
              dplyr::arrange(.data$Reference)

            t_i <- t_i_tibble %>% dplyr::pull(.data$Reference)
          }
          if (!common_hp) {
            i_v <- draw(int_i_v)
            i_l <- draw(int_i_l)
            i_sigma <- draw(int_i_sigma)
          }
          mean_i <- m_0 %>%
            dplyr::filter(.data$Reference %in% t_i) %>%
            dplyr::pull(.data$Output)

          if (K > 1) {
            ID <- paste0("ID", i, "-Clust", k)
          } else {
            ID <- as.character(i)
          }

          simu_indiv_se(
            ID = ID,
            input = t_i_tibble,
            mean = mean_i,
            v = i_v,
            l = i_l,
            sigma = i_sigma
          ) %>%
            return()
        }

        db <- sapply(seq_len(M), floop_i, k = k,simplify = F, USE.NAMES = T) %>%
          dplyr::bind_rows() %>%
          dplyr::select(-.data$Reference)
        if (!add_hp) {
          db <- db %>%
            dplyr::select(-c("se_variance", "se_lengthscale", "noise"))
        }

        if (add_clust) {
          db <- db %>% dplyr::mutate("Cluster" = k)
        }
        return(db)
      }

      lapply(seq_len(K), floop_k) %>%
        dplyr::bind_rows() %>%
        return()
  }else{
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

        if (K > 1) {
          ID <- paste0("ID", i, "-Clust", k)
        } else {
          ID <- as.character(i)
        }

        simu_indiv_se(
          ID = ID,
          input = t_i,
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

      if (add_clust) {
        db <- db %>% dplyr::mutate("Cluster" = k)
      }
      return(db)
    }

    lapply(seq_len(K), floop_k) %>%
      dplyr::bind_rows() %>%
      return()
  }
}
