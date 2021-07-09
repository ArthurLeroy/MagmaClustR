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
#' @param kern A kernel function, associated with the GP.
#' @param v A number. The variance hyper-parameter of the SE kernel.
#' @param l A number. The lengthscale hyper-parameter of the SE kernel.
#' @param sigma A number. The noise hyper-parameter.
#'
#' @return A tibble containing a batch of output data along with input and
#'  additional information for a simulated individual.
#'
#' @importFrom rlang .data
#'
#' @examples
#' simu_indiv_se("A", 1:10, 0, rep(0, 10), "SE", 2, 1, 0.5)
#' simu_indiv_se("B", 1:10, 2:11, 3:12, "SE", 1, 1, 1)
#' simu_indiv_se("C", 1:10, 5, rep(0, 10), "SE", 2, 1, 0.5)
simu_indiv_se <- function(ID, input, covariate, mean, kern, v, l, sigma) {
  db <- tibble::tibble(
    "ID" = ID,
    "Output" = mvtnorm::rmvnorm(
      1,
      mean + covariate,
      kern_to_cov(input, kern, tibble::tibble(variance = v, lengthscale = l)) +
        diag(sigma, length(input), length(input))
    ) %>% as.vector(),
    "Input" = input,
    "Covariate" = covariate,
    "variance" = v,
    "lengthscale" = l,
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
#' @examples
#' draw(c(1, 2))
draw <- function(int){
  stats::runif(1, int[1], int[2]) %>%
    round(2) %>%
    return()
}

#' Simulate a complete training dataset
#'
#' Simulate a complete training dataset, which may be coherent in many different
#' settings. The various flexible arguments allow to adjust the number of
#' individuals, of observed input, and the values of many parameters
#' controlling the data generation.
#'
#' @param M An integer. The number of individual.
#' @param N An integer. The number of observations per individual.
#' @param grid A vector of numbers defining a grid of observations
#' (i.e. the reference inputs).
#' @param common_input A logical value indicating whether the reference inputs are
#'  common to all individual.
#' @param common_hp  A logical value indicating whether the hyper-parameters are
#'  common to all individual.
#' @param add_hp A logical value indicating whether the values of
#' hyper-parameters should be added as columns of the dataset.
#' @param kern_0 A kernel function, associated with the mean process.
#' @param kern_i A kernel function, associated with the individual processes.
#' @param int_mu_v A vector of 2 numbers, defining an interval of admissible
#' values for the variance hyper-parameter of the mean process' kernel.
#' @param int_mu_l A vector of 2 numbers, defining an interval of admissible
#' values for the lengthscale hyper-parameter of the mean process' kernel.
#' @param int_i_v A vector of 2 numbers, defining an interval of admissible
#' values for the variance hyper-parameter of the individual process' kernel.
#' @param int_i_l A vector of 2 numbers, defining an interval of admissible
#' values for the lengthscale hyper-parameter of the individual process' kernel.
#' @param int_i_sigma A vector of 2 numbers, defining an interval of admissible
#' values for the noise hyper-parameter.
#' @param m0_intercept A vector of 2 numbers, defining an interval of admissible
#' values for the intercept of m0.
#' @param m0_slope A vector of 2 numbers, defining an interval of admissible
#' values for the slope of m0.
#' @param covariate A vector of 2 numbers, defining an interval of admissible
#' values for the covariate inputs.
#'
#' @return A full dataset of simulated training data
#' @export
#'
#' @examples
#' simu_db(M = 5, N = 3)
#' simu_db(M = 5, N = 3, common_input = F)
#' simu_db(M = 5, N = 3, common_hp = F, add_hp = T)
#' simu_db(M = 5, N = 3, common_input = F, common_hp = F)
simu_db <- function(M = 10,
                    N = 10,
                    grid = seq(0, 10, 0.05),
                    common_input = T,
                    common_hp = T,
                    add_hp = F,
                    kern_0 = "SE",
                    kern_i = "SE",
                    int_mu_v = c(0, 5),
                    int_mu_l = c(0, 2),
                    int_i_v = c(0, 5),
                    int_i_l = c(0, 2),
                    int_i_sigma = c(0, 1),
                    m0_slope =  c(-2, 2),
                    m0_intercept = c(0, 10),
                    covariate = c(-5, 5)) {
  m_0 <- draw(m0_intercept) + draw(m0_slope) * grid
  mu_v <- draw(int_mu_v)
  mu_l <- draw(int_mu_l)

  db_0 <- simu_indiv_se(
    ID = "0", input = grid, covariate = 0, mean = m_0,
    kern = kern_0, v = mu_v, l = mu_l, sigma = 0
  )
  if (common_input) {
    t_i <- sample(grid, N, replace = F) %>% sort()
  }
  if (common_hp) {
    i_v <- draw(int_i_v)
    i_l <- draw(int_i_l)
    i_sigma <- draw(int_i_sigma)
  }

  floop <- function(i) {
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

    covariate_i <- stats::runif(N, covariate[1], covariate[2]) %>% round(2)

    simu_indiv_se(
      as.character(i),
      t_i,
      covariate_i,
      mean_i,
      kern_i,
      i_v,
      i_l,
      i_sigma
    ) %>%
      return()
  }
  db = sapply(seq_len(M), floop, simplify = F, USE.NAMES = T) %>%
    dplyr::bind_rows()
  if(!add_hp){db = db %>% dplyr::select(- c('variance', 'lengthscale', 'noise'))}

  return(db)
}
