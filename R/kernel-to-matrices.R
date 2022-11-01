#' Create covariance matrix from a kernel
#'
#' \code{kern_to_cov()} creates a covariance matrix between input values (that
#'    could be either scalars or vectors) evaluated within a kernel function,
#'    which is characterised by specified hyper-parameters. This matrix is
#'    a finite-dimensional evaluation of the infinite-dimensional covariance
#'    structure of a GP, defined thanks to this kernel.
#'
#' @param input A vector, matrix, data frame or tibble containing all inputs for
#'    one individual. If a vector, the elements are used as reference, otherwise
#'    , one column should be named 'Input' to indicate that it represents the
#'    reference (e.g. 'Input' would contain the timestamps in time-series
#'    applications). The other columns are considered as being covariates. If
#'    no column is named 'Input', the first one is used by default.
#' @param kern A kernel function. Several popular kernels
#'    (see \href{https://www.cs.toronto.edu/~duvenaud/cookbook/}{The Kernel
#'    Cookbook}) are already implemented and can be selected within the
#'    following list:
#'    - "SE": (default value) the Squared Exponential Kernel (also called
#'        Radial Basis Function or Gaussian kernel),
#'    - "LIN": the Linear kernel,
#'    - "PERIO": the Periodic kernel,
#'    - "RQ": the Rational Quadratic kernel.
#'    Compound kernels can be created as sums or products of the above kernels.
#'    For combining kernels, simply provide a formula as a character string
#'    where elements are separated by whitespaces (e.g. "SE + PERIO"). As the
#'    elements are treated sequentially from the left to the right, the product
#'    operator '*' shall always be used before the '+' operators (e.g.
#'    'SE * LIN + RQ' is valid whereas 'RQ + SE * LIN' is  not).
#' @param hp A list, data frame or tibble containing the hyper-parameters used
#'    in the kernel. The name of the elements (or columns) should correspond
#'    exactly to those used in the kernel definition. If \code{hp} contains an
#'    element or a column 'Noise', its value will be added on the diagonal of
#'    the covariance matrix.
#' @param deriv A character, indicating according to which hyper-parameter the
#'    derivative should be computed. If NULL (default), the function simply
#'    returns the covariance matrix.
#' @param input_2 (optional) A vector, matrix, data frame or tibble under the
#'    same format as \code{input}. This argument should be used only when the
#'    kernel needs to be evaluated between two different sets of inputs,
#'    typically resulting in a non-square matrix.
#'
#' @export
#'
#' @return A covariance matrix, where elements are evaluations of the associated
#' kernel for each pair of reference inputs.
#'
#' @keywords internal
#'
#' @examples
#' TRUE
kern_to_cov <- function(input,
                        kern = "SE",
                        hp,
                        deriv = NULL,
                        input_2 = NULL) {
  ## If a second set of inputs is not provided, only 'input' against itself
  if (input_2 %>% is.null()) {
    input_2 <- input
  }

  ## Test whether some input values are duplicated
  if (!(unique(input) %>% identical(input))) {
    warning(
      "Some inputs are duplicated. This will result in a singular ",
      "matrix and an unexpected behaviour of the algorithm."
    )
  }

  ## Process the character string defining the covariance structure
  if (is.character(kern)) {
    kernel <- function(deriv = deriv, vectorized = TRUE, ...) {
      ## Read the character string
      str_kern <- strsplit(kern, " +")[[1]]
      ## Check whether the string correctly starts with a kernel
      if (str_kern[1] == "+" | str_kern[1] == "*") {
        stop("The string in 'kern' should start with a kernel not an operator.")
      }
      ## Initialise the first and last element of the character string
      str <- c("start", str_kern, "stop")
      ## Initialise the combined kernel and the operator
      out <- 0
      operator <- get("+")
      ## Tracking whether the kernel associated with 'deriv' has been founded
      true_deriv <- F
      past_true_deriv <- F

      for (i in 2:(length(str) - 1))
      {
        s_m <- str[i - 1]
        s <- str[i]
        s_p <- str[i + 1]
        ## Detect whether we compute kernel or its derivatives
        if (deriv %>% is.null()) {
          ## Detect the adequate kernel and combine to others
          if (s == "SE") {
            out <- operator(out, se_kernel(deriv = deriv, vectorized = T, ...))
          } else if (s == "PERIO") {
            out <- operator(
              out,
              perio_kernel(deriv = deriv, vectorized = T, ...)
            )
          } else if (s == "RQ") {
            out <- operator(out, rq_kernel(deriv = deriv, vectorized = T, ...))
          } else if (s == "LIN") {
            out <- operator(out, lin_kernel(deriv = deriv, vectorized = T, ...))
          } else if ((s == "+") | (s == "*")) {
            ## Combine the kernel with the adequate operator
            operator <- get(s)
          } else {
            stop("Incorrect character string specified in the 'kern' argument.")
          }
        } else {
          ## Check whether the element is a kernel or an operator
          if (any(s %in% c("SE", "PERIO", "RQ", "LIN"))) {
            ## Choose the correct kernel
            if (s == "SE") {
              temp_kern <- se_kernel
              if (any(deriv %in% c("se_variance", "se_lengthscale"))) {
                true_deriv <- T
                past_true_deriv <- T
              }
            } else if (s == "PERIO") {
              temp_kern <- perio_kernel
              if (any(deriv %in%
                      c("perio_variance", "perio_lengthscale", "period"))) {
                true_deriv <- T
                past_true_deriv <- T
              }
            } else if (s == "RQ") {
              temp_kern <- rq_kernel
              if (any(deriv %in%
                      c("rq_variance", "rq_lengthscale", "rq_scale"))) {
                true_deriv <- T
                past_true_deriv <- T
              }
            } else if (s == "LIN") {
              temp_kern <- lin_kernel
              if (any(deriv %in% c("lin_slope", "lin_offset"))) {
                true_deriv <- T
                past_true_deriv <- T
              }
            }

            ## Detect the appropriate situation for the derivatives calculus
            if (s_m == "start") {
              if (s_p == "stop") {
                out <- temp_kern(deriv = deriv, vectorized = T, ...)
              } else if (s_p == "+") {
                if (true_deriv) {
                  out <- operator(
                    out,
                    temp_kern(deriv = deriv, vectorized = T, ...)
                  )
                }
              } else if (s_p == "*") {
                if (true_deriv) {
                  out <- operator(
                    out,
                    temp_kern(deriv = deriv, vectorized = T, ...)
                  )
                } else {
                  out <- operator(
                    out,
                    temp_kern(deriv = NULL, vectorized = T, ...)
                  )
                }
              } else {
                stop("Incorrect character string specified in 'kern'.")
              }
            } else if (s_m == "+") {
              if (s_p == "stop") {
                if (true_deriv) {
                  out <- temp_kern(deriv = deriv, vectorized = T, ...)
                }
              } else if (s_p == "+") {
                if (true_deriv) {
                  out <- temp_kern(deriv = deriv, vectorized = T, ...)
                  break
                }
              } else if (s_p == "*") {
                stop(
                  "Incorrect character string specified in 'kern'. The ",
                  "'*' operators should come before '+' operators."
                )
              } else {
                stop("Incorrect character string specified in 'kern'.")
              }
            } else if (s_m == "*") {
              if (s_p == "stop") {
                if (true_deriv) {
                  out <- operator(
                    out,
                    temp_kern(deriv = deriv, vectorized = T, ...)
                  )
                } else {
                  out <- operator(
                    out,
                    temp_kern(deriv = NULL, vectorized = T, ...)
                  )
                }
              } else if (s_p == "+") {
                if (true_deriv) {
                  out <- operator(
                    out,
                    temp_kern(deriv = deriv, vectorized = T, ...)
                  )
                  break
                } else {
                  if (past_true_deriv) {
                    out <- operator(
                      out,
                      temp_kern(deriv = NULL, vectorized = T, ...)
                    )
                    break
                  } else {
                    out <- 0
                  }
                }
              } else if (s_p == "*") {
                if (true_deriv) {
                  out <- operator(
                    out,
                    temp_kern(deriv = deriv, vectorized = T, ...)
                  )
                } else {
                  out <- operator(
                    out,
                    temp_kern(deriv = NULL, vectorized = T, ...)
                  )
                }
              } else {
                stop("Incorrect character string specified in 'kern'.")
              }
            } else {
              stop("Incorrect character string specified in 'kern'.")
            }
          } else if ((s == "+") | (s == "*")) {
            ## Combine the kernel with the adequate operator
            operator <- get(s)
          } else {
            stop("Incorrect character string specified in the 'kern' argument.")
          }
        }
        true_deriv <- F
      }
      return(out)
    }
  } else if (is.function(kern)) {
    kernel <- kern
  } else {
    stop("Error in the 'kern' argument: please use a valid character string, or",
         "provide a valid custom kernel function"
    )
  }

  # Transform the batches of input into lists
  if (input %>% is.vector()) {
    list_input <- input
    list_input_2 <- input_2
    reference <- as.character(input)
    reference_2 <- as.character(input_2)
  } else {
    if (("Reference" %in% colnames(input)) &
        ("Reference" %in% colnames(input_2))) {
      ## If the Reference column exists, store it to name rows and columns
      reference <- input$Reference %>% as.character()

      ## Only retain the actual input columns
      input <- input %>% dplyr::select(-.data$Reference)

      ## Format inputs to be used in a subsequent 'outer()' function
      list_input <- split(t(input),
                          rep(1:nrow(input),each = ncol(input))
      )

      ## If the Reference column exists, store it to name rows and columns
      reference_2 <- input_2$Reference %>% as.character()

      ## Only retain the actual input columns
      input_2 <- input_2 %>% dplyr::select(-.data$Reference)

      ## Format inputs to be used in a subsequent 'outer()' function
      list_input_2 <- split(t(input_2),
                            rep(1:nrow(input_2), each = ncol(input_2))
      )
    } else {
      ## Create the Reference column if absent
      reference <- tidyr::unite(
        as.data.frame(input),
        'Reference',
        tidyselect::everything(),
        sep = ':') %>%
        dplyr::pull(.data$Reference) %>%
        as.character()

      reference_2 <- tidyr::unite(
        as.data.frame(input_2),
        'Reference',
        tidyselect::everything(),
        sep = ':') %>%
        dplyr::pull(.data$Reference) %>%
        as.character()

      ## Format inputs to be used in a subsequent 'outer()' function
      list_input <- split(t(input),
                          rep(1:nrow(input), each = ncol(input))
      )
      list_input_2 <- split(t(input_2),
                            rep(1:nrow(input_2), each = ncol(input_2))
      )
    }
  }

  ## Return the derivative of the noise if required
  if (!is.null(deriv)) {
    if (deriv == "noise") {
      mat <- cpp_noise(as.matrix(input),
                       as.matrix(input_2),
                       hp[["noise"]]
      ) %>%
        `rownames<-`(reference) %>%
        `colnames<-`(reference_2)
      return(mat)
    }
  }

  ## Detect whether the custom kernel provides derivatives
  if ("deriv" %in% methods::formalArgs(kernel)) {
    ## Detect whether speed-up vectorised computation is provided
    if ("vectorized" %in% methods::formalArgs(kernel)) {
      mat <- kernel(
        x = input,
        y = input_2,
        hp = hp,
        deriv = deriv,
        vectorized = T
      )
    } else { ## Compute the matrix element by element
      mat <- outer(
        list_input, list_input_2,
        Vectorize(function(x, y) kernel(x, y, hp, deriv = deriv))
      )
    }
  } else {
    ## Detect whether speed-up vectorised computation is provided
    if ("vectorized" %in% methods::formalArgs(kernel)) {
      mat <- kernel(x = input, y = input_2, hp = hp, vectorized = TRUE)
    } else { ## Compute the matrix element by element
      mat <- outer(
        list_input, list_input_2,
        Vectorize(function(x, y) kernel(x, y, hp))
      )
    }
  }

  ## Add noise on the diagonal if provided (and if not computing gradients)
  if (("noise" %in% names(hp)) & is.null(deriv)) {
    mat <- mat + cpp_noise(as.matrix(input), as.matrix(input_2), hp[["noise"]])
  }

  mat %>%
    `rownames<-`(reference) %>%
    `colnames<-`(reference_2) %>%
    return()
}

#' Create inverse of a covariance matrix from a kernel
#'
#' \code{kern_to_inv()} creates the inverse of a covariance matrix between
#'    input values (that could be either scalars or vectors) evaluated within
#'    a kernel function, which is characterised by specified hyper-parameters.
#'    This matrix is a finite-dimensional evaluation of the
#'    infinite-dimensional covariance structure of a GP, defined thanks to this
#'    kernel.
#'
#' @param input A vector, matrix, data frame or tibble containing all inputs for
#'    one individual. If a vector, the elements are used as reference, otherwise
#'    ,one column should be named 'Input' to indicate that it represents the
#'    reference (e.g. 'Input' would contain the timestamps in time-series
#'    applications). The other columns are considered as being covariates. If
#'    no column is named 'Input', the first one is used by default.
#' @param kern A kernel function. Several popular kernels
#'    (see \href{https://www.cs.toronto.edu/~duvenaud/cookbook/}{The Kernel
#'    Cookbook}) are already implemented and can be selected within the
#'    following list:
#'    - "SE": (default value) the Squared Exponential Kernel (also called
#'        Radial Basis Function or Gaussian kernel),
#'    - "LIN": the Linear kernel,
#'    - "PERIO": the Periodic kernel,
#'    - "RQ": the Rational Quadratic kernel.
#'    Compound kernels can be created as sums or products of the above kernels.
#'    For combining kernels, simply provide a formula as a character string
#'    where elements are separated by whitespaces (e.g. "SE + PERIO"). As the
#'    elements are treated sequentially from the left to the right, the product
#'    operator '*' shall always be used before the '+' operators (e.g.
#'    'SE * LIN + RQ' is valid whereas 'RQ + SE * LIN' is  not).
#' @param hp A list, data frame or tibble containing the hyper-parameters used
#'   in the kernel. The name of the elements (or columns) should correspond
#'   exactly to those used in the kernel definition.
#' @param pen_diag A jitter term that is added to the covariance matrix to avoid
#'    numerical issues when inverting, in cases of nearly singular matrices.
#' @param deriv A character, indicating according to which hyper-parameter the
#'  derivative should be computed. If NULL (default), the function simply returns
#'  the inverse covariance matrix.
#'
#' @return The inverse of a covariance matrix, which elements are evaluations of
#'    the associated kernel for each pair of reference inputs.
#'
#' @export
#'
#' @keywords internal
#'
#' @examples
#' TRUE
kern_to_inv <- function(input, kern, hp, pen_diag = 1e-10, deriv = NULL) {

  mat_cov <- kern_to_cov(input = input, kern = kern, hp = hp, deriv = deriv)
  reference <- row.names(mat_cov)

  inv <- mat_cov %>%
    chol_inv_jitter(pen_diag = pen_diag) %>%
    `rownames<-`(reference) %>%
    `colnames<-`(reference) %>%
    return()
}

#' Compute a covariance matrix for multiple individuals
#'
#' Compute the covariance matrices associated with all individuals in the
#' database, taking into account their specific inputs and hyper-parameters.
#'
#' @param data A tibble or data frame of input data. Required column: 'ID'.
#'   Suggested column: 'Input' (for indicating the reference input).
#' @param kern A kernel function.
#' @param hp A tibble or data frame, containing the hyper-parameters associated
#' with each individual.
#' @param deriv A character, indicating according to which hyper-parameter the
#'  derivative should be computed. If NULL (default), the function simply returns
#'  the list of covariance matrices.
#'
#' @return A named list containing all of the inverse covariance matrices.
#'
#' @export
#'
#' @keywords internal
#'
#' @examples
#' TRUE
list_kern_to_cov <- function(data, kern, hp, deriv = NULL) {
  floop <- function(i) {
    db_i <- data %>%
      dplyr::filter(.data$ID == i) %>%
      dplyr::select(-.data$ID)

    ## To avoid throwing an error if 'Output' has already been removed
    if ("Output" %in% names(db_i)) {
      db_i <- db_i %>% dplyr::select(-.data$Output)
    }

    hp_i <- hp %>%
      dplyr::filter(.data$ID == i) %>%
      dplyr::select(-.data$ID)

    kern_to_cov(db_i, kern, hp_i, deriv = deriv) %>%
      return()
  }
  sapply(unique(data$ID), floop, simplify = F, USE.NAMES = T) %>%
    return()
}

#' Compute an inverse covariance matrix for multiple individuals
#'
#' Compute the inverse covariance matrices associated with all individuals
#' in the database, taking into account their specific inputs and
#' hyper-parameters.
#'
#' @param db A tibble or data frame of input data. Required column: 'ID'.
#'   Suggested column: 'Input' (for indicating the reference input).
#' @param kern A kernel function.
#' @param hp A tibble or data frame, containing the hyper-parameters associated
#' with each individual.
#' @param pen_diag A number. A jitter term, added on the diagonal to prevent
#' numerical issues when inverting nearly singular matrices.
#' @param deriv A character, indicating according to which hyper-parameter the
#'  derivative should be computed. If NULL (default), the function simply returns
#'  the list of covariance matrices.
#'
#' @return A named list containing all of the inverse covariance matrices.
#'
#' @export
#'
#' @keywords internal
#'
#' @examples
#' TRUE
list_kern_to_inv <- function(db, kern, hp, pen_diag, deriv = NULL) {
  floop <- function(i) {
    db_i <- db %>%
      dplyr::filter(.data$ID == i) %>%
      dplyr::select(-.data$ID)
    ## To avoid throwing an error if 'Output' has already been removed
    if ("Output" %in% names(db_i)) {
      db_i <- db_i %>% dplyr::select(-.data$Output)
    }

    hp_i <- hp %>%
      dplyr::filter(.data$ID == i) %>%
      dplyr::select(-.data$ID)

    kern_to_inv(db_i, kern, hp_i, pen_diag, deriv = deriv) %>%
      return()
  }
  sapply(unique(db$ID), floop, simplify = F, USE.NAMES = T) %>%
    return()
}

#' Inverse a matrix using an adaptive jitter term
#'
#' Inverse a matrix from its Choleski decomposition. If (nearly-)singular,
#' increase the order of magnitude of the jitter term added to the diagonal
#' until the matrix becomes non-singular.
#'
#' @param mat A matrix, possibly singular.
#' @param pen_diag A number, a jitter term to add on the diagonal.
#'
#' @return A matrix, inverse of \code{mat} plus an adaptive jitter term
#'    added on the diagonal.
#'
#' @keywords internal
#'
#' @examples
#' TRUE
chol_inv_jitter <- function(mat, pen_diag){
  ## Add a jitter term to the diagonal
  diag(mat) <- diag(mat) + pen_diag
  ## Recursive pattern for the adaptive jitter (if error, increase jitter)
  tryCatch(
    mat %>% chol() %>% chol2inv(),
    error = function(e) {
      chol_inv_jitter(mat, 10*pen_diag)
      }
    )
}

#' @importFrom Rcpp sourceCpp
#' @useDynLib MagmaClustR, .registration = TRUE
NULL
