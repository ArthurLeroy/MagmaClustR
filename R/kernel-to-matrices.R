#' Create covariance matrix from a kernel
#'
#' \code{kern_to_cov()} creates a covariance matrix between input values (that
#'    could be either scalars or vectors) evaluated within a kernel function,
#'    which is characterised by specified hyper-parameters. This matrix is
#'    a finite-dimensional evaluation of the infinite-dimensional covariance
#'    structure of a GP, defined thanks to this kernel.
#'
#' @param input A vector, matrix, data frame or tibble containing all inputs for
#'    one task. If a vector, the elements are used as reference, otherwise
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
#'    - convolution_kernel : the Convolution kernel used to manage MO scenario.
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

  ## Treat the convolution kernel appart from other kernels
  if (kern %>% rlang::is_function() && length(input$Output_ID %>% unique()) > 1
                                    && is.null(deriv)) {
    # Need a unique dataframe, containing all observed inputs of all outputs.
    # The convolution_kernel() function generates the whole multi-output
    # covariance matrix.

    # Call convolution_kernel() with vectorized mode to generate the whole MO
    # covariance matrix.
    mat <- convolution_kernel(x = input, y = input_2, hp = hp, vectorized = TRUE)

    # We select ONLY the columns that uniquely define a point
    # (coordinates + Output_ID) to ensure consistent naming across different
    # function calls.

    # Dynamically find all coordinate columns (starting with "Input")
    coord_cols <- grep("^Input", names(input), value = TRUE)

    # Paste the coordinate values together for each row
    pasted_coords <- do.call(paste, c(input[coord_cols], sep = ";"))

    # Prepend the output ID to create the final name
    reference <- paste0("o", input$Output_ID, ";", pasted_coords)

    # Do the same for the second set of inputs if it's different
    if (identical(input, input_2)) {
      reference_2 <- reference
    } else {
      pasted_coords_2 <- do.call(paste, c(input_2[coord_cols], sep = ";"))
      reference_2 <- paste0("o", input_2$Output_ID, ";", pasted_coords_2)
    }

    ## Add a 'noise' term on the diagonal if provided
    if (("noise" %in% names(hp)) && is.null(deriv)) {
      # Join to ensure that the noise of each point corresponds to the right
      # output
      noise_info <- input %>%
        dplyr::select(Output_ID) %>%
        dplyr::left_join(hp, by = "Output_ID")

      full_noise_vector <- exp(noise_info$noise)

      # Add the noise vector to the diagonal
      mat <- mat + diag(full_noise_vector)

    } else if (is.null(deriv)) {
      warning("Noise parameter is not provided. Data is then supposed to be ",
      "noiseless.")
    }

    return(mat %>% `rownames<-`(reference) %>% `colnames<-`(reference_2))
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
      # input <- input %>% dplyr::select(-c(Output_ID, Reference))
      input <- input %>% dplyr::select(-Reference)

      ## Format inputs to be used in a subsequent 'outer()' function
      list_input <- split(t(input),
                          rep(1:nrow(input),each = ncol(input))
      )

      ## If the Reference column exists, store it to name rows and columns
      reference_2 <- input_2$Reference %>% as.character()

      ## Only retain the actual input columns
      # input_2 <- input_2 %>% dplyr::select(-c(Output_ID, Reference))
      input_2 <- input_2 %>% dplyr::select(-Reference)

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
    if ("Output_ID" %in% names(hp)){
      if (any(startsWith(deriv, "noise"))){
        # Extract the ID of the hyper-parameter (ex: "noise_2" -> 2)
        deriv_id_str <- stringr::str_extract(deriv, "\\d+$")
        if (is.na(deriv_id_str)) {
          stop("The derivative name should contain an ID, ex: 'noise_1'")
        }
        deriv_id <- as.integer(deriv_id_str)

        if ("Task_ID" %in% colnames(input)) {
          input <- input %>% dplyr::select(-Task_ID)
        }
        input_2 <- input_2 %>% dplyr::arrange(Output_ID)
        if ("Task_ID" %in% colnames(input_2)) {
          input_2 <- input_2 %>% dplyr::select(-Task_ID)
        }

        unique_ids <- unique(input$Output_ID)
        list_of_blocks <- list()

        for (id in unique_ids) {
          subset_input <- input %>% dplyr::filter(Output_ID == id)
          subset_input_2 <- input_2 %>% dplyr::filter(Output_ID == id)

          # Compute the block only if Output_ID matches
          if (id == deriv_id) {
            # Compute the matrix derivative according to the current HP
            current_noise_hp <- hp %>%
              dplyr::filter(Output_ID == id) %>%
              dplyr::pull(noise)

            if (length(current_noise_hp) == 0) {
              stop(paste("'Noise' parameter not found for Output_ID :", id))
            }


            # Create a sub-block of the whole noise matrix
            block_matrix <- cpp_noise(
              as.matrix(dplyr::select(subset_input, Input_1)),
              as.matrix(dplyr::select(subset_input_2, Input_1)),
              current_noise_hp
            )

          } else {
            # Derivative is zero
            block_matrix <- matrix(0,
                                   nrow = nrow(subset_input),
                                   ncol = nrow(subset_input_2)) %>%
              `rownames<-` (subset_input$Input_1) %>%
              `colnames<-` (subset_input$Input_1)
          }

          list_of_blocks[[as.character(id)]] <- block_matrix
        }

        # Aggregate blocks into the complete block-diagonal matrix
        mat <- Matrix::bdiag(unname(list_of_blocks))
        mat <- as.matrix(mat) %>%
          `rownames<-`(input$Input_1) %>%
          `colnames<-` (input_2$Input_1)
        return(mat)
      }
    }

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

  # Add noise to the diagonal if provided (for standard, single-output kernels)
  if (("noise" %in% names(hp)) & is.null(deriv)) {
    # This generic noise logic is only for single-output kernels where noise is a single value.
    # The multi-output convolution_kernel has its own noise logic in the block above.
    if (length(hp[["noise"]]) == 1) {
      mat <- mat + cpp_noise(as.matrix(input), as.matrix(input_2), hp[["noise"]])
    }
  }

  ## Sanity check: detect NaN/Inf in the resulting covariance matrix
  if (any(is.nan(mat)) || any(is.infinite(mat))) {
    warning(
      "kern_to_cov: the computed covariance matrix contains NaN or Inf. ",
      "This often indicates that hyper-parameters have diverged ",
      "(e.g. extremely large variance or very small lengthscale). ",
      "Current hp: ", paste(names(hp), "=", hp, collapse = ", ")
    )
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
#'    one task. If a vector, the elements are used as reference, otherwise
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

  ## Sanity check: warn early if the covariance matrix is degenerate
  if (any(is.nan(mat_cov)) || any(is.infinite(mat_cov))) {
    warning(
      "kern_to_inv: the covariance matrix contains NaN/Inf. ",
      "This likely indicates divergent hyper-parameters."
    )
  }

  inv <- mat_cov %>%
    chol_inv_jitter(pen_diag = pen_diag) %>%
    `rownames<-`(reference) %>%
    `colnames<-`(reference) %>%
    return()
}

#' Compute a covariance matrix for multiple tasks
#'
#' Compute the covariance matrices associated with all tasks in the
#' database, taking into account their specific inputs and hyper-parameters.
#'
#' @param data A tibble or data frame of input data. Required column: 'Task_ID'.
#'   Suggested column: 'Input' (for indicating the reference input).
#' @param kern A kernel function.
#' @param hp A tibble or data frame, containing the hyper-parameters associated
#' with each task.
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
  floop <- function(t) {
    db_t <- data %>%
      dplyr::filter(Task_ID == t) %>%
      dplyr::select(-Task_ID)

    ## To avoid throwing an error if 'Output' has already been removed
    if ("Output" %in% names(db_t)) {
      db_t <- db_t %>% dplyr::select(-Output)
    }

    hp_t <- hp %>%
      dplyr::filter(Task_ID == t) %>%
      dplyr::select(-Task_ID)

    kern_to_cov(db_t, kern, hp_t, deriv = deriv) %>%
      return()
  }
  sapply(unique(data$Task_ID), floop, simplify = F, USE.NAMES = T) %>%
    return()
}

#' Compute an inverse covariance matrix for multiple tasks
#'
#' Compute the inverse covariance matrices associated with all tasks
#' in the database, taking into account their specific inputs and
#' hyper-parameters.
#'
#' @param db A tibble or data frame of input data. Required column: 'Task_ID'.
#'   Suggested column: 'Input' (for indicating the reference input).
#' @param kern A kernel function.
#' @param hp A tibble or data frame, containing the hyper-parameters associated
#' with each task.
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
  floop <- function(t) {
    db_t <- db %>%
      dplyr::filter(Task_ID == t) %>%
      dplyr::select(-Task_ID)
    ## To avoid throwing an error if 'Output' has already been removed
    if ("Output" %in% names(db_t)) {
      db_t <- db_t %>% dplyr::select(-Output)
    }

    hp_t <- hp %>%
      dplyr::filter(Task_ID == t) %>%
      dplyr::select(-Task_ID)

    kern_to_inv(db_t, kern, hp_t, pen_diag, deriv = deriv) %>%
      return()
  }
  sapply(unique(db$Task_ID), floop, simplify = F, USE.NAMES = T) %>%
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
#' @param go_one_more Logical. If TRUE and Cholesky succeeds, recurse once
#'   more with jitter multiplied by 10 (and go_one_more = FALSE) to add a
#'   safety margin. Default FALSE preserves the standard behaviour.
#' @param max_jitter A number, the maximum jitter allowed before falling
#'   back to a pseudo-inverse.
#'
#' @return A matrix, inverse of \code{mat} plus an adaptive jitter term
#'    added on the diagonal. An attribute \code{"effective_jitter"} is
#'    attached, recording the actual jitter value that was used.
#'
#' @keywords internal
#'
#' @examples
#' TRUE
chol_inv_jitter <- function(mat, pen_diag, go_one_more = FALSE,
                            max_jitter = 1e-1, inv_tol = 0.01) {

  ## Sanity check: abort early if the matrix contains NaN or Inf
  if (any(is.nan(mat)) || any(is.infinite(mat))) {
    warning(
      "The covariance matrix contains NaN or Inf values. ",
      "Returning a matrix of NaN instead of attempting inversion."
    )
    result <- matrix(NaN, nrow = nrow(mat), ncol = ncol(mat),
                     dimnames = dimnames(mat))
    attr(result, "effective_jitter") <- pen_diag
    return(result)
  }

  ## Add the jitter to the diagonal
  mat_jittered <- mat
  diag(mat_jittered) <- diag(mat_jittered) + pen_diag

  ## Attempt Cholesky decomposition and inversion
  inv <- tryCatch(
    chol(mat_jittered) %>% chol2inv(),
    error = function(e) NULL
  )

  if (!is.null(inv)) {
    ## Cholesky succeeded — verify the quality of the inverse.
    ## diag(K * K^{-1}) should be ≈ 1.  rowSums(A * B) = diag(A %*% B)
    ## for symmetric matrices, at O(n^2) cost instead of O(n^3).
    diag_product <- rowSums(mat_jittered * inv)
    max_err <- max(abs(diag_product - 1))

    if (max_err > inv_tol) {
      ## Inverse is numerically unreliable despite Cholesky succeeding.
      ## Treat as failure: escalate jitter.
      inv <- NULL
    } else if (go_one_more) {
      ## Quality OK but caller asked for an extra safety margin:
      ## redo with jitter * 10, go_one_more = FALSE.
      return(chol_inv_jitter(mat, pen_diag * 10,
                             go_one_more = FALSE,
                             max_jitter = max_jitter,
                             inv_tol = inv_tol))
    } else {
      ## Quality OK, return
      attr(inv, "effective_jitter") <- pen_diag
      return(inv)
    }
  }

  ## Cholesky failed or quality check failed — escalate jitter
  new_jitter <- pen_diag * 10

  if (new_jitter > max_jitter) {
    ## Cap reached: fall back to pseudo-inverse
    warning(
      "chol_inv_jitter: Cholesky inversion failed (jitter reached ",
      format(pen_diag, scientific = TRUE),
      "). Falling back to a pseudo-inverse (MASS::ginv). ",
      "Results may be less precise."
    )
    diag(mat) <- diag(mat) + max_jitter
    result <- MASS::ginv(mat)
    attr(result, "effective_jitter") <- max_jitter
    return(result)
  }

  ## Recurse with increased jitter (keep go_one_more as-is)
  chol_inv_jitter(mat, new_jitter, go_one_more = go_one_more,
                  max_jitter = max_jitter, inv_tol = inv_tol)
}


#' Round a matrix to make if symmetric
#'
#' If a matrix is non-symmetric due to numerical errors, round with a decreasing
#' number of digits until the matrix becomes symmetric.
#'
#' @param mat A matrix, possibly non-symmetric.
#' @param digits A number, the starting number of digits to round from if
#'    \code{mat} is not symmetric
#'
#' @return A matrix, rounded approximation of \code{mat} that is symmetric.
#'
#' @keywords internal
#'
#' @examples
#' TRUE
check_symmetric <- function(mat, digits = 10){

  if(mat %>% isSymmetric()) {
    return(mat)
  }
  else {
    ## Round matrix to remove numerical errors and make it symmetric
    mat <- round(x = mat, digits = digits)

    ## Recursive loop
    check_symmetric(mat, digits-1)
  }
}

#' @importFrom Rcpp sourceCpp
#' @useDynLib MagmaClustR, .registration = TRUE
NULL
