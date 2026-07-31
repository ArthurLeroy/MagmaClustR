## Finite-difference checks of the multi-output ELBO gradients.

.eb <- local({
  set.seed(314)
  outs <- c("A", "B"); clusters <- c("K1", "K2")
  inp <- tibble::tibble(Output_ID = rep(outs, each = 3),
                        Input = rep(c(0.2, 0.5, 0.9), 2))
  inp$Reference <- paste(inp$Output_ID, inp$Input, sep = ":")
  n <- nrow(inp); refs <- inp$Reference
  mk_cov <- function() {
    A <- matrix(stats::rnorm(n * n), n)
    C <- crossprod(A) / n + diag(1e-2, n)
    dimnames(C) <- list(refs, refs); C
  }
  hyperpost <- list(
    mean = stats::setNames(
      lapply(clusters, function(k) dplyr::mutate(inp, Output = stats::rnorm(n))),
      clusters),
    cov = stats::setNames(lapply(clusters, function(k) mk_cov()), clusters),
    mixture = tibble::tibble(ID = c("i1", "i2"), K1 = c(0.6, 0.3), K2 = c(0.4, 0.7))
  )
  db <- dplyr::bind_rows(
    dplyr::mutate(inp, ID = "i1", Output = stats::rnorm(n)),
    dplyr::mutate(inp, ID = "i2", Output = stats::rnorm(n))
  )
  hp_tbl <- tibble::tibble(Output_ID = outs, l_t = c(-0.3, 0.2), S_t = c(0.1, -0.2),
                           l_u_t = c(0.25, 0.25), noise = c(-3, -4))
  list(hyperpost = hyperpost, db = db, flat = flatten_hp_mo(hp_tbl),
       mean_k = stats::setNames(lapply(clusters, function(k) rep(0, n)), clusters))
})

.fd_grad <- function(fn, flat, eps = 1e-6) {
  vapply(names(flat), function(p) {
    hp_p <- flat; hp_m <- flat
    hp_p[[p]] <- hp_p[[p]] + eps; hp_m[[p]] <- hp_m[[p]] - eps
    (fn(hp_p) - fn(hp_m)) / (2 * eps)
  }, numeric(1))
}

test_that("gr_clust_multi_GP matches finite differences (multi-output)", {
  db_i <- dplyr::filter(.eb$db, ID == "i1")
  ana <- gr_clust_multi_GP(.eb$flat, db_i, .eb$hyperpost, convolution_kernel, 1e-10)
  fd <- .fd_grad(function(h)
    elbo_clust_multi_GP(h, db_i, .eb$hyperpost, convolution_kernel, 1e-10), .eb$flat)
  expect_equal(as.numeric(ana), as.numeric(fd[names(ana)]), tolerance = 1e-3)
})

test_that("gr_clust_multi_GP_shared_hp_i matches finite differences (multi-output)", {
  ana <- gr_clust_multi_GP_shared_hp_i(.eb$flat, .eb$db, .eb$hyperpost,
                                       convolution_kernel, 1e-10)
  fd <- .fd_grad(function(h)
    elbo_clust_multi_GP_shared_hp_i(h, .eb$db, .eb$hyperpost, convolution_kernel, 1e-10),
    .eb$flat)
  expect_equal(as.numeric(ana), as.numeric(fd[names(ana)]), tolerance = 1e-3)
})

test_that("gr_GP_mod_shared_hp_k matches finite differences (multi-output)", {
  ana <- gr_GP_mod_shared_hp_k(.eb$flat, .eb$hyperpost$mean, .eb$mean_k,
                               convolution_kernel, .eb$hyperpost$cov, 1e-10)
  fd <- .fd_grad(function(h)
    elbo_GP_mod_shared_hp_k(h, .eb$hyperpost$mean, .eb$mean_k, convolution_kernel,
                            .eb$hyperpost$cov, 1e-10), .eb$flat)
  expect_equal(as.numeric(ana), as.numeric(fd[names(ana)]), tolerance = 1e-3)
})
