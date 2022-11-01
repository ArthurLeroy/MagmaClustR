test_that("kern_to_cov() works for scalar inputs", {
  hp <- tibble::tibble(se_variance = 2, se_lengthscale = 1)
  input <- c(2, 3, 4)

  res <- matrix(NA, ncol = 3, nrow = 3)
  for (i in 1:3)
  {
    for (j in 1:3) {
      res[i, j] <- se_kernel(input[i], input[j], hp)
    }
  }
  res <- res %>%
    `rownames<-`(as.character(input)) %>%
    `colnames<-`(as.character(input))

  kern_to_cov(input, "SE", hp) %>% expect_equal(res)
})

test_that("kern_to_cov() works for vector inputs", {
  hp <- tibble::tibble(se_variance = 2, se_lengthscale = 1)
  input <- data.frame(Input = c(1, 2, 3), Cov1 = c(2, 3, 4))

  res <- matrix(NA, ncol = 3, nrow = 3)
  for (i in 1:3)
  {
    for (j in 1:3) {
      res[i, j] <- se_kernel(input[i, ], input[j, ], hp)
    }
  }

  ref = paste(input$Input, Cov1 = input$Cov1, sep=":")
  res <- res %>%
    `rownames<-`(ref) %>%
    `colnames<-`(ref)

  kern_to_cov(input, "SE", hp) %>% expect_equal(res)
})

test_that("matrix, dataframe and tibble work the same", {
  hp <- tibble::tibble(se_variance = 2, se_lengthscale = 1)
  df <- data.frame(Input = c(1, 2, 3), Cov1 = c(2, 3, 4))
  tib <- df %>% tibble::as_tibble()
  mat <- df %>% as.matrix()

  kern_to_cov(df, "SE", hp) %>% expect_equal(kern_to_cov(mat, "SE", hp))
  kern_to_cov(df, "SE", hp) %>% expect_equal(kern_to_cov(tib, "SE", hp))
  kern_to_cov(tib, "SE", hp) %>% expect_equal(kern_to_cov(mat, "SE", hp))
})

test_that("1D-matrix and vector work the same", {
  hp <- tibble::tibble(se_variance = 2, se_lengthscale = 1)
  vec <- c(1, 2, 3)
  mat <- vec %>% as.matrix()

  kern_to_cov(vec, "SE", hp) %>% expect_equal(kern_to_cov(mat, "SE", hp))
})

test_that("dimension names are correct", {
  hp <- tibble::tibble(se_variance = 2, se_lengthscale = 1)
  df <- data.frame(Input = c(5, 6, 7), Cov1 = c(2, 3, 4),
                   Reference = c('5:2', '6:3', '7:4'))
  df2 <- data.frame(Cov1 = c(2, 3, 4), Reference = c('5:2', '6:3', '7:4'),
                    Input = c(5, 6, 7))
  df3 <- data.frame(c(5, 6, 7), c(2, 3, 4))
  df4 <- data.frame(fu = c(5, 6, 7), blob = c(2, 3, 4))

  dimnames(kern_to_cov(df, "SE", hp)) %>%
    expect_equal(dimnames(kern_to_cov(df2, "SE", hp)))
  dimnames(kern_to_cov(df, "SE", hp)) %>%
    expect_equal(dimnames(kern_to_cov(df3, "SE", hp)))
  dimnames(kern_to_cov(df, "SE", hp)) %>%
    expect_equal(dimnames(kern_to_cov(df4, "SE", hp)))
})

test_that("kern_to_cov() works for custom kernels", {
  hp_se <- tibble::tibble(se_variance = 2, se_lengthscale = 1)
  hp_perio <- tibble::tibble(
    perio_variance = 2, perio_lengthscale = 1,
    period = 1
  )
  hp_rq <- tibble::tibble(rq_variance = 2, rq_lengthscale = 1, rq_scale = 1)
  df <- data.frame(Input = c(5, 6, 7), Cov1 = c(2, 3, 4))

  kern_to_cov(df, "SE", hp_se) %>%
    expect_equal(kern_to_cov(df, se_kernel, hp_se))
  kern_to_cov(df, "RQ", hp_rq) %>%
    expect_equal(kern_to_cov(df, rq_kernel, hp_rq))
  kern_to_cov(df, "PERIO", hp_perio) %>%
    expect_equal(kern_to_cov(df, perio_kernel, hp_perio))
})

test_that("kern_to_cov() works for derivative matrices", {
  hp_se <- tibble::tibble(se_variance = 2, se_lengthscale = 1)
  hp_perio <- tibble::tibble(
    perio_variance = 2, perio_lengthscale = 1,
    period = 1
  )
  hp_rq <- tibble::tibble(rq_variance = 2, rq_lengthscale = 1, rq_scale = 1)
  df <- data.frame(Input = c(5, 6, 7), Cov1 = c(2, 3, 4))

  kern_to_cov(df, "SE", hp_se) %>%
    expect_equal(kern_to_cov(df, "SE", hp_se, "se_variance"))
  kern_to_cov(df, "RQ", hp_rq) %>%
    expect_equal(kern_to_cov(df, "RQ", hp_rq, "rq_variance"))
  kern_to_cov(df, "PERIO", hp_perio) %>%
    expect_equal(kern_to_cov(df, "PERIO", hp_perio, "perio_variance"))
  ## Test for custom kernel
  kern_to_cov(df, "PERIO", hp_perio) %>%
    expect_equal(kern_to_cov(df, perio_kernel, hp_perio, "perio_variance"))
})

test_that("kern_to_cov() works for compound kernels", {
  hp_se <- tibble::tibble(se_variance = 2, se_lengthscale = 1)
  hp_perio <- tibble::tibble(
    perio_variance = 2, perio_lengthscale = 1,
    period = 1
  )
  hp_rq <- tibble::tibble(rq_variance = 2, rq_lengthscale = 1, rq_scale = 1)
  hp <- hp_se %>%
    dplyr::bind_cols(hp_perio) %>%
    dplyr::bind_cols(hp_rq)
  df <- data.frame(Input = c(5, 6, 7), Cov1 = c(2, 3, 4))

  kern_to_cov(df, "SE + RQ + PERIO", hp) %>%
    expect_equal(kern_to_cov(df, "SE", hp_se) +
      kern_to_cov(df, "RQ", hp_rq) +
      kern_to_cov(df, "PERIO", hp_perio))

  kern_to_cov(df, "SE * RQ * PERIO", hp) %>%
    expect_equal(kern_to_cov(df, "SE", hp_se) *
      kern_to_cov(df, "RQ", hp_rq) *
      kern_to_cov(df, "PERIO", hp_perio))

  kern_to_cov(df, "SE * RQ + PERIO", hp) %>%
    expect_equal(kern_to_cov(df, "SE", hp_se) *
      kern_to_cov(df, "RQ", hp_rq) +
      kern_to_cov(df, "PERIO", hp_perio))
})

test_that("kern_to_cov() works for compound kernels' derivatives", {
  hp_se <- tibble::tibble(se_variance = 2, se_lengthscale = 1)
  hp_perio <- tibble::tibble(
    perio_variance = 2, perio_lengthscale = 1,
    period = 1
  )
  hp_rq <- tibble::tibble(rq_variance = 2, rq_lengthscale = 1, rq_scale = 1)
  hp_lin <- tibble::tibble(lin_slope = 2, lin_offset = 1)
  hp <- hp_se %>% dplyr::bind_cols(hp_perio, hp_rq, hp_lin)

  df <- data.frame(Input = c(5, 6, 7), Cov1 = c(2, 3, 4))

  kern_to_cov(df, "SE + RQ + PERIO", hp, "se_variance") %>%
    expect_equal(kern_to_cov(df, "SE", hp_se, "se_variance"))

  kern_to_cov(df, "SE * RQ * PERIO * LIN", hp, "rq_lengthscale") %>%
    expect_equal(kern_to_cov(df, "SE", hp_se, ) *
      kern_to_cov(df, "RQ", hp_rq, "rq_lengthscale") *
      kern_to_cov(df, "PERIO", hp_perio) *
      kern_to_cov(df, "LIN", hp_lin))

  kern_to_cov(df, "SE * RQ + PERIO", hp, "period") %>%
    expect_equal(kern_to_cov(df, "PERIO", hp_perio, "period"))

  kern_to_cov(df, "LIN", hp, "lin_offset") %>%
    expect_equal(kern_to_cov(df, "LIN", hp, NULL) -
      kern_to_cov(df, "LIN", hp, "lin_slope"))
})
