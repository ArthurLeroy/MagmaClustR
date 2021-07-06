test_that("kern_to_cov() works for scalar inputs", {
  hp <- tibble::tibble(variance = 2, lengthscale = 1)
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
  expect_equal(kern_to_cov(input, "SE", hp), res)
})

test_that("kern_to_cov() works for vector inputs", {
  hp <- tibble::tibble(variance = 2, lengthscale = 1)
  input <- data.frame(Input = c(1, 2, 3), Cov1 = c(2, 3, 4))

  res <- matrix(NA, ncol = 3, nrow = 3)
  for (i in 1:3)
  {
    for (j in 1:3) {
      res[i, j] <- se_kernel(input[i, ], input[j, ], hp)
    }
  }
  res = res %>%
    `rownames<-`(input$Input) %>%
    `colnames<-`(input$Input)
  expect_equal(kern_to_cov(input, "SE", hp), res)
})

test_that("matrix, dataframe and tibble work the same", {
  hp <- tibble::tibble(variance = 2, lengthscale = 1)
  df <- data.frame(Input = c(1, 2, 3), Cov1 = c(2, 3, 4))
  tib <- df %>% tibble::as_tibble()
  mat <- df %>% as.matrix()

  expect_equal(kern_to_cov(df, "SE", hp), kern_to_cov(mat, "SE", hp))
  expect_equal(kern_to_cov(df, "SE", hp), kern_to_cov(tib, "SE", hp))
  expect_equal(kern_to_cov(tib, "SE", hp), kern_to_cov(mat, "SE", hp))
})

test_that("1D-matrix and vector work the same", {
  hp <- tibble::tibble(variance = 2, lengthscale = 1)
  vec <- c(1, 2, 3)
  mat <- vec %>% as.matrix()

  expect_equal(kern_to_cov(vec, "SE", hp), kern_to_cov(mat, "SE", hp))
})

test_that("dimension names are correct", {
  hp <- tibble::tibble(variance = 2, lengthscale = 1)
  df <- data.frame(Input = c(5, 6, 7), Cov1 = c(2, 3, 4))
  df2 <- data.frame(Cov1 = c(2, 3, 4), Input = c(5, 6, 7))
  df3 <- data.frame(c(5, 6, 7), c(2, 3, 4))
  df4 <- data.frame(fu = c(5, 6, 7), blob = c(2, 3, 4))

  expect_equal(dimnames(kern_to_cov(df, "SE", hp)),
               dimnames(kern_to_cov(df2, "SE", hp)))
  expect_equal(dimnames(kern_to_cov(df, "SE", hp)),
               dimnames(kern_to_cov(df3, "SE", hp)))
  expect_equal(dimnames(kern_to_cov(df, "SE", hp)),
               dimnames(kern_to_cov(df4, "SE", hp)))
})
