test_that("kern_to_inv() works for scalar inputs", {
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
    solve() %>%
    `rownames<-`(as.character(input)) %>%
    `colnames<-`(as.character(input))
  expect_equal(kern_to_inv(input, "SE", hp), res)
})

test_that("kern_to_inv() works for vector inputs", {
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
    solve() %>%
    `rownames<-`(ref) %>%
    `colnames<-`(ref)

  kern_to_inv(input, "SE", hp) %>% expect_equal(res)
})

test_that("dimension names are correct", {
  hp <- tibble::tibble(se_variance = 2, se_lengthscale = 1)
  df <- data.frame(Input = c(5, 6, 7), Cov1 = c(2, 3, 4),
                   Reference = c('5:2', '6:3', '7:4'))
  df2 <- data.frame(Cov1 = c(2, 3, 4), Reference = c('5:2', '6:3', '7:4'),
                    Input = c(5, 6, 7))
  df3 <- data.frame(c(5, 6, 7), c(2, 3, 4))
  df4 <- data.frame(fu = c(5, 6, 7), blob = c(2, 3, 4))

  dimnames(kern_to_inv(df, "SE", hp)) %>%
    expect_equal(list(c('5:2', '6:3', '7:4'), c('5:2', '6:3', '7:4')))
  dimnames(kern_to_inv(df, "SE", hp)) %>%
    expect_equal(dimnames(kern_to_inv(df2, "SE", hp)))
  dimnames(kern_to_inv(df, "SE", hp)) %>%
    expect_equal(dimnames(kern_to_inv(df3, "SE", hp)))
  dimnames(kern_to_inv(df, "SE", hp)) %>%
    expect_equal(dimnames(kern_to_inv(df4, "SE", hp)))
})
