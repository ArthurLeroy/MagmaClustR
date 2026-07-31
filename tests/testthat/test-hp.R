test_that("hp() convolution returns character identifier columns", {
  h <- hp(convolution_kernel, list_task_ID = c("a", "b"),
          list_output_ID = c("A", "B"), noise = TRUE)
  expect_type(h$Task_ID, "character")
  expect_type(h$Output_ID, "character")
  expect_setequal(unique(h$Output_ID), c("A", "B"))
})

test_that("hp() convolution: sharing flags control l_t/S_t granularity", {
  args <- list(kern = convolution_kernel, list_task_ID = c("a", "b"),
               list_output_ID = c("A", "B"), noise = TRUE)
  h_TT <- do.call(hp, c(args, shared_hp_tasks = TRUE,  shared_hp_outputs = TRUE))
  h_TF <- do.call(hp, c(args, shared_hp_tasks = TRUE,  shared_hp_outputs = FALSE))
  h_FT <- do.call(hp, c(args, shared_hp_tasks = FALSE, shared_hp_outputs = TRUE))
  h_FF <- do.call(hp, c(args, shared_hp_tasks = FALSE, shared_hp_outputs = FALSE))
  expect_equal(dplyr::n_distinct(h_TT$l_t), 1)  # shared everywhere
  expect_equal(dplyr::n_distinct(h_TF$l_t), 2)  # per output
  expect_equal(dplyr::n_distinct(h_FT$l_t), 2)  # per task
  expect_equal(dplyr::n_distinct(h_FF$l_t), 4)  # per (task, output)
})

test_that("hp() convolution: l_u_t is constant across outputs within a task", {
  h <- hp(convolution_kernel, list_task_ID = c("a", "b"),
          list_output_ID = c("A", "B", "C"), shared_hp_tasks = FALSE, noise = TRUE)
  nd <- h %>% dplyr::group_by(Task_ID) %>%
    dplyr::summarise(nd = dplyr::n_distinct(l_u_t), .groups = "drop")
  expect_true(all(nd$nd == 1))
  expect_equal(dplyr::n_distinct(h$l_u_t), 2)  # one per task
})

test_that("hp() convolution: n_latent structure and constances", {
  h <- hp(convolution_kernel, list_task_ID = c("a", "b"),
          list_output_ID = c("A", "B"), n_latent = 2,
          shared_hp_tasks = TRUE, noise = TRUE)
  expect_true("Latent_ID" %in% names(h))
  expect_type(h$Latent_ID, "character")
  expect_equal(dplyr::n_distinct(h$Latent_ID), 2)
  expect_equal(nrow(h), 2 * 2 * 2)  # tasks * outputs * latents
  ndl <- h %>% dplyr::group_by(Task_ID, Latent_ID) %>%
    dplyr::summarise(nd = dplyr::n_distinct(l_u_t), .groups = "drop")
  expect_true(all(ndl$nd == 1))  # l_u_t constant across outputs within a latent
  ndn <- h %>% dplyr::group_by(Task_ID, Output_ID) %>%
    dplyr::summarise(nd = dplyr::n_distinct(noise), .groups = "drop")
  expect_true(all(ndn$nd == 1))  # noise per output, constant across latents
})

test_that("hp() string kernel keeps single-output behaviour", {
  h <- hp("SE", noise = TRUE)
  expect_setequal(names(h), c("se_variance", "se_lengthscale", "noise"))
  expect_equal(nrow(h), 1)
})
