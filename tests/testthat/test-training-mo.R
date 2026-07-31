## End-to-end multi-output smoke tests (string output labels).
## Kept intentionally tiny (small data, n_iter_max = 2) to run fast.

.outs <- c("A", "B")
.mo_data <- local({
  set.seed(42)
  raw <- simu_data(n_tasks = 3, n_points = 4, n_outputs = 2,
                   grid_input = seq(-1, 1, 0.25))
  raw %>% dplyr::transmute(
    ID = .data$Task_ID,
    Output_ID = .outs[as.integer(.data$Output_ID)],
    Input = .data$Input, Output = .data$Output
  )
})
.ids <- unique(.mo_data$ID)
.hp0 <- hp(convolution_kernel, list_task_ID = "0", list_output_ID = .outs,
           noise = TRUE) %>% dplyr::select(-.data$Task_ID)
.hpi <- hp(convolution_kernel, list_task_ID = .ids, list_output_ID = .outs,
           shared_hp_tasks = TRUE, noise = TRUE) %>%
  dplyr::rename(ID = .data$Task_ID) %>% dplyr::mutate(ID = as.character(.data$ID))
.hpk <- hp(convolution_kernel, list_task_ID = c("K1", "K2"), list_output_ID = .outs,
           shared_hp_tasks = TRUE, noise = TRUE) %>%
  dplyr::rename(ID = .data$Task_ID) %>% dplyr::mutate(ID = as.character(.data$ID))
.grid <- tidyr::crossing(Output_ID = .outs, Input = seq(-1, 1, 0.5))
.d1 <- .mo_data %>% dplyr::filter(.data$ID == .ids[1]) %>% dplyr::select(-.data$ID)

test_that("train_magma + pred_magma run end-to-end with string output labels", {
  invisible(utils::capture.output(suppressWarnings(
    model <- train_magma(.mo_data, kern_0 = convolution_kernel,
      kern_i = convolution_kernel, ini_hp_0 = .hp0, ini_hp_i = .hpi,
      shared_hp = TRUE, n_iter_max = 2)
  )))
  expect_setequal(as.character(model$hp_0$Output_ID), .outs)
  expect_type(model$hp_i$Output_ID, "character")
  expect_gt(utils::tail(model$seq_loglikelihood, 1), model$seq_loglikelihood[1])

  invisible(utils::capture.output(suppressWarnings(
    pred <- pred_magma(.d1, trained_model = model, grid_inputs = .grid, plot = FALSE)
  )))
  expect_setequal(sort(unique(pred$Output_ID)), .outs)
  expect_true(all(pred$Var >= 0))
  expect_true(all(is.finite(pred$Mean)))
})

test_that("train_magmaclust + pred_magmaclust run end-to-end with string output labels", {
  invisible(utils::capture.output(suppressWarnings(
    mc <- train_magmaclust(.mo_data, nb_cluster = 2, kern_k = convolution_kernel,
      kern_i = convolution_kernel, ini_hp_k = .hpk, ini_hp_i = .hpi,
      shared_hp_k = TRUE, shared_hp_i = TRUE, n_iter_max = 2)
  )))
  expect_setequal(as.character(unique(mc$hp_k$Output_ID)), .outs)
  expect_true(all(c("ID", "K1", "K2") %in% names(mc$hyperpost$mixture)))

  invisible(utils::capture.output(suppressWarnings(
    pred <- pred_magmaclust(.d1, trained_model = mc, grid_inputs = .grid, plot = FALSE)
  )))
  expect_true(all(c("ID", "K1", "K2") %in% names(pred$mixture)))
  expect_setequal(sort(unique(pred$mixture_pred$Output_ID)), .outs)
  expect_true(all(pred$mixture_pred$Var >= 0))
})
