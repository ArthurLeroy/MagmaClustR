# Characterization tests for simu_data(), written to protect the internal
# refactoring of the n_outputs == 1 / n_outputs > 1 code paths against
# unintended behaviour changes.
#
# `fixtures/simu_data_reference_pre_refactor.rds` stores the exact output of
# simu_data() (with fixed seeds) captured *before* the refactoring. The
# multi-output refactor intentionally changed the identifier columns
# ('Task_ID', 'Output_ID', 'Input_ID') from factor to character; that type
# change is accommodated below so the tests still guard the (unchanged)
# numeric behaviour and the overall structure.

reference <- readRDS(
  testthat::test_path("fixtures", "simu_data_reference_pre_refactor.rds")
)

.id_cols <- c("Task_ID", "Output_ID", "Input_ID", "ID")
.chr_ids <- function(d) {
  dplyr::mutate(d, dplyr::across(dplyr::any_of(.id_cols), as.character))
}

test_that("simu_data() output structure matches the pre-refactor reference", {
  for (nm in names(reference$configs)) {
    cfg <- reference$configs[[nm]]
    set.seed(reference$seeds[[nm]])
    new_db <- do.call(simu_data, cfg)
    old_db <- reference$results[[nm]]

    expect_s3_class(new_db, "data.frame")
    expect_named(new_db, names(old_db), info = nm)
    expect_equal(nrow(new_db), nrow(old_db), info = nm)

    ## Identifier columns are now character (intended); every other column
    ## keeps its original class.
    expected_cls <- vapply(old_db, function(x) class(x)[1], character(1))
    expected_cls[intersect(.id_cols, names(expected_cls))] <- "character"
    expect_equal(vapply(new_db, function(x) class(x)[1], character(1)),
                 expected_cls, info = nm)

    expect_false(anyNA(new_db$Output), info = nm)
  }
})

test_that("simu_data() output is numerically unchanged for a fixed seed", {
  for (nm in names(reference$configs)) {
    cfg <- reference$configs[[nm]]
    set.seed(reference$seeds[[nm]])
    new_db <- do.call(simu_data, cfg)
    old_db <- reference$results[[nm]]

    ## Compare with identifier columns coerced to character on both sides
    ## (the only intended change); numeric content must be identical.
    expect_equal(.chr_ids(new_db), .chr_ids(old_db), tolerance = 1e-8, info = nm)
  }
})
