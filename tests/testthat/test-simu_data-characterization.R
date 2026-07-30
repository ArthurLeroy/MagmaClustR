# Characterization tests for simu_data(), written to protect the internal
# refactoring of the n_outputs == 1 / n_outputs > 1 code paths against
# unintended behaviour changes.
#
# `fixtures/simu_data_reference_pre_refactor.rds` stores the exact output of
# simu_data() (with fixed seeds) captured *before* the refactoring. If the
# refactor preserves the exact sequence of random draws, the outputs should
# remain bit-identical; if not, at least the structural invariants below must
# still hold, and any numerical discrepancy should be investigated.

reference <- readRDS(
  testthat::test_path("fixtures", "simu_data_reference_pre_refactor.rds")
)

test_that("simu_data() output structure matches the pre-refactor reference", {
  for (nm in names(reference$configs)) {
    cfg <- reference$configs[[nm]]
    set.seed(reference$seeds[[nm]])
    new_db <- do.call(simu_data, cfg)
    old_db <- reference$results[[nm]]

    expect_s3_class(new_db, "data.frame")
    expect_named(new_db, names(old_db), info = nm)
    expect_equal(nrow(new_db), nrow(old_db), info = nm)
    expect_equal(
      vapply(new_db, class, character(1)),
      vapply(old_db, class, character(1)),
      info = nm
    )
    expect_false(anyNA(new_db$Output), info = nm)
  }
})

test_that("simu_data() output is numerically unchanged for a fixed seed", {
  for (nm in names(reference$configs)) {
    cfg <- reference$configs[[nm]]
    set.seed(reference$seeds[[nm]])
    new_db <- do.call(simu_data, cfg)
    old_db <- reference$results[[nm]]

    expect_equal(new_db, old_db, tolerance = 1e-8, info = nm)
  }
})
