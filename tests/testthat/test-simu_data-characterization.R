# Golden-output (characterization) tests for simu_data().
#
# `fixtures/simu_data_reference.rds` stores the exact output of simu_data() for
# a set of fixed configs and seeds, guarding against unintended changes to the
# simulation logic. simu_data() draws its samples with a Cholesky factorisation
# of the covariance (a canonical, sign-fixed matrix square root), so the exact
# Output values are reproducible across BLAS/LAPACK backends, OSes and R
# versions. This makes the numeric test below portable -- unlike an eigen-based
# sampler, whose eigenvector signs vary across environments.

reference <- readRDS(
  testthat::test_path("fixtures", "simu_data_reference.rds")
)

# simu_data() returns the ID columns (Task_ID, Input_ID, Output_ID) as
# character. Normalise any factor columns to character before comparing, so the
# comparison is robust to that representation.
coerce_factor_to_chr <- function(d) {
  is_fac <- vapply(d, is.factor, logical(1))
  d[is_fac] <- lapply(d[is_fac], as.character)
  d
}

test_that("simu_data() output structure matches the reference", {
  for (nm in names(reference$configs)) {
    cfg <- reference$configs[[nm]]
    set.seed(reference$seeds[[nm]])
    new_db <- coerce_factor_to_chr(do.call(simu_data, cfg))
    old_db <- coerce_factor_to_chr(reference$results[[nm]])

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

test_that("simu_data() output is numerically reproducible for a fixed seed", {
  for (nm in names(reference$configs)) {
    cfg <- reference$configs[[nm]]
    set.seed(reference$seeds[[nm]])
    new_db <- coerce_factor_to_chr(do.call(simu_data, cfg))
    old_db <- coerce_factor_to_chr(reference$results[[nm]])

    expect_equal(new_db, old_db, tolerance = 1e-6, info = nm)
  }
})
