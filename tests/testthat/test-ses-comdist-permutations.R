test_that("ses.comdistnt2 restores dimnames after stratified permatswap", {
  samp <- matrix(
    c(
      3, 0, 1,
      0, 2, 1,
      1, 1, 0
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(
      c("sample_a", "sample_b", "sample_c"),
      c("taxon_1", "taxon_2", "taxon_3")
    )
  )
  dis <- matrix(
    c(
      0, 1, 2,
      1, 0, 1.5,
      2, 1.5, 0
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(colnames(samp), colnames(samp))
  )

  local_mocked_bindings(
    permatswap = function(...) list(perm = list(unname(samp))),
    .package = "MicEco"
  )

  expect_no_error(
    result <- ses.comdistnt2(
      samp,
      dis,
      method = "quasiswap",
      strata = c(1, 1, 2),
      runs = 1,
      cores = 1
    )
  )
  expect_equal(rownames(result$comdistnt.obs), rownames(samp))
  expect_equal(colnames(result$comdistnt.obs), rownames(samp))
})

test_that("ses.comdist2 restores dimnames after stratified permatswap", {
  samp <- matrix(
    c(
      3, 0, 1,
      0, 2, 1,
      1, 1, 0
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(
      c("sample_a", "sample_b", "sample_c"),
      c("taxon_1", "taxon_2", "taxon_3")
    )
  )
  dis <- matrix(
    c(
      0, 1, 2,
      1, 0, 1.5,
      2, 1.5, 0
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(colnames(samp), colnames(samp))
  )

  local_mocked_bindings(
    permatswap = function(...) list(perm = list(unname(samp))),
    .package = "MicEco"
  )

  expect_no_error(
    result <- ses.comdist2(
      samp,
      dis,
      method = "quasiswap",
      strata = c(1, 1, 2),
      runs = 1,
      cores = 1
    )
  )
  expect_equal(rownames(result$comdist.obs), rownames(samp))
  expect_equal(colnames(result$comdist.obs), rownames(samp))
})
