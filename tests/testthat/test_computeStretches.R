context("Computing stretches of nucleotides")

test_that("function computes stretches correctly", {
    expect_equal(computeStretches(c(1, 2, 3, 4, 7, 8, 9, 15, 17, 18)),
                                  list(c(1, 4),
                                       c(7, 9),
                                       c(17, 18)))
    expect_equal(computeStretches(c(1, 2)), list(c(1, 2)))
    expect_equal(computeStretches(c(1, 3)), list())
    expect_equal(computeStretches(c(1, 3, 4)), list(c(3, 4)))
})

test_that("function requires at least two indices", {
    expect_error(computeStretches(c(1)),
                 "There should be at least 2 nucleotide positions in the list.")
})
