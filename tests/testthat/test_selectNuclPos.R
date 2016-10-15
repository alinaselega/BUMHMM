context("Selecting nucleotide positions")

covFile = matrix(c(1,2,4,6,0,3,7,0,2,5,1,4), nrow = 3, ncol = 4)
docFile = matrix(c(1,0,1,5,0,1,5,0,1,4,0,0), nrow = 3, ncol = 4)

test_that("function selects correct nucleotides for control-control pairs", {
    expect_equal(selectNuclPos(covFile, docFile, 2, 2, 1)$analysedC,
                 list(c(1, 3)))
})

test_that("function selects correct nucleotides for treatment-control pairs", {
    expect_equal(selectNuclPos(covFile, docFile, 2, 2, 1)$analysedCT,
                 list(c(1, 3),
                      c(1, 3),
                      c(1),
                      c(1)))
})

test_that("function selects correct nucleotides for posteriors", {
    expect_equal(selectNuclPos(covFile, docFile, 2, 2, 1)$computePosteriors,
                 list(c(1, 3),
                      c(1)))
})

test_that("function does not allow less than two replicates", {
    expect_error(selectNuclPos(covFile, docFile, 1, 2, 1),
                 "The number of experimental replicates must be at least 2.")
})

test_that("function does not allow negative coverage threshold", {
    expect_error(selectNuclPos(covFile, docFile, 2, 2, -1),
                 "The minumum coverage threshold must be non-negative.")
})

covFile[1,1] <- NA

test_that("function does not allow NA in data matrices", {
    expect_error(selectNuclPos(covFile, docFile, 2, 2, 1),
               "The coverage and drop-off count matrices should not have NA
              entries.")
})

covFile = matrix(0, nrow = 3, ncol = 4)
docFile = matrix(0, nrow = 3, ncol = 4)

test_that("function works without error if all entries in data matrix are 0", {
    expect_equal(selectNuclPos(covFile, docFile, 2, 2, 1)$analysedC,
               list(integer(0)))

    expect_equal(selectNuclPos(covFile, docFile, 2, 2, 1)$analysedCT,
               list(integer(0),
                    integer(0),
                    integer(0),
                    integer(0)))

    expect_equal(selectNuclPos(covFile, docFile, 2, 2, 1)$computePosteriors,
               list(integer(0),
                    integer(0)))
})

