context("Finding patterns in sequence")

sequence <- 'TGACGTG'

patterns <- nuclPerm(2)

test_that("function correctly finds patterns", {
    expect_equal(all(findPatternPos(patterns, sequence, '+')[[2]] == c(3, 4)),
                 TRUE)
    expect_equal(all(findPatternPos(patterns, sequence, '+')[[15]][1,]
                     == c(1, 2)), TRUE)
    expect_equal(all(findPatternPos(patterns, sequence, '+')[[15]][2,]
                     == c(6, 7)), TRUE)
    expect_equal(all(is.na(findPatternPos(patterns, sequence, '+')[[1]])), TRUE)
})

test_that("function requires strand to be '+' or '-'", {
    expect_error(findPatternPos(patterns, sequence, 'plus'),
              "Strand should be either plus or minus, specified with a sign.")
})

test_that("function requires a non-empty sequence", {
    expect_error(findPatternPos(patterns, '', '+'),
                 "The sequence should be non-empty.")
})

test_that("function requires a non-empty list of patterns", {
    expect_error(findPatternPos(list(), sequence, '+'),
                 "The list of patterns should be non-empty.")
})
