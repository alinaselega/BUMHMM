context("Generating nucleobase patterns")

test_that("function generates nucleobase patterns of length n", {
    expect_equal(nchar(nuclPerm(1)[1]), 1)
    expect_equal(nchar(nuclPerm(2)[1]), 2)
    expect_equal(all(sapply(nuclPerm(2), function(x) nchar(x) == 2)), TRUE)
})

test_that("function asks for positive length", {
    expect_error(nuclPerm(0),
                 "The length of patterns provided is not a positive number.")
    expect_error(nuclPerm(-1),
                 "The length of patterns provided is not a positive number.")
})
