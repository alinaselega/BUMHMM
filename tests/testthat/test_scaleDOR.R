# BUM-HMM
# Copyright (C) 2016 Alina Selega

context("Scaling drop-off rates")

covFile <- matrix(c(10,15,20,8,0,18,12,0,14,14,10,20), nrow = 3, ncol = 4)
docFile <- matrix(c(2,0,4,4,0,9,3,0,7,14,8,12), nrow = 3, ncol = 4)
dorFile <- docFile / covFile
dorFile[is.na(dorFile)] <- 0

se <- SummarizedExperiment(
  list(
    coverage=as.matrix(covFile),
    dropoff_count=as.matrix(docFile),
    dropoff_rate=as.matrix(dorFile)
  ), colData=DataFrame(
    replicate=rep(c("control", "treatment"), each=2)
  )
)

nuclSelection <- selectNuclPos(se, 2, 2, 1)

test_that("function correctly scales drop-off rates", {
    expect_equal(scaleDOR(se, nuclSelection, 2, 2)[1,1],
                 0.35, tolerance=1e-7)

    expect_equal(scaleDOR(se, nuclSelection, 2, 2)[3,2],
                 0.35, tolerance=1e-7)
})

test_that("function expects more than 2 replicates", {
    expect_error(scaleDOR(se, nuclSelection, 1, 2),
              "Number of control and treatment replicates must be at least 2.")
    expect_error(scaleDOR(se, nuclSelection, 2, 1),
              "Number of control and treatment replicates must be at least 2.")
})

test_that("function expects nuclSelection with 2 elements", {
    expect_error(scaleDOR(se, list(), 2, 2),
                 "Nucleotide selection should have two elements.")
})

test_that("function expects nuclSelection to have non-empty elements", {
    expect_error(scaleDOR(se, list(nuclSelection[[1]],
                                      list()), 2, 2),
               "All lists of positions selected for pair-wise comparisons should be non-empty.")
})

dorFile[1,1] <- NA
assay(se, "dropoff_rate") <- as.matrix(dorFile)

test_that("function does not allow any NAs in the data matrix", {
    expect_error(scaleDOR(se, nuclSelection, 2, 2),
                 "Drop-off rate matrix should not have any NA entries.")
})

dorFile <- matrix(0, nrow = 3, ncol = 4)
assay(se, "dropoff_rate") <- as.matrix(dorFile)

test_that("function works without error if all drop-off rates are 0", {
    expect_equal(all(is.na(scaleDOR(se, nuclSelection, 2, 2))), TRUE)
})
