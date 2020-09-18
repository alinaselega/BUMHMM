# BUM-HMM
# Copyright (C) 2016 Alina Selega

context("Computing stretches of nucleotides")

covFile <- matrix(c(1,2,4,6,1,3,7,0,2,5,1,4,1,2,4,6,1,2,4,6,1,2,4,6),
                  nrow=6, ncol=4)
docFile <- matrix(c(1,0,1,5,1,1,5,0,1,4,1,1,1,1,1,5,1,1,1,5,1,1,1,5),
                  nrow=6, ncol=4)
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

test_that("function computes stretches correctly", {
    expect_equal(computeStretches(se, 1), IRanges(start = 3, end = 6))
})

test_that("function does not allow negative coverage threshold", {
    expect_error(computeStretches(se, -1),
                 "The minumum coverage threshold must be non-negative.")
})
