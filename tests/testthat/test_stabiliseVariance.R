# BUM-HMM
# Copyright (C) 2016 Alina Selega

context("Coverage bias, stabilising variance")

covFile <- matrix(c(1,2,4,7,0,3,7,0,2,5,1,4), nrow = 3, ncol = 4)
docFile <- matrix(c(1,0,1,5,0,1,5,0,1,4,0,0), nrow = 3, ncol = 4)
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

assay(se, "dropoff_rate") <- scaleDOR(se, nuclSelection, 2, 2)

test_that("function throws error when all quantiles are 0", {
    expect_error(stabiliseVariance(se, nuclSelection, 2, 2),
               "Unable to fit the model for correcting the coverage bias.")
})


covFile <- matrix(c(100,38,22,15,14,117,20,37,12,11,7,0,2,3,4,5,1,4,8,7),
                 nrow = 5, ncol = 4)
docFile <- matrix(c(1,1,1,1,1,5,1,1,3,3,5,0,1,1,1,4,0,0,1,1), nrow = 5, ncol = 4)
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

assay(se, "dropoff_rate") <- scaleDOR(se, nuclSelection, 2, 2)

LDR_CT <- stabiliseVariance(se, nuclSelection, 2, 2)$LDR_CT

test_that("function rescales LDRs correctly", {
    expect_equal(stabiliseVariance(se, nuclSelection, 2, 2)$LDR_C[1,1],
               -7.697804, tolerance=1e-7)
    expect_equal(LDR_CT[1,2], 3.94272, tolerance=1e-7)
})

notSelected <- sapply(1:length(nuclSelection$analysedCT),
          function(i) setdiff(1:dim(covFile)[1], nuclSelection$analysedCT[[i]]))

## Nc = 2 and Nt = 2
indexT <- t(matrix(c(rep((2+1):(2+2), each=2), rep(1:2, 2)),
                   2, byrow=TRUE))

test_that("function only rescales LDRs from selected positions", {
    expect_equal(all(is.na(LDR_CT[notSelected[[1]], 1])), TRUE)
    expect_equal(all(is.na(LDR_CT[notSelected[[2]], 2])), TRUE)
    expect_equal(all(is.na(LDR_CT[notSelected[[3]], 3])), TRUE)
    expect_equal(all(is.na(LDR_CT[notSelected[[4]], 4])), TRUE)
})

test_that("function requires non-empty lists of selected positions", {
    expect_error(stabiliseVariance(se, list(nuclSelection$analysedC, list()),
                                   2, 2),
                 "All lists of positions selected for pair-wise comparisons should
             be non-empty.")

    expect_error(stabiliseVariance(se, list(nuclSelection$analysedC,
                                            list(
                                                c(1,3), list()
                                                )),
                                   2, 2),
                 "All lists of positions selected for pair-wise comparisons should
             be non-empty.")
})
