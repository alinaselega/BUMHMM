context("Compute probabilities")

library(SummarizedExperiment)

covFile <-
  matrix(c(100,50,60,50,70,80,60,50,60,50,80,80,70,70,60,70,100,80,80,70),
                 nrow = 5, ncol = 4)
docFile <- matrix(c(10,2,5,10,12,8,3,5,8,7,60,60,50,50,50,50,60,50,50,40),
                 nrow = 5, ncol = 4)

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

LDR_C <- stabiliseVariance(se, nuclSelection, 2, 2)$LDR_C
LDR_CT <- stabiliseVariance(se, nuclSelection, 2, 2)$LDR_CT

nuclPosition <- list()
nuclPosition[[1]] <- 1:dim(covFile)[1]

stretches <- computeStretches(se, 1)

test_that("function requires a non-empty list of considered nucleotides", {
    expect_error(computeProbs(LDR_C, LDR_CT, 2, 2, '+', list(),
                            nuclSelection$analysedC,
                            nuclSelection$analysedCT, stretches),
               "The list of considered nucleotide positions should be non-empty.
              If patterns are not used, a vector with all positions should be
              provided in the first element of the list.")
})

test_that("function expects LDR_CT to have as many columns as CT comparisons", {
    expect_error(computeProbs(LDR_C, LDR_CT[,1:2], 2, 2, '+', nuclPosition,
                            nuclSelection$analysedC,
                            nuclSelection$analysedCT, stretches),
                 "The matrix of treatment-control LDRs should have as many columns
                as there are treatment-control comparisons.")
})

test_that("function computes probabilities correctly", {
    expect_equal(computeProbs(LDR_C, LDR_CT, 2, 2, '+', nuclPosition,
                            nuclSelection$analysedC,
                            nuclSelection$analysedCT, stretches)[2,2], 1,
               tolerance=1e-5)
    expect_equal(computeProbs(LDR_C, LDR_CT, 2, 2, '+', nuclPosition,
                            nuclSelection$analysedC,
                            nuclSelection$analysedCT, stretches)[3,2], 0,
               tolerance=1e-7)
})

test_that("function requires a tolerance if optimisation is to be used", {
    expect_error(computeProbs(LDR_C, LDR_CT, 2, 2, '+', nuclPosition,
                              nuclSelection$analysedC,
                              nuclSelection$analysedCT, stretches, TRUE),
                 "Please provide a tolerance if shape parameters are to be optimised
             with EM algorithm.")
})
