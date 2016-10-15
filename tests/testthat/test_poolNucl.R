context("Extracting nucleotides selected for comparisons")

covFile = matrix(c(1,2,4,6,0,3,7,0,2,5,1,4), nrow = 3, ncol = 4)
docFile = matrix(c(1,0,1,5,0,1,5,0,1,4,0,0), nrow = 3, ncol = 4)

dorFile = docFile / covFile
dorFile[is.na(dorFile)] <- 0

nuclSelection <- selectNuclPos(covFile, docFile, 2, 2, 1)

comparisons <- t(combn(2, 2))
comparisons2 <- t(matrix(c(rep((2+1):(2+2), each=2), rep(1:2, 2)), 2,
                         byrow=TRUE))

test_that("function correctly extracts positions in control replicates that are
          selected for comparisons", {
    expect_equal(poolNucl(1, comparisons, nuclSelection$analysedC),
                 c(1, 3))

    expect_equal(poolNucl(1, comparisons2, nuclSelection$analysedCT),
                 c(1, 3))
})

test_that("function correctly extracts positions in treatment replicates that
          are selected for comparisons", {
    expect_equal(poolNucl(4, comparisons2, nuclSelection$analysedCT),
                         c(1))
})

