# BUM-HMM
# Copyright (C) 2016 Alina Selega, Christel Sirocchi, Sander Granneman, Guido Sanguinetti

## This function is for internal use by computeProbs.R

hmmFwbw <- function(pVals, trans, initialProb, alpha, beta) {

    if (is.null(dim(pVals))) {
        dim(pVals) <- c(1, length(pVals))
    }

    ## Number of experiments (treatment-control comparisons)
    nexp <- length(pVals[, 1])

    ## Number of nucleotides
    nBins <- length(pVals[1, ])

    ## Number of hidden states
    nStates <- length(trans[, 1])

    ## Log-likelihood of observations given state
    obsLike <- matrix(1, ncol=nBins, nrow=nStates)

    fwdMessage <- matrix(0, ncol=nBins, nrow=nStates)
    bwdMessage <- matrix(0, ncol=nBins, nrow=nStates)

    ## Calculation of likelihoods (sum over replicates of likelihoods of each
    ## experiment)
    for (index in 1:nexp) {
        for (index2 in 1:nBins) {
            if (is.na(pVals[index, index2])) {
                obsLike[, index2] <- obsLike[, index2]
            } else {
                obsLike[, index2] <- obsLike[, index2] *
                stats::dbeta(pVals[index, index2], c(1, alpha), c(1, beta))
            }
        }
    }

    ## Calculation of the forward messages
    fwdMessage[, 1] <- initialProb * obsLike[, 1]
    fwdMessage[, 1] <- fwdMessage[, 1] / sum(fwdMessage[, 1])

    for (index in 2:nBins) {
        fwdMessage[, index] <- (trans %*% fwdMessage[, index - 1]) *
                               obsLike[, index]
        fwdMessage[, index] <- fwdMessage[, index] / sum(fwdMessage[, index])
    }

    ## Calculation of the backward message
    bwdMessage[, nBins] <- 1

    for (index in (nBins - 1):1) {
        bwdMessage[, index] <- (trans %*% bwdMessage[, index + 1]) *
                               obsLike[, index]
        bwdMessage[, index] <- bwdMessage[, index] / sum(bwdMessage[, index])
    }

    ## Calculation of posteriors
    posterior <- fwdMessage * bwdMessage
    posterior <- posterior / (matrix(1, nrow=length(posterior[, 1]))
                              %*% colSums(posterior))
    return(posterior)
}
