computeProbs <- function(LDR_C, LDR_CT, Nc, Nt, strand, nuclPosition, analysedC,
                         analysedCT, stretches, optimise=NULL) {

    if ((Nc < 2) | (Nt < 2)) {
        stop('Number of control and treatment replicates must be at least 2.')
    }
    else if (any(sapply(analysedC, function(x) length(x) == 0))) {
        stop('All lists of positions selected for pair-wise comparisons should
              be non-empty.')
    }
    else if (any(sapply(analysedCT, function(x) length(x) == 0))) {
        stop('All lists of positions selected for pair-wise comparisons should
              be non-empty.')
    }
    if ((strand != '+') & (strand != '-')) {
        stop('Strand should be either plus or minus, specified with a sign.')
    }
    if (length(nuclPosition) == 0) {
        stop('The list of considered nucleotide positions should be non-empty.
              If patterns are not used, a vector with all positions should be
              provided in the first element of the list.')
    }
    if (dim(LDR_C)[2] != length(analysedC)) {
        stop('The matrix of control-control LDRs should have as many columns
              as there are control-control comparisons.')
    }
    if (dim(LDR_CT)[2] != length(analysedCT)) {
      stop('The matrix of treatment-control LDRs should have as many columns
                as there are treatment-control comparisons.')
    }
    if (!is.null(optimise) & !is.numeric(optimise)) {
        stop('Please provide a tolerance if shape parameters are to be optimised
             with EM algorithm.')
    }
    else {

        ## Number of nucleotides
        nNucl <- dim(LDR_C)[1]

        quantiles <- list()

        ## Intervals at which to compute quantiles
        intervals <- c(seq(0.01, 0.89, 0.01), seq(0.9, 0.999, 0.001))

        ## Number of intervals
        intervalNum <- length(intervals)

        message('Computing quantiles of null distributions...')

        ## Precompute quantiles of LDR distributions for different patterns
        for (i in 1:length(nuclPosition)) {

            subDistr <- array()
            for (j in 1:dim(LDR_C)[2]) {
                ## Extract a subset of the null distribution, corresponding to
                ## the current pattern
                observedPos <- intersect(unlist(nuclPosition[i]),
                                         analysedC[[j]])
                subDistr <- c(LDR_C[observedPos, j], subDistr)
                if (j == 1) {
                    subDistr <- subDistr[1:(length(subDistr)-1)]
                }
            }

            ## Compute the quantiles
            quantiles[[i]] = as.numeric(stats::quantile(subDistr, intervals))
        }

        ## Annotate with the names of patterns
        names(quantiles) <- names(nuclPosition)

        ## Matrix for p-values
        empPvals <- matrix(, nrow=nNucl, ncol=Nc * Nt)

        ## Compute empirical p-values

        message('Computing empirical p-values...')

        ## For each pattern
        for (i in 1:length(nuclPosition))  {

            ## Collect positions for all treatment-control comparisons
            all_analysedCT <- Reduce(union, analysedCT)

            ## Of those, find which nucleotides correspond to the current
            ## pattern
            patternPos <- unlist(nuclPosition[i])
            patternPos <- intersect(patternPos, all_analysedCT)

            ## Extract treatment-control LDRs at these positions from all
            ## comparisons
            subDist <- LDR_CT[patternPos, ]

            if (is.null(dim(subDist))) {
                dim(subDist) <- c(length(subDist), 1)
            }

            ## Compare each LDR to the quantiles of the null distribution
            ## for that pattern
            pValPattern <- matrix(, nrow=length(patternPos), ncol=Nc * Nt)

            ## Select the closest quantile to the LDR
            pValPattern <- apply(subDist, c(1, 2),
                                 function(x)
                                 if (is.na(x)) {NA}
                                 else {
                                   which.min(abs(unlist(quantiles[[i]]) - x))
                                 })

            ## Compute the p-value as 1 - quantile
            pValPattern <- apply(pValPattern, c(1, 2),
                                 function(x)
                                 1 - intervals[x])

            ## Store the empirical p-values for these positions
            empPvals[patternPos, ] <- pValPattern
        }

        empPvals <- t(empPvals)

        ## Settings for HMM
        ## Set the transition matrix
        trans <- matrix(c(0.95, 0.2, 0.05, 0.8), nrow=2, ncol=2)
        trans <- t(trans)

        ## Set the values for Beta shape parameters in the emission mixture
        ## model
        alpha <- 1
        beta <- 10

        ## Set initial probability to 0.5 in each state
        in_prob <- c(0.5, 0.5)

        ## Matrix for posterior probabilities of modification
        posteriors <- matrix(, nrow=nNucl, ncol=2)

        message('Computing posteriors...')

        for (i in 1:length(stretches)) {

            ## Extract start and end of a current stretch
            stretchStart <- start(stretches)[i]
            stretchEnd <- end(stretches)[i]

            ## Run HMM on the stretch

            ## For plus strand
            if (strand == '+') {
                posterior <- hmmFwbw(empPvals[, stretchStart:stretchEnd],
                                     trans, in_prob, alpha, beta)
                posteriors[stretchStart:stretchEnd, ] <- t(posterior)
            }

            ## For minus strand
            if (strand == '-') {
                posterior = hmmFwbw(empPvals[, stretchEnd:stretchStart],
                                    trans, in_prob, alpha, beta)
                posteriors[stretchStart:stretchEnd, ] <- t(posterior)
            }
        }

        if (!is.null(optimise)) {
            tolerance <- optimise
            optimised_beta_parameters <- betaParamsEM(posteriors, empPvals,
                                                      alpha, beta, tolerance,
                                                      strand, stretches, trans,
                                                      in_prob)
            posteriors <- optimised_beta_parameters$posteriors
        }

        return(posteriors)
    }
}
