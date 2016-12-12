\name{computeProbs}
\alias{computeProbs}
\title{
    Function to compute posterior probabilities of modification with the BUM-HMM
    model for all nucleotides.
}
\description{
    This function computes posterior probabilities of chemical modification in
    a structure probing experiment. It uses empirical p-values, which quantify
    the difference between LDRs (log-ratios of drop-off rates) from
    treatment-control and control-control experimental comparisons, as
    observations in a HMM with an emission model defined as a Beta-Uniform
    mixture.
}
\usage{
    computeProbs(LDR_C, LDR_CT, Nc, Nt, strand, nuclPosition, analysedC,
                         analysedCT, stretches, optimise)
}
\arguments{
    \item{LDR_C}{
        A matrix of transformed LDRs for control-control comparisons as returned
        by \code{\link[BUMHMM]{stabiliseVariance}}. The matrix rows correspond
        to nucleotide positions and columns to control-control comparisons.
    }
    \item{LDR_CT}{
        A matrix of transformed LDRs for treatment-control comparisons as
        returned by \code{\link[BUMHMM]{stabiliseVariance}}. The matrix rows
        correspond to nucleotide positions and columns to treatment-control
        comparisons.
    }
    \item{Nc}{
        Number of control experimental replicates. Must be at least 2.
    }
    \item{Nt}{
        Number of treatment experimental replicates. Must be at least 2.
    }
    \item{strand}{
        A character, indicating the plus (\code{+}) or minus strand (\code{-}).
    }
    \item{nuclPosition}{
        A list where each component corresponds to a pattern (indicated by the
        field \code{names}) and contains indices of the middle nucleotides of
        that pattern's occurrences within the sequence, as returned by
        \code{\link[BUMHMM]{findPatternPos}}. If the sequence bias is not
        addressed and all nucleotide positions are used together (e.g. for short
        transcripts), the only component of the list should contain all
        nucleotide positions of the sequence.
    }
    \item{analysedC}{
        The first element of the list returned by
        \code{\link[BUMHMM]{selectNuclPos}}, containing the positions of
        nucleotides selected for all control-control comparisons.
    }
    \item{analysedCT}{
        The second element of the list returned by
        \code{\link[BUMHMM]{selectNuclPos}}, containing the positions of
        nucleotides selected for all treatment-control comparisons.
    }
    \item{stretches}{
        An \code{\link{IRanges}} object storing each uninterrupted stretch of
        nucleotides for which posterior probabilities are to be computed, as
        returned by \code{\link[BUMHMM]{computeStretches}}.
    }
    \item{optimise}{
        An indicator variable to use the EM optimisation of shape parameters of
        the Beta distribution (emission model of the HMM for the modified
        state). If optimisation is to be used, the variable should be set to the
        desired tolerance for parameter estimation (once the estimates stop
        changing within this tolerance, the algorithm stops). Set to \code{NULL}
        by default.
    }
}
\details{
    The function first computes percentiles of LDR null distributions, for each
    percent from 1 to 89 and for a tenth of a percent from 90 to 100. If
    multiple nucleotide patterns are used, then null distributions are computed
    separately for each pattern. Then, LDRs from treatment-control comparisons
    are compared to the percentiles of a corresponding null distribution via
    computing the p-values, defined as \code{1 - closest percentile}. Thus, if
    the closest percentile to an LDR is 99th, then its p-value will be 0.01,
    reflecting a small probability for it to be from this null distribution.
    The p-values are then used as observations in a hidden Markov model.
    P-values corresponding to multiple comparisons at the same nucleotide
    position are considered to be independent measurements.

    HMM has a binary latent state corresponding to the true state of a
    nucleotide, modified or unmodified. Transition probabilities are set
    according to the empirically derived expectations of unmodified nucleotide
    stretches of average length 20 and modified stretches of average length 5.
    Emission probabilites come from a Beta-Uniform mixture, where the Uniform
    component represents the conditional model for the unmodifed state, and the
    Beta component represents the conditional model for the modified state. The
    shape parameters of Beta distribution are set to (1, 10) to approximately
    match the likelihood under the null hypothesis (unmodified state) for the
    LDRs in the top quintile of the distribution.

    We also provide an option to optimise the shape parameters of Beta
    distribution with the EM algorithm. To do this, make the function call with
    the last parameter \code{optimise} set to the desired tolerance, within
    which the parameter estimates are allowed not to change at the next
    iteration. The EM algorithm is run for 10 iterations and each M-step is
    performed with Newton's optimisation method and is allowed 20 iterations. We
    remark that during our experiments, this optimisation appeared vulnerable to
    local minima. Thus, the default version of the pipeline does not use this
    optimisation.
}
\value{
    An n-by-2 matrix of posterior probabilities for all selected nucleotides.
    First column corresponds to the probability of the nucleotide being
    unmodified (given the observations) and the second column - to the
    probability of the nucleotide being modified. Those nucleotides not
    belonging to stretches for which the HMM was run are assigned an `NA'.
}
\references{
    Selega et al. "Robust statistical modeling improves sensitivity of
    high-throughput RNA structure probing experiments", Nature Methods (2016).
}
\author{
    Alina Selega, Sander Granneman, Guido Sanguinetti
}
\section{Error }{
    The following errors are returned if:

    "Number of control and treatment replicates must be at least 2."
    the number of control or treatment experimental replicates is less than 2;

    "All lists of positions selected for pair-wise comparisons should
              be non-empty."
    a list of nucleotides for control-control or treatment-control comparisons
    is empty;

    "Strand should be either plus or minus, specified with a sign."
    strand is not specified as "+" or "-";

    "The list of considered nucleotide positions should be non-empty. If
    patterns are not used, a vector with all positions should be provided in the
    first element of the list."
    the list of considered nucleotides (which should either correspond to
    patterns if addressing sequence bias or list all positions) is empty;

    "The matrix of control-control (treatment-control) LDRs should have as many
    columns as there are control-control (treatment-control) comparisons."
    dimensionalities of matrices with LDRs do not agree with the number of
    pair-wise comparisons;

    "Please provide a tolerance if shape parameters are to be optimised with EM
    algorithm."
    the \code{optimise} parameter was not set to a (real) tolerance value.

}

\seealso{
    See Also \code{\link[BUMHMM]{stabiliseVariance}},
    \code{\link[BUMHMM]{selectNuclPos}}, \code{\link[BUMHMM]{findPatternPos}},
    and \code{\link[BUMHMM]{computeStretches}}.
}
\examples{
    library(SummarizedExperiment)

    Nc <- 3
    Nt <- 3
    t <- 1

    nuclSelection <- selectNuclPos(se, Nc, Nt, t)

    assay(se, "dropoff_rate") <- scaleDOR(se, nuclSelection, Nc, Nt)

    stretches <- computeStretches(se, t)

    varStab <- stabiliseVariance(se, nuclSelection, Nc, Nt)

    LDR_C <- varStab$LDR_C
    LDR_CT <- varStab$LDR_CT

    sequence <- subject(rowData(se)$nucl)
    nuclPosition <- list()
    nuclPosition[[1]] <- 1:nchar(sequence)

    posteriors <- computeProbs(LDR_C, LDR_CT, Nc, Nt, '+', nuclPosition,
                             nuclSelection$analysedC, nuclSelection$analysedCT,
                             stretches)
}

\keyword{models}
\keyword{htest}
\keyword{classif}
