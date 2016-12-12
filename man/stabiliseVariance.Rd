\name{stabiliseVariance}
\alias{stabiliseVariance}
\title{
    Function to reduce variance of LDRs as a function of coverage.
}
\description{
     This function computes the log drop-off rate ratios (LDRs) for
     control-control and treatment-control comparisons and applies a
     transformation to them that reduces their variance as a function of
     coverage.
}
\usage{
    stabiliseVariance(se, nuclSelection, Nc, Nt)
}
\arguments{
  \item{se}{
        A \code{\link{SummarizedExperiment}} object storing structure probing
        data and the associated genomic sequence. The documentation for the
        example data set provided with the package \code{\link[BUMHMM]{se}}
        outlines how the object should be defined. \code{stabiliseVariance} uses
        the assays "coverage" and "dropoff_rate".
  }
  \item{nuclSelection}{
        A list returned by \code{\link[BUMHMM]{selectNuclPos}}, containing the
        positions of nucleotides selected for all control-control and
        treatment-control comparisons.
  }
  \item{Nc}{
      Number of control experimental replicates. Must be at least 2.
  }
  \item{Nt}{
      Number of treatment experimental replicates. Must be at least 2.
  }
}
\details{
    The variance is reduced by sorting all LDRs in the null distribution by the
    average coverage and splitting the data in bins. Each bin spans a coverage
    of 100; or, if maximum coverage is not larger by 100 than the minimum
    coverage, the range is set to their difference divided by 10.

    For each bin, the 95th quantile of LDRs with subtracted mean and the average
    coverage are computed. Then non-linear least squares are used to fit the
    following model (with parameters k, b):
    f = 1/sqrt(n) * k + b,
    for f - quantiles, n - mean coverage in the bin.

    All LDRs are then rescaled by this model according to their corresponding
    average coverage in the pair of replicates.
}
\value{

    \item{LDR_C }{
        A matrix of transformed LDRs for control-control comparisons. The matrix
        rows correspond to nucleotide positions and columns to a control-control
        comparison. Only those positions selected for a pair-wise comparison
        will be assigned a value; the rest will be left as an \code{NA}.
    }
    \item{LDR_CT }{
        A matrix of transformed LDRs for treatment-control comparisons. The
        matrix rows correspond to nucleotide positions and columns to a
        treatment-control comparison. Only those positions selected for a
        pair-wise comparison will be assigned a value; the rest will be left as
        an \code{NA}.
    }
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

    "All lists of positions selected for pair-wise comparisons should be
    non-empty."
    a list of nucleotides for control-control or treatment-control comparisons
    is empty;

    "The coverage and drop-off count matrices should not have NA entries."
    the coverage and drop-off count matrices have NA entries;

    "Unable to fit the model for correcting the coverage bias."
    The function \code{\link{nls}} could not execute successfully. The 95th
    quantiles of the LDR distribution in each bin could be equal to 0 or not
    have enough elements. This would happen if not enough nucleotides end up in
    a bin; e.g. one nucleotide per bin.
}

\seealso{
    See Also \code{\link[BUMHMM]{selectNuclPos}}.
}
\examples{
    library(SummarizedExperiment)
    Nc <- 3
    Nt <- 3
    t <- 1
    nuclSelection <- selectNuclPos(se, Nc, Nt, t)
    assay(se, "dropoff_rate") <- scaleDOR(se, nuclSelection, Nc, Nt)
    varStab <- stabiliseVariance(se, nuclSelection, Nc, Nt)
}

\keyword{models&nonlinear}
