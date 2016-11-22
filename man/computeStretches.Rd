\name{computeStretches}
\alias{computeStretches}
\title{
    Function to compute continuous stretches of nucleotide positions.
}
\description{
    This function computes continuous stretches of nucleotide positions, for
    which posterior probabilities are to be computed. Such positions should have
    the minimum allowed coverage (as defined with the parameter \code{t}) in
    all experimental replicates and a non-zero drop-off count in at least one of
    the treatment replicates. The returned stretches are at least two
    nucleotides long.
}
\usage{
    computeStretches(se, t)
}
\arguments{
    \item{se}{
        A \code{\link{SummarizedExperiment}} object storing structure probing
        data and the associated genomic sequence. The documentation for the
        example data set provided with the package \code{\link[BUMHMM]{se}}
        outlines how the object should be defined. \code{computeStretches} uses
        the assays "coverage" and "dropoff_count".
    }
    \item{t}{
        Threshold for the minimum allowed coverage. Must be non-negative.
    }
}
\value{
    An \code{\link{IRanges}} object storing each stretch.
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

    "The minumum coverage threshold must be non-negative."
    the threshold for the minimum considered coverage is negative.
}
\examples{
    t <- 1
    stretches <- computeStretches(se, t)
}
\keyword{manip}
\keyword{arith}
