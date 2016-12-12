\name{selectNuclPos}
\alias{selectNuclPos}
\title{
    Function to compute selections of nucleotide positions.
}
\description{
    This function selects pairs of nucleotide positions for computing log-ratios
    of control-control and treatment-control replicate comparisons.
}
\usage{
    selectNuclPos(se, Nc, Nt, t)
}
\arguments{
    \item{se}{
        A \code{\link{SummarizedExperiment}} object storing structure probing
        data and the associated genomic sequence. The documentation for the
        example data set provided with the package \code{\link[BUMHMM]{se}}
        outlines how the object should be defined. \code{selectNuclPos} uses the
        assays "coverage" and "dropoff_count".
    }
    \item{Nc}{
        Number of control experimental replicates. Must be at least 2.
    }
    \item{Nt}{
        Number of treatment experimental replicates. Must be at least 2.
    }
    \item{t}{
        Threshold for the minimum allowed coverage. Must be non-negative.
    }
}
\details{
    This function uses \code{\link[utils]{combn}}.
}
\value{
    \item{analysedC }{
        List where each element corresponds to a control-control replicate
        comparison. Each element holds indices of nucleotides that have coverage
        >= \code{t} and a drop-off count > 0 in both replicates of that
        comparison.
    }
    \item{analysedCT }{
        List where each element corresponds to a treatment-control replicate
        comparison. Each element holds indices of nucleotides that have coverage
        >= \code{t} and a drop-off count > 0 in both replicates of that
        comparison.
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

    "The number of experimental replicates must be at least 2."
    the number of control or treatment experimental replicates is less than 2;

    "The minumum coverage threshold must be non-negative."
    the threshold for the minimum considered coverage is negative;

    "The coverage and drop-off count matrices should not have NA entries."
    the coverage and drop-off count matrices have NA entries.
}
\examples{
    selectNuclPos(se, 3, 3, 1)
}
\keyword{manip}
