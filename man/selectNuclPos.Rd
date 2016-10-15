\name{selectNuclPos}
\alias{selectNuclPos}
\title{
    Function to compute selections of nucleotide positions.
}
\description{
    This function selects pairs of nucleotide positions for computing log-ratios
    of control-control and treatment-control replicate comparisons, and the
    nucleotide positions for which the posterior probabilities of modification
    will be computed.
}
\usage{
    selectNuclPos(covFile, docFile, Nc, Nt, t)
}
\arguments{
    \item{covFile}{
        A file with coverage information at each nucleotide position (rows).
        Experimental replicates should be arranged in different columns. Control
        replicates come first, followed by treatment replicates.
    }
    \item{docFile}{
        A file of drop-off counts at each nucleotide position. Same arrangement
        as for \code{covFile}.
    }
    \item{Nc}{
        Number of control experimental replicates. Must be at least 2.
    }
    \item{Nt}{
        Number of treatment experimental replicates. Must be at least 2.
    }
    \item{t}{
        Threshold for the minimum allowable coverage. Must be non-negative.
    }
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
    \item{computePosteriors }{
        List where each element corresponds to a treatment replicate. Each
        element holds indices of nucleotides that have coverage >= \code{t} in
        all replicates and a drop-off count > 0 in that experimental replicate.
    }
}
\references{
    Selega et al. "Robust statistical modeling improves sensitivity of
    high-throughput RNA structure probing experiments" (in press)
}
\author{
    Alina Selega, Sander Granneman, Guido Sanguinetti
}
\section{Error }{
    The following rrrors are returned if:

    "The number of experimental replicates must be at least 2."
    the number of control or treatment experimental
    replicates is less than 2;

    "The minumum coverage threshold must be non-negative."
    the threshold for the minimum considered coverage is negative;

    "The coverage and drop-off count matrices should not have NA entries."
    the coverage and drop-off count matrices have NA entries.
}
\examples{
    selectNuclPos(covFile, docFile, 3, 3, 1)
}
\keyword{manip}
