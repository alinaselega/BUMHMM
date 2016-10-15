\name{computeStretches}
\alias{computeStretches}
\title{
    Function to compute continuous stretches of nucleotide positions.
}
\description{
    This function computes continuous stretches from the list of nucleotide
    position indices, for which posterior probabilities are to be computed,
    returned in the third element of the output of the
    \code{\link[BUMHMM]{selectNuclPos}} function. The stretch is at least two
    nucleotides long.
}
\usage{
    computeStretches(indices)
}
\arguments{
    \item{indices}{
        A sorted list of nucleotide positions for which posterior probabilities
        are to be computed. This list can be obtained from the third element of
        the list returned by the \code{\link[BUMHMM]{selectNuclPos}} function.
    }
}
\details{
    For example, if the posterior probabilities are to be computed for
    nucleotides with positions 1, 2, 3, 4, 7, 8, 9, 15, 17, 18, the function
    will return three stretches: [1 4], [7 9], and [17 18].
}
\value{
    A list where each element corresponds to one stretch and has two elements:
    the start index and the end index of the stretch.
}
\references{
    Selega et al. "Robust statistical modeling improves sensitivity of
    high-throughput RNA structure probing experiments" (in press)
}
\author{
    Alina Selega, Sander Granneman, Guido Sanguinetti
}
\section{Error }{
    The following errors are returned if:

    "There should be at least 2 nucleotide positions in the list."
    less than two position indices are provided in the list.
}

\seealso{
    See Also \code{\link[BUMHMM]{selectNuclPos}}.
}
\examples{
    nuclSelection <- selectNuclPos(covFile, docFile, 3, 3, 1)

    ## Pool together all positions in a single list and sort it
    nuclForPosteriors <- sort(Reduce(union, nuclSelection$computePosteriors))
    stretches <- computeStretches(nuclForPosteriors)
}
\keyword{manip}
\keyword{arith}
