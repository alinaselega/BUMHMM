\name{nuclPerm}
\alias{nuclPerm}
\title{
    Function to create nucleobase patterns.
}
\description{
    This function creates all permutations of 4 nucleobases (A, T, G, C) of
    length \code{n}.
}
\usage{
    nuclPerm(n)
}
\arguments{
    \item{n}{
        Length of the pattern.
    }
}
\details{
    This function uses \code{\link[gtools]{permutations}}.
}
\value{
    A vector of characters with each element being a nucleobase pattern of
    length \code{n}.
}
\references{
    Selega et al. "Robust statistical modeling improves sensitivity of
    high-throughput RNA structure probing experiments" (in press)
}
\author{
    Alina Selega, Sander Granneman, Guido Sanguinetti
}
\section{Error }{
    The following error is returned if:

    "The length of patterns provided is not a positive number."
    the provided length of patterns to be generated is not positive.
}

\seealso{
    \code{\link[gtools]{permutations}}
}
\examples{
    nuclPerm(3)
}
\keyword{character}
