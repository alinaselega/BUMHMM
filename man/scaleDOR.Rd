\name{scaleDOR}
\alias{scaleDOR}
\title{
    Function to scale drop-off rates across replicates to a common median.
}
\description{
    The function extracts all nucleotides selected for control-control and
    treatment-control comparisons and scales them to have a common median across
    all replicates. This value is computed as the median of all drop-off rates
    of the selected nucleotides in control replicate samples. The function
    requires the output of \code{\link[BUMHMM]{selectNuclPos}}, which holds
    lists of nucleotide positions selected for pair-wise comparisons.
}
\usage{
    scaleDOR(se, nuclSelection, Nc, Nt)
}
\arguments{
    \item{se}{
        A \code{\link{SummarizedExperiment}} object storing structure probing
        data and the associated genomic sequence. The documentation for the
        example data set provided with the package \code{\link[BUMHMM]{se}}
        outlines how the object should be defined. \code{scaleDOR} uses the
        assay "dropoff_rate".
    }
    \item{nuclSelection}{
        A list returned by \code{\link[BUMHMM]{selectNuclPos}}, containing the
        positions of nucleotides selected for all control-control and
        treatment-control comparisons.
    }
    \item{Nc}{
        Number of control replicate samples. Must be at least 2.
    }
    \item{Nt}{
        Number of treatment replicate samples. Must be at least 2.
    }
}
\value{
    Returns a modified \code{n}-by-\code{m} matrix of drop-off rates, scaled
    per replicate to have the same median value, where \code{n} is the number of
    nucleotides and \code{m} is the total number of replicate samples.
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

    "Nucleotide selection should have two elements."
    "All lists of positions selected for pair-wise comparisons or for which
    posteriors will be computed should be non-empty."
    the list containing positions of nucleotides selected for control-control
    and treatment-control comparisons does not have 2 elements or any of the
    elements is empty;

    "Drop-off rate matrix should not have any NA entries."
    the drop-off rate matrix has NA entries.
}

\seealso{
    See Also \code{\link[BUMHMM]{selectNuclPos}}.
}
\examples{
    Nc <- 3
    Nt <- 3
    t <- 1
    nuclSelection <- selectNuclPos(se, Nc, Nt, t)
    dorFile <- scaleDOR(se, nuclSelection, Nc, Nt)

}
\keyword{manip}
\keyword{univar}
