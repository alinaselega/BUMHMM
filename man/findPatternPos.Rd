\name{findPatternPos}
\alias{findPatternPos}
\title{
    Function to find positions of the nucleotide patterns in the sequence.
}
\description{
    This function finds all occurences of a nucleotide pattern in the sequence.
    The positions of the start and end of the occurrence are returned for each
    pattern. The function supports data for the plus and minus DNA strands;
    for the minus strand, all patterns are turned to complementary sequence.
}
\usage{
    findPatternPos(patterns, sequence, strand)
}
\arguments{
    \item{patterns}{
        A list of nucleotide permutations of length \code{n}, as returned by
        \code{\link[BUMHMM]{nuclPerm}}.
    }
    \item{sequence}{
        A string with the sequence to search for the patterns in. The sequence
        corresponding to plus strand is expected.
    }
    \item{strand}{
        A character, indicating the plus (\code{+}) or minus strand (\code{-}).
        For the minus strand, the occurrences found for a particular pattern
        will be attributed to the pattern with complementary sequence.
    }
}
\details{
    This function uses \code{\link[stringi]{stri_locate_all_fixed}}.
}
\value{
    This function returns a list where each component corresponds to a pattern
    (indicated by the field \code{names}) and contains the start and end
    positions of the pattern's occurrences within the sequence.
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

    "Strand should be either plus or minus, specified with a sign."
    strand is not specified as "+" or "-";

    "The sequence should be non-empty."
    provided sequence is empty;

    "The list of patterns should be non-empty."
    the list of patterns to search for in the sequence is empty.
}

\seealso{
    See also \code{\link[BUMHMM]{nuclPerm}}.
}
\examples{
    sequence <- 'TGACGTGGAGAGAGATTTTTCCCCCCTTAAAGC'
    patterns <- nuclPerm(3)
    nuclPosition <- findPatternPos(patterns, sequence, '+')
}

\keyword{classes&character}
\keyword{manip}
