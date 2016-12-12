\name{findPatternPos}
\alias{findPatternPos}
\title{
    Function to find positions of the nucleotide patterns in the sequence.
}
\description{
    This function finds all occurrences of a nucleotide pattern in the sequence.
    For each occurrence, the function returns the index of the middle
    nucleotide, computed as: \code{ceiling(length(pattern) / 2)}. The function
    supports data for the plus and minus DNA strands; for the minus strand, all
    patterns are turned to complementary sequence.
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
        A \code{\link{DNAString}} object storing the reference genomic sequence
        to search for the patterns in. The sequence corresponding to plus strand
        is expected.
    }
    \item{strand}{
        A character, indicating the plus (\code{+}) or minus strand (\code{-}).
        For the minus strand, the occurrences found for a particular pattern
        will be attributed to the pattern with complementary sequence.
    }
}
\details{
    This function uses \code{\link[stringi]{stri_locate_all_fixed}}.

    This function aims to assist with addressing sequence bias in structure
    probing data. The sequence in the neighbourhood of a nucleotide is assumed
    to have an effect on its structural state. By considering sequence patterns
    of a certain length (specified by the user), this function finds indices of
    the middle nucleotide of each pattern's occurrences within the sequence. We
    then separately analyse the nucleotides occurring in the middle of each
    pattern, taking into account sequence dependency.
}
\value{
    This function returns a list where each component corresponds to a pattern
    (indicated by the field \code{names}) and contains indices of the middle
    nucleotides of that pattern's occurrences within the sequence.
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
    library(SummarizedExperiment)

    ## Extract the DNA sequence from se
    sequence <- subject(rowData(se)$nucl)

    ## Generate patterns of length 3
    n <- 3
    patterns <- nuclPerm(n)

    ## Find positions of pattern occurrences
    nuclPosition <- findPatternPos(patterns, sequence, '+')
}

\keyword{classes&character}
\keyword{manip}
