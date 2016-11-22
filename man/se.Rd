\docType{data}
\name{se}
\alias{se}
\title{
    Example RNA structure probing data set.
}
\format{An instance of the \code{\link{SummarizedExperiment}} class storing data
    obtained in a structure probing experiment for the yeast rRNA 18S using the
    DMS chemical probing agent. The data can be accessed in the Gene Expression
    Omnibus with the code GSE52878.

    The available data types are summarised in three assays: coverage, drop-off
    count of the reverse transcriptase, and drop-off rate at each nucleotide
    position. The coverage and drop-off count data are obtained in a
    random-primed paired-end sequencing structure probing experiment. The
    drop-off rate is computed as the ratio between the drop-off count and
    coverage at each nucleotide position. Control (no reagent added) and
    treatment experiments are available in triplicates.

    The data from each assay is summarised in a rectangular matrix. The rows
    represent nucleotide position and the columns represent experimental
    replicates. The replicates are labelled as "control" or "treatment", which
    is stored as column data and can be accessed through the field
    \code{replicate}. The row data stores the associated DNA sequence,
    represented with a \code{DNAString} object.}
\usage{
    se
}
\value{
    RNA structure probing data.
}
\description{
    A \code{\link{SummarizedExperiment}} container storing the genomic sequence
    of the yeast rRNA 18S and three matrices with coverage, drop-off count, and
    drop-off rate information for each nucleotide position in three replicates
    of control and treatment DMS structure probing experiments.
}
\details{
    The \code{\link{SummarizedExperiment}} object can be created from individual
    matrices and the sequence string as follows. Below:
    \describe{
        \item{\code{seq}}{a string storing the genomic sequence}
        \item{\code{covFile}}{a matrix storing coverage information, with rows
        as nucleotide positions and columns as experimental replicates (control
        replicates come first, followed by treatment replicates)}
        \item{\code{docFile}}{a matrix storing drop-off count information,
        arranged as \code{covFile}}
        \item{\code{dorFile}}{a matrix storing drop-off rate information,
        arranged as \code{covFile}}}

    The \code{seq}, \code{covFile}, \code{docFile}, and \code{dorFile}
    are also provided as accompanying data sets for reference.

    Example code:

    library(BUMHMM)
    library(Biostrings)
    library(SummarizedExperiment)

    dna = DNAString(seq)
    se <- SummarizedExperiment(
        list(
            coverage=as.matrix(covFile),
            dropoff_count=as.matrix(docFile),
            dropoff_rate=as.matrix(dorFile)
        ), colData=DataFrame(
            replicate=rep(c("control", "treatment"), each=3)
        ), rowData=DataFrame(
            nucl=Views(dna, successiveIRanges(rep(1, nchar(dna))))
        ))
    colnames(se) <- c('C1', 'C2', 'C3', 'T1', 'T2', 'T3')
}
\references{
    Hector, R. D. et al. "Snapshots of pre-rRNA structural flexibility reveal
    eukaryotic 40S assembly dynamics at nucleotide resolution." Nucleic acids
    research (2014).
}
\keyword{datasets}
