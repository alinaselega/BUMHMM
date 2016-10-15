\docType{data}
\name{covFile}
\alias{covFile}
\title{
    Example coverage data set.
}
\format{A matrix containing the coverage information for each nucleotide
    position of a ribosomal RNA 18S obtained in a structure probing experiment
    using DMS as a probe, where:
    \describe{
        \item{1st column}{control experimental replicate}
        \item{2nd column}{control experimental replicate}
        \item{3rd column}{control experimental replicate}
        \item{4th column}{treatment experimental replicate}
        \item{5th column}{treatment experimental replicate}
        \item{6th column}{treatment experimental replicate}
}}
\usage{
    covFile
}
\value{
    Coverage information per nucleotide.
}
\description{
    A matrix containing 1800x6 entries of per nucleotide coverage information
    from structure probing experiments.
}
\references{
    Hector, R. D. et al. "Snapshots of pre-rRNA structural flexibility reveal
    eukaryotic 40S assembly dynamics at nucleotide resolution." Nucleic acids
    research (2014).
}
\keyword{datasets}
