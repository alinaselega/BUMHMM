## This is a function for internal use by scaleDOR.R
poolNucl <- function(c, comparisons, indices) {

    ## Pool nucleotide positions selected for comparisons including replicate c
    rep <- which(comparisons == c, arr.ind=TRUE)[, 1]
    return(Reduce(union, indices[rep]))
}
