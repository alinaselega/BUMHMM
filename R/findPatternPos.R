findPatternPos <- function(patterns, sequence, strand) {

    if ((strand != '+') & (strand != '-')) {
        stop("Strand should be either plus or minus, specified with a sign.")
    }
    else if (nchar(sequence) == 0) {
        stop('The sequence should be non-empty.')
    }
    else if (length(patterns) == 0) {
        stop('The list of patterns should be non-empty.')
    }
    else {

        indices <- list()

        ## Find all occurences of patterns in the full sequence
        for (i in 1:length(patterns)) {
            p <- patterns[i]
            indices[i] <- stringi::stri_locate_all_fixed(sequence, p,
                                                         overlap=TRUE)
        }

        ## Annotate found indices by the corresponding pattern
        if (strand == '+') {
            names(indices) <- patterns
        }

        ## Turn patterns to complementary if working with the minus strand
        if (strand == "-") {
            complPatterns = list()
            for (i in 1:length(patterns)) {
                chars = array('')
                for (j in 1:nchar(patterns[i])) {
                    chars[j] <- switch(substr(patterns[i], j, j),
                                'A' = 'T', 'T' = 'A', 'G' = 'C', 'C' = 'G')
                }
                complPatterns[[i]] <- paste(chars, collapse='')
            }
            ## Annotate found indices by the corresponding pattern
            names(indices) <- complPatterns
        }

        ## Return them in the order of original patterns
        position <- indices[patterns]

        return(position)
    }
}
