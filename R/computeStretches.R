computeStretches <- function(indices) {

    if (length(indices) < 2) {
        stop('There should be at least 2 nucleotide positions in the list.')
    }
    else {

      if (length(indices) > 10000) {
          message('Computing stretches... This might take a few minutes.')
      }

      ## Compute continuous stretches of nucleotides for which posterior
      ## posterior probabilities will be computed
      stretches <- list()
      k <- 1
      stretchStart <- 1

      for (i in 2:length(indices)) {
          if ((indices[i] - indices[i-1]) > 1) {
              if (stretchStart < i - 1) {
                  stretches[[k]] <- c(indices[stretchStart],
                                      indices[i - 1])
                  k <- k + 1
              }
              stretchStart <- i
          }
      }

      if (stretchStart != length(indices)) {
          stretches[[k]] <- c(indices[stretchStart],
                              indices[length(indices)])
      }

      return(stretches)
    }
}
