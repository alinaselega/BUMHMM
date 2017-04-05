# About

BUM-HMM (Beta-Uniform Mixture Hidden Markov Model) is a probabilistic modelling pipeline for computing per-nucleotide posterior probabilities of modification from structure probing data. The model supports multiple experimental replicates and empirically corrects coverage- and sequence-dependent biases. The model utilises a measure of "drop-off rate" for each nucleotide, which is compared between replicates through a log-ratio (LDR). The LDRs between control replicates define a null distribution of variability in drop-off rate observed by chance and LDRs between treatment and control replicates gets compared to this distribution. Resulting empirical p-values (probability of being "drawn" from the null distribution) are used as observations in a Hidden Markov Model with a Beta-Uniform Mixture model used as an emission model. The resulting posterior probabilities indicate the probability of a nucleotide of having being modified in a structure probing experiment.

# Software

BUM-HMM is available as a Bioconductor package at the following URL:

https://bioconductor.org/packages/BUMHMM/

# Contact

If you have any questions, please contact me at alina.selega@gmail.com.

# References

Selega, A., Sirocchi, C., Iosub, I., Granneman, S., & Sanguinetti, G. (2016). Robust statistical modeling improves sensitivity of high-throughput RNA structure probing experiments. *Nature Methods*. doi:10.1038/nmeth.4068.
