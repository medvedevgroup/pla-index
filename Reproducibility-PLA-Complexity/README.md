# Reproducibility Information

`Refseq_dataset.csv` file contains the RefSeq IDs, Kingdom, alpha and beta values of `Table 5` from our paper.
`RefseqID_7Genomes.pdf` contains the RefSeq IDs of the genomes from `Table 6`.
`RefSeq_ID_List.txt` contains the RefSeq IDs of the 549 genomes that were used for analyzing the pla-complexity.

## Requirements

- R (for finding the alpha and beta value (from pla-complexity) from segments vs eps curve fitting)

## Fitting Segments

To find the fitted degree (alpha) and fitted constant (beta) values, `fit_segments.R` can be used.
First, the epsilon values vs number of segments needs to be listed in a text file. (one such example is `tests/Gorilla.segments.txt`).
Then, list this file path and the genome length (each at a different line) in `fit_segment_file.txt`. (similar to `Reproducibility/fit_segment_file.txt`).
Then, running `fit_segments.R` will output the goodness of fits as well as the alpha and beta values.