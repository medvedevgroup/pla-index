# Read Aligner Application
To use the read aligner application, run the following command:
```shell
git checkout strobealign-application
```

# Exact Version
To use `PLA-index-exact` version, check out to the `pla-index-exact` branch, and follow the `README` instructions there.
```shell
git checkout pla-index-exact
```

# Reproducibility Information
`Refseq_dataset.csv` file contains the RefSeq IDs, Kingdom, alpha and beta values of `Table 5` from our paper. 
`RefseqID_7Genomes.pdf` contains the RefSeq IDs of the genomes from `Table 6`.
`RefSeq_ID_List.txt` contains the RefSeq IDs of the 549 genomes that were used for analyzing the pla-complexity.

`Snakemake` pipeline can be used to generate the results shown in our paper.
To do this, run `cmake` as before, as this will genereate the necessary executables to process the data corresponding to the `Snakemake` input.

# Requirements
In addition to the requirements already stated, one will also need the followings:
- [Snakemake](https://snakemake.readthedocs.io/en/stable/) (for automating all computations)
- R (for finding the alpha and beta value (from pla-complexity) from segments vs eps curve fitting)

## Auxilary Files
To build and query pla-index, we need the whole genome 'N' filtered and concatenated, the binary suffix array file and the query files as before. 
To produce these files, the executables files inside the `pla-index/build` folder can be used.

To filter 'N' from a genome and concatenate the genome to a string:
```
./process_fasta FASTA-FILE OUTPUT-PREFIX
```
Parameter description:

| Parameter Name | Description |
|----------|----------|
| FASTA-FILE | Fasta file on which pla-index will be built |
| OUTPUT-PREFIX | Prefix of the output file. The file name will be `OUTPUT-PREFIX.processed.fasta` |

The output fasta will have a single header line and a single string containing all the bases (without 'N') concatenated.

To construct suffix array on the concatenated genome:

```
./mksary PROCESSED-GENOME-FASTA SUFFIX-ARRAY-FILE
```
Parameter description:
| Parameter Name | Description |
|----------|----------|
| PROCESSED-GENOME-FASTA | Fasta file name of the concatenated genome without any 'N' character at the bases (one single string). Essentially, this file has only two lines: 1) Header line 2) All bases concatenated single string |
| SUFFIX-ARRAY-FILE | Name of the binary file where suffix array will be stored |

We use the suffix array format of `mksary` executable output. This is basically a binary file of the suffix array (each index value being a 64 bit integer).

To construct query files:
```
./create_queries PROCESSED-GENOME-FILE SUFFIX-ARRAY-FILE PREFIX KMER-SIZE NO-OF-QUERIES NO-OF-FILES
```

Parameter description:

Parameter description:
| Parameter Name | Description |
|----------|----------|
|  PROCESSED-GENOME-FASTA | Fasta file name of the concatenated genome without any 'N' character at the bases (one single string). Essentially, this file has only two lines: 1) Header line 2) All bases concatenated single string |
| SUFFIX-ARRAY-FILE | Name of the binary file where suffix array is stored |
| PREFIX | Prefix of the query file name |
| KMER-SIZE | Kmer size of the query kmers |
| NO-OF-QUERIES | How many query k-mers to randomly construct |
| NO-OF-FILES | How many query files to construct |

Thus, for example, if we have an `ecoli.fasta` file, we can remove 'N' and concatenate the genome by:
```shell
./process_fasta ecoli.fasta ecoli
```
The output file is written at `ecoli.processed.fasta`

Next, we can build the suffix array of this genome by:
```shell
./mksary ecoli.processed.fasta ecoli.sa.bin
```
The output suffix array is written at `ecoli.sa.bin`

Next, we can create one query file of 5,000,000 queries of 21 size k-mer for this genome by:
```shell
./create_queries ecoli.processed.fasta ecoli.sa.bin ecoli 21 5000000 1
```
The output query file is written at `ecoli.1.query.txt`

## Snakemake

To automate the index building and querying using snakemake, create a folder containing the genome fasta file, binary file of the suffix array and query files.
The snakemake file assumes the following naming conventions:
- Genome folder name without any white-spaces(example: GENOME)
- Fasta file with one entry and only ACGT characters inside the genome folder needs to be named as GENOME/GENOME.processed.fasta
- Binary suffix array file (64 bit integers as suffix array values) inside the genome folder as GENOME/GENOME.sa.bin
- Query files (one kmer at each line) with their tag number as GENOME/GENOME.QUERY_TAG.query.txt

Thus, if one have a fasta file at hand, output from `process_fasta`, `mksary` and `create_queries` (all are inside the `build` folder) can be used to format the files necessary for `snakemake`.

For example, let the genome folder name be `ecoli`. 
Then, we will have three files inside `ecoli` folder: `ecoli/ecoli.processed.fasta`, `ecoli/ecoli.sa.bin`, and `ecoli/ecoli.1.query.txt` (`1` being the tag number of the query)

Then, create a config file using CreateConfigFile.py. To see the config options:
```shell
python3 CreateConfigFile.py -h
```
Config options:
```
usage: CreateConfigFile.py [-h] --genome_folder GENOME_FOLDER --epsilon EPSILON [--index_type INDEX_TYPE]
                           [--fast_rank FAST_RANK] [--query_type QUERY_TYPE] [--kmer_size KMER_SIZE]
                           [--code_path CODE_PATH] [--exec_path EXEC_PATH] [--l_bits L_BITS] [--query_tag QUERY_TAG]
                           [--sdsl_lib_path SDSL_LIB_PATH] [--sdsl_inc_path SDSL_INC_PATH]

PLA-Index config file builder

options:
  -h, --help            show this help message and exit
  --genome_folder GENOME_FOLDER
                        Name of the folder containing the genome, suffix array and query files [Required]
  --epsilon EPSILON     Epsilon value to be used to create pla-index [Required]
  --index_type INDEX_TYPE
                        Whether to build "basic-pla" or "repeat-pla" index. "repeat-pla" is generally the better option
                        but can slow down rank queries unless the -r option is also used.(default: basic-pla)
  --fast_rank FAST_RANK
                        Construct an extra bitvector with the same length as the suffix array. Use this to speed up rank
                        queries. Possible values: "y" for yes and "n" for no. (default: n)
  --query_type QUERY_TYPE
                        Whether to do "search" or "rank" query. "rank" gives you the value of the first position in the
                        suffix array where the query is found. "search" on the other hand, also returns the value where
                        the query appears, but it might not be the very first spot. (default: search)
  --kmer_size KMER_SIZE
                        Kmer size for which pla-index will be calculated (default: 21)
  --code_path CODE_PATH
                        Path where the source code is (default: ../src/)
  --exec_path EXEC_PATH
                        Path where the executables will be stored (default: ../executables/)
  --l_bits L_BITS       How many elements to store in the shortcut array, D. |D| = 2^l (default: 16)
  --query_tag QUERY_TAG
                        QUERY_TAG to identify which query file to use from GENOME.QUERY_TAG.query.txt (default: 1)
  --sdsl_lib_path SDSL_LIB_PATH
                        Path to the SDSL library folder (default: ~/lib)
  --sdsl_inc_path SDSL_INC_PATH
                        Path to the SDSL include folder (defaul: ~/include)
```

Here is an example using the `ecoli` folder inside the `tests` folder:

```
cd tests
mkdir -p ../executables
python3 ../CreateConfigFile.py --genome_folder ecoli --epsilon 15
snakemake -s ../Snakefile --cores=1 -p
```

## Fitting Segments

To find the fitted degree (alpha) and fitted constant (beta) values, `fit_segments.R` can be used.
First, the epsilon values vs number of segments needs to be listed in a text file. (one such example is `tests/Gorilla.segments.txt`).
Then, list this file path and the genome length (each at a different line) in `fit_segment_file.txt`. (similar to `Reproducibility/fit_segment_file.txt`).
Then, running `fit_segments.R` will output the goodness of fits as well as the alpha and beta values.
