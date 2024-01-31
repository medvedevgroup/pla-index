# Reproducibility Information
`Snakemake` pipeline can be used to generate the results shown in our paper.
To do this, run `cmake` as before, as this will genereate the necessary executables to process the data corresponding to the `Snakemake` input.

## Auxilary Files
To build and query pla-index, we need the whole genome 'N' filtered and concatenated, the binary suffix array file and the query files as before. 
To produce these files, the executables files inside the `pla-index/build` folder can be used.

To filter 'N' from a genome, and concatenate afterwards:
```
./filterNConcat FASTA-FILE OUTPUT-PREFIX

Parameter description:
FASTA-FILE: [string], Fasta file on which pla-index will be built
OUTPUT-PREFIX: [string], The output file name format is OUTPUT-PREFIX.ConcatenatedGenome.txt
```

To construct suffix array on the concatenated genome:

```
./mksary GENOME SUFFIX-ARRAY

Parameter description:
GENOME: [string], Concatenated genome without 'N' (one single string)
SUFFIX-ARRAY: [string], Name of the binary file where suffix array will be stored
```

We use the suffix array format of `mksary` executable output. This is basically a binary file of the suffix array(each index value being a 64 bit integer).

To construct query files:
```
./create_queries GENOME SUFFIX-ARRAY PREFIX KMER-SIZE NO-OF-QUERIES NO-OF-FILES

Parameter description:

GENOME: [string], Whole genome string concatenated filtered by 'N'
SUFFIX-ARRAY: [string], Suffix array for GENOME in a binary file
PREFIX: [string], File name prefix of the query file
KMER-SIZE: [int], Kmer size to be used to construct the index
NO-OF-QUERIES: [int], How many query k-mers to randomly construct
NO-OF-FILES: [int], How many query files to construct
```

Thus, for example, if we have an `ecoli.fasta` file, we can remove 'N' and concatenate the genome by:
```shell
./filterNConcat ecoli.fasta ecoli
```
The output file is written at `ecoli.ConcatenatedGenome.txt`

Next, we can build the suffix array of this genome by:
```shell
./mksary ecoli.ConcatenatedGenome.txt ecoli.sa.bin
```
The output suffix array is written at `ecoli.sa.bin`

Next, we can create one query file of 5,000,000 queries of 21 size k-mer for this genome by:
```shell
./create_queries ecoli.ConcatenatedGenome.txt ecoli.sa.bin ecoli 21 5000000 1
```
The output query file is written at `ecoli.1.query.txt`

## Snakemake

To automate the index building and querying using snakemake, create a folder containing the genome, binary file of the suffix array and query files.
The snakemake file assumes the following naming conventions:
- Genome folder name without any white-spaces(example: GENOME)
- Concatenated genome without N inside the genome folder needs to be named as GENOME/GENOME.ConcatenatedGenome.txt
- Binary suffix array file (64 bit integers as suffix array values) inside the genome folder as GENOME/GENOME.sa.bin
- Query files (one kmer at each line) with their tag number as GENOME/GENOME.QUERY_TAG.query.txt

Thus, if one have a fasta file at hand, output from `filterNConcat`, `mksary` and `create_queries` can be used to format the files necessary for `snakemake`.

For example, let the genome folder name be `ecoli`. 
Then, we will have three files inside `ecoli` folder: `ecoli/ecoli.ConcatenatedGenome.txt`, `ecoli/ecoli.sa.bin`, and `ecoli/ecoli.1.query.txt` (`1` being the tag number of the query)

Then, create a config file using CreateConfigFile.py

Required Parameters:

| Parameter  | Type    | Description    |
|-------------|-------------|-------------|
|--genome_folder | [string] |The name of the folder containing the genome, suffix array and query files|
|--epsilon |  [int]   |Epsilon value to be used to create pla-index|

Optional parameters with required argument:

| Parameter  | Type    | Description    |
|-----------------|-------------|-------------|
|--index_type     |[string] | What kind of index type to construct and/or use (default: basic-pla). Possible values: "basic-pla" or "repeat-pla"|
|--query_type     |[string] | What kind of query to do (default: search). Possible values: "search" or "rank"|
|--kmer_size      |[int] | Kmer size for which pla-index will be calculated (default: 21)|
|--code_path      |[string] | Path where the source code is (default: ../src/)|
|--exec_path      |[string] | Path where the executables will be stored (default: ../executables/)|
|--l_bits         |[int]  | How many elements to store in the shortcut array, D. &#124;D&#124; = 2<sup>l</sup> (default: 16)|
|--query_tag      |[int] | Which query file to use (default: 1)|
|--sdsl_lib_path  |[string] | Path to the SDSL library folder (default: ~/lib)|
|--sdsl_inc_path  |[string] | Path to the SDSL include folder (defaul: ~/include)|

Here is an example assuming we have the `ecoli` folder inside the `tests` folder:

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
