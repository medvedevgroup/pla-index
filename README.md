# pla-index-exact
pla-index-exact allows faster query of a k-mer rank function by storing errors from piece-wise linear approximation using an MPHF (pthash).

## Requirements
- Linux (64 bit)
- C++17
- [SDSL](https://github.com/simongog/sdsl-lite/tree/master) (>=2.1.1)
- cmake (>= 3.16)
- [Snakemake](https://snakemake.readthedocs.io/en/stable/) (for automating all computations)

## Quick Start

### Installation

```shell
git clone --recursive https://github.com/medvedevgroup/pla-index.git
cd pla-index
git checkout pla-index-exact
./update_cmake_lists.sh SDSL_INCLUDE_PATH SDSL_LIB_PATH
mkdir build
cd build
cmake ..
make -j 4
```

The two parameters for `update_cmake_lists.sh` are
| Parameter Name | Description |
|----------|----------|
| SDSL_INCLUDE_PATH  | Path to the SDSL include folder [typically ~/include/]   |
| SDSL_LIB_PATH  | Path to the SDSL library folder [typically ~/lib/]  |

In case you want to construct your own suffix array, you may use [libdivsufsort](https://github.com/hasin-abrar/libdivsufsort) tool. We modify it slightly to allow 64 bit integers.

To create the suffix array builder executable mksary from the current `build` folder:

```shell
cd ../libdivsufsort
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE="Release" \
-DCMAKE_INSTALL_PREFIX="/usr/local" ..
sed -i 's/int32_t/int64_t/g' include/divsufsort.h
make
cp examples/mksary ../../build/
```

### Usage

Now, all the executables are inside the `pla-index/build` folder.

The followings steps are almost similar to the `main` branch.
The only difference is that one does not need to specify the INDEX-TYPE or QUERY-TYPE anywhere. (It will always be `exact-pla` and `rank` operation)

```shell
./build_pla_index -h
```

would display the command line interface:

```text
 build_pla_index {OPTIONS}

    build_pla_index

OPTIONS:

    -h, --help    Print help and exit
    -g [STRING], --genome_fasta=[STRING]
                Fasta file with one entry and only ACGT characters. [Required]
    -s [STRING], --suffix_array=[STRING]
                Suffix array of the genome in a binary file. [Required]
    -k [INT], --kmer_size=[INT]
                Kmer size to be used to construct the index. [default: 21]
    -e [INT], --eps=[INT]
                Epsilon value to be used for constructing the pla-index. [default: 15]
    -o [STRING], --index=[STRING]
                File name where to save the index. [default: genome_fasta.index]
    -l [INT], --lookup=[INT]
                On average on how many elements the binary search on X array will take
                place. Used to determine the prefix lookup table. [default: 16]
```

Similarly, to see the query command line interface:

```shell
./query_pla_index -h
```

and the interface:

```text
query_pla_index {OPTIONS}

    query_pla_index

OPTIONS:

    -h, --help    Print help and exit
    -g [STRING], --genome_fasta=[STRING]
                Fasta file with one entry and only ACGT characters [Required]
    -s [STRING], --suffix_array=[STRING]
                Suffix array of the genome in a binary file [Required]
    -i [STRING], --index=[STRING]
                Name of the stored index file [Required]
    -q [STRING], --query=[STRING]
                Name of a text file containing queries. Each line of the file should
                contain a single k-mer. The length of the kmer should be the same as
                the one provided while building the index. [Required]
```

## Suffix Array

To build the suffix array:
```shell
./mksary INFILE OUTFILE
```
Parameter description:

| Parameter Name | Description |
|----------|----------|
| INFILE | Fasta file with one entry and only ACGT characters|
| OUTFILE  | Name of the file where the suffix array for INFILE will be stored |


## Small Example

We will use the formatted fasta file: `tests/ecoli/ecoli.processed.fasta` for the example.
From now on, we will be using the executables that are inside the pla-index/build folder.
To build suffix array:

```shell
./mksary ../tests/ecoli/ecoli.processed.fasta ../tests/ecoli/ecoli.sa.bin
```
The built suffix array is written on `tests/ecoli/ecoli.sa.bin`.

To build an `exact-pla` index with `21` size k-mer, default epsilon value (`15`), and default average number of elements (`16`) on which the binary search on X array will take place per query, we will use the formatted fasta file and the built suffix array:

```shell
./build_pla_index -g ../tests/ecoli/ecoli.processed.fasta -s ../tests/ecoli/ecoli.sa.bin -o ../tests/ecoli/ecoli.index
```

The built index will be saved on `../tests/ecoli/ecoli.index` file.

To do `rank` query with the provided `tests/ecoli/ecoli.1.query.txt` file on the built `tests/ecoli/ecoli.index` one can do the following:

```shell
./query_pla_index -g ../tests/ecoli/ecoli.processed.fasta -s ../tests/ecoli/ecoli.sa.bin -i ../tests/ecoli/ecoli.index -q ../tests/ecoli/ecoli.1.query.txt
```

The query time will be shown on console.

## Snakemake

To automate the index building and querying using snakemake, create a folder containing the genome fasta file, binary file of the suffix array and query files.
The snakemake file assumes the following naming conventions:

- Genome folder name without any white-spaces(example: GENOME)
- Fasta file with one entry and only ACGT characters inside the genome folder needs to be named as GENOME/GENOME.processed.fasta
- Binary suffix array file (64 bit integers as suffix array values) inside the genome folder as GENOME/GENOME.sa.bin
- Query files (one kmer at each line) with their tag number as GENOME/GENOME.QUERY_TAG.query.txt

Thus, if one have a fasta file at hand, output from `process_fasta`, `mksary` and `create_queries` (all are inside the `build` folder) can be used to format the files necessary for `snakemake`.
(To see how to use `process_fasta` and `create_queries`, follow the instructions from the `Reproducibility/README` on the `main` branch.)

For example, let the genome folder name be `ecoli`.
Then, we will have three files inside `ecoli` folder: `ecoli/ecoli.processed.fasta`, `ecoli/ecoli.sa.bin`, and `ecoli/ecoli.1.query.txt` (`1` being the tag number of the query)

Then, create a config file using CreateConfigFile.py. To see the config options:

```shell
python3 CreateConfigFile.py -h
```

Config options:

```text
usage: CreateConfigFile.py [-h] --genome_folder GENOME_FOLDER --epsilon EPSILON
                           [--kmer_size KMER_SIZE] [--code_path CODE_PATH] [--exec_path EXEC_PATH]
                           [--l_bits L_BITS] [--query_tag QUERY_TAG] [--sdsl_lib_path SDSL_LIB_PATH]
                           [--sdsl_inc_path SDSL_INC_PATH]

PLA-Index config file builder

options:
  -h, --help            show this help message and exit
  --genome_folder GENOME_FOLDER
                        Name of the folder containing the genome, suffix array and query files
                        [Required]
  --epsilon EPSILON     Epsilon value to be used to create pla-index [Required]
  --kmer_size KMER_SIZE
                        Kmer size for which pla-index will be calculated (default: 21)
  --code_path CODE_PATH
                        Path where the source code is (default: ../src/)
  --exec_path EXEC_PATH
                        Path where the executables will be stored (default: ../executables/)
  --l_bits L_BITS       How many elements to store in the shortcut array, D. |D| = 2^l (default: 16)
  --query_tag QUERY_TAG
                        QUERY_TAG to identify which query file to use from GENOME.QUERY_TAG.query.txt
                        (default: 1)
  --sdsl_lib_path SDSL_LIB_PATH
                        Path to the SDSL library folder (default: ~/lib)
  --sdsl_inc_path SDSL_INC_PATH
                        Path to the SDSL include folder (defaul: ~/include)
```

Here is an example using the `ecoli` folder inside the `tests` folder:

```shell
cd tests
mkdir -p ../executables
python3 ../CreateConfigFile.py --genome_folder ecoli --epsilon 15
snakemake -s ../Snakefile --cores=1 -p
```

## Citation

If you use `pla-index`, please cite:

```latex
@article {Abrar2024.02.08.579510,
    author = {Md Hasin Abrar and Paul Medvedev},
    title = {PLA-complexity of k-mer multisets},
    elocation-id = {2024.02.08.579510},
    year = {2024},
    doi = {10.1101/2024.02.08.579510}
}
```