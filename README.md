# pla-index
pla-index builds a lightweight index from a FASTA file and its corresponding suffix array. The index supports k-mer membership queries and suffix array searches. It is substantially faster than a binary search on the suffix array.

# Requirements
- Linux (64 bit)
- C++17
- [SDSL](https://github.com/simongog/sdsl-lite/tree/master)
- cmake (>= 3.16)

# Quick Start

## Installation
```shell
git clone --recursive https://github.com/medvedevgroup/pla-index.git
cd pla-index
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

To use `PLA-index-exact`, check out to the `pla-index-exact` branch, and follow the `README` instructions there.

## Usage

Now, all the executables are inside the `pla-index/build` folder. 
```shell
./build_index -h
```
would display the command line interface:
```
build_pla_index {OPTIONS}

build_pla_index

OPTIONS:

    -h, --help    Print help and exit
    -g [STRING], --genome_fasta=[STRING]
                Fasta file with one entry and only ACGT characters [Required]
    -s [STRING], --suffix_array=[STRING]
                Suffix array of the genome in a binary file [Required]
    -k [INT], --kmer_size=[INT]
                Kmer size to be used to construct the index [21]
    -e [INT], --eps=[INT]
                Epsilon value to be used for constructing the pla-index [15]
    -i [STRING], --index=[STRING]
                File name where to save the index [genome_fasta.index]
    -l [INT], --l_val=[INT]
                To determine the size of the shortcut array, D. |D| = 2^l. [16]
    -t [STRING], --index_type=[STRING]
                Whether to build "basic-pla" or "repeat-pla" index [basic-pla]
    -r            Whether to build bit vector to support fast rank query. [disabled]
```

<!-- Parameter description:

| Parameter Name | Description |
|----------|----------|
| GENOME-FASTA-FILE | Fasta file with one entry and only ACGT characters|
| SUFFIX-ARRAY-FILE  |  Suffix array for GENOME in a binary file |
| KMER-SIZE | Kmer size to be used to construct the index |
| EPS | Epsilon value to be used for constructing the pla-index |
| INDEX-NAME | File name where to save the index |
| L-VALUE | To determine the size of the shortcut array, D. &#124;D&#124; = 2<sup>l</sup> |
| INDEX-TYPE | Whether to build "basic-pla" or "repeat-pla" index |
| ENABLE-FAST-RANK | Whether to build bit vector to support fast rank query. Provide either "y" or "n" for yes and  no respectively | -->

For example, to build a `basic-pla` index with `21` size k-mer, default epsilon value (15), default shortcut array size (|D| = 2<sup>16</sup>) and not creating a bit vector for faster rank query afterwards on the `ecoli` genome inside the `tests/ecoli/` folder:
```shell
./build_index -g ../tests/ecoli/ecoli.processed.fasta -s ../tests/ecoli/ecoli.sa.bin -i ../tests/ecoli/ecoli.index
```
The built index will be saved on `../tests/ecoli/ecoli.index` file.

Similarly, query command line interface can be seen by:
```shell
./query_index -h
```
and the interface:
```
query_pla_index {OPTIONS}

query_pla_index

OPTIONS:

    -h, --help    Print help and exit
    -g [STRING], --genome_fasta=[STRING]
                Fasta file with one entry and only ACGT characters [Required]
    -s [STRING], --suffix_array=[STRING]
                Suffix array of the genome in a binary file [Required]
    -i [STRING], --index=[STRING]
                File name where to save the index [genome_fasta.index] [Required]
    -q [STRING], --query=[STRING]
                Name of the file containing all the queries (one kmer per line)
                [Required]
    --query_type=[STRING]
                Whether to do "search" or "rank" query [search]
```
<!-- 
To query the index:
```
./query_index GENOME-FASTA-FILE SUFFIX-ARRAY-FILE QUERY-FILE INDEX-NAME QUERY-TYPE
```
Parameter description:

| Parameter Name | Description |
|----------|----------|
| GENOME-FASTA-FILE | Fasta file with one entry and only ACGT characters|
| SUFFIX-ARRAY-FILE  |  Suffix array for GENOME in a binary file |
| QUERY-FILE | Name of the file containing all the queries (one kmer per line) |
| INDEX-NAME | Index file to use |
| QUERY-TYPE | Whether to do "search" or "rank" query | -->

For example, to do `search` query with the provided `tests/ecoli/ecoli.1.query.txt` on the built `tests/ecoli/ecoli.index` one can do the following:
```shell
./query_index -g ../tests/ecoli/ecoli.processed.fasta -s ../tests/ecoli/ecoli.sa.bin -i ../tests/ecoli/ecoli.index -q ../tests/ecoli/ecoli.1.query.txt
```
The query time information will be showed on the console.

We follow the following formats at the input for both building and querying pla-index:
- The genome file is a fasta file with only one header line. There is no 'N' character at the nucleotide bases in the fasta file.
- The suffix array file is a binary file of consecutive 64 bit integers (values of the suffix array indices)
- Index type can be either "basic-pla" or "repeat-pla"
- Query file contains same length kmers (one per line). The kmer-size of the query kmers is the same as the provided one while building the index

To construct suffix array, we use [libdivsufsort](https://github.com/hasin-abrar/libdivsufsort) tool. 
We modify it slightly to allow 64 bit integers. 

To construct the suffix array builder executable `mksary`:
```
cd ../libdivsufsort
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE="Release" \
-DCMAKE_INSTALL_PREFIX="/usr/local" ..
sed -i 's/int32_t/int64_t/g' include/divsufsort.h
make
cp examples/mksary ../../build/
```

# Reproducibility

To reproduce the results shown in our paper, one can follow the `README` file inside the `Reproducibility` folder.


