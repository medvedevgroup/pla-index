# pla-index

pla-index builds a lightweight index from a FASTA file and its corresponding suffix array. The index supports k-mer membership queries and suffix array searches. It is substantially faster than a binary search on the suffix array.

## Requirements

- Linux (64 bit)
- C++17
- [SDSL](https://github.com/simongog/sdsl-lite/tree/master) (>=2.1.1)
- cmake (>= 3.16)

## Quick Start

### Installation

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

In case you want to construct your own suffix array, you may use [libdivsufsort](https://github.com/hasin-abrar/libdivsufsort) tool.
We modify it slightly to allow 64 bit integers.

To create the suffix array builder executable `mksary` from the current `build` folder:

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
    -t [STRING], --index_type=[STRING]
                Whether to build "basic-pla" or "repeat-pla" index. "repeat-pla" is
                generally the better option but can slow down rank queries unless the
                -r option is also used. [default: basic-pla]
    -k [INT], --kmer_size=[INT]
                Kmer size to be used to construct the index. [default: 21]
    -e [INT], --eps=[INT]
                Epsilon value to be used for constructing the pla-index. [default: 15]
    -o [STRING], --index=[STRING]
                File name where to save the index. [default: genome_fasta.index]
    -l [INT], --lookup=[INT]
                The average number of elements that will be searched on the X array. It is used to determine the size of the prefix lookup table. [default: 16]
    -r            Construct an extra bitvector with the same length as the suffix array.
                Use this to speed up rank queries. [default: disabled]
```

Similarly, query command line interface can be seen by:

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
    --query_type=[STRING]
                Whether to do "search" or "rank" query. "rank" gives you the value of
                the first position in the suffix array where the query is found.
                "search" on the other hand, also returns the value where the query
                appears, but it might not be the very first spot. [default: search]
```

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
From now on, we will be using the executables that are inside the `pla-index/build` folder.
To build suffix array:

```shell
cd build/
./mksary ../tests/ecoli/ecoli.processed.fasta ../tests/ecoli/ecoli.sa.bin
```

The built suffix array is written on `tests/ecoli/ecoli.sa.bin`.

To build a `basic-pla` index with default parameters (k=`21`, eps = `15`), we will use the formatted fasta file and the built suffix array:

```shell
./build_pla_index -g ../tests/ecoli/ecoli.processed.fasta -s ../tests/ecoli/ecoli.sa.bin -o ../tests/ecoli/ecoli.index
```

The built index will be saved on `../tests/ecoli/ecoli.index` file.

To do a `search` query with the provided `tests/ecoli/ecoli.1.query.txt` on the built `tests/ecoli/ecoli.index` one can do the following:

```shell
./query_pla_index -g ../tests/ecoli/ecoli.processed.fasta -s ../tests/ecoli/ecoli.sa.bin -i ../tests/ecoli/ecoli.index -q ../tests/ecoli/ecoli.1.query.txt
```

The query time will be shown on console.

## Reproducibility


To reproduce the results related to the evaluation of the PLA-index, please follow the `README` file inside the [Reproducibility](Reproducibility/README.md) folder.
To reproduce the results related to the PLA complexity of various datasets, please follow the `README` file inside the [Reproducibility-PLA-Complexity](Reproducibility-PLA-Complexity/README.md) folder.


## Citation

If you use `pla-index`, please cite:

```text
@article{Abrar2024.02.08.579510,
    author = {Md Hasin Abrar and Paul Medvedev},
    title = {PLA-complexity of k-mer multisets},
    elocation-id = {2024.02.08.579510},
    year = {2024},
    doi = {10.1101/2024.02.08.579510}
}
```
