# pla-index
pla-index allows faster query of a k-mer rank function using an index that is created using piece-wise linear approximation.

# Requirements
- Linux (64 bit)
- C++17
- [SDSL](https://github.com/simongog/sdsl-lite/tree/master)
- cmake (>= 3.16)
- [Snakemake](https://snakemake.readthedocs.io/en/stable/) (for automating all computations)

# Quick Start

## Installation

Clone the repo using:

```shell
git clone --recursive https://github.com/medvedevgroup/pla-index.git
```

This repository has two other branches for read aligner and exact-pla.\
To use the read aligner application, run the following command:
```shell
git checkout strobealign-application
```

To use the exact-pla, run the following command:
```shell
git checkout pla-index-exact
```

To compile from sources, at first update the CMakeLists.txt with SDSL library and include path.
This can be done by running the following script:

```
cd pla-index
./update_cmake_lists.sh SDSL_INCLUDE_PATH SDSL_LIB_PATH

Parameter descriptions:
SDSL_INCLUDE_PATH: Path to the SDSL include folder [typically ~/include/]
SDSL_LIB_PATH: Path to the SDSL library folder [typically ~/lib/]
```

After updating CMakeLists.txt, all files can be compiled using cmake.

```shell
mkdir build; cd build; cmake ..; make -j 4
```

To construct suffix array we use [libdivsufsort](https://github.com/hasin-abrar/libdivsufsort) tool. 
We modify it slightly to allow 64 bit integers. 
To construct the suffix array builder executable `mksary`:
```
cd ../libdivsufsort
mkdir build_sa
cd build_sa
cmake -DCMAKE_BUILD_TYPE="Release" \
-DCMAKE_INSTALL_PREFIX="/usr/local" ..
sed -i 's/int32_t/int64_t/g' include/divsufsort.h
make
cp examples/mksary ../../build/
```

## Usage

Now, all the executables are inside the `build` folder. 
Here on, we assume we are inside this folder.

To filter 'N' from a genome:
```
./filterNConcat FASTA-FILE OUTPUT-PREFIX

Parameter description:
FASTA-FILE: Fasta file on which pla-index will be built
OUTPUT-PREFIX: The output file name format is OUTPUT-PREFIX.ConcatenatedGenome.txt
```

To construct suffix array on the concatenated genome:

```
./mksary GENOME SUFFIX-ARRAY

Parameter description:
GENOME: Concatenated genome without 'N'
SUFFIX-ARRAY: Name of the binary file where suffix array will be stored
```

To build the index:
```
./build_index GENOME SUFFIX-ARRAY KMER-SIZE EPS INDEX-NAME L-VALUE INDEX-TYPE QUERY-TYPE

Parameter description:

GENOME: Whole genome string concatenated filtered by 'N'
SUFFIX-ARRAY: Suffix array for GENOME in a binary file
KMER-SIZE: Kmer size to be used to construct the index
EPS: Epsilon value to be used for constructing the pla-index
INDEX-NAME: File name where to save the index
L-VALUE: To determine the size of the shortcut array, D. |D| = 2^l
INDEX-TYPE: Whether to build "basic-pla" or "repeat-pla"
QUERY-TYPE: Whether to do "search" or "rank" query afterwards
```

To construct query files:
```
./create_queries GENOME SUFFIX-ARRAY PREFIX KMER-SIZE NO-OF-QUERIES NO-OF-FILES

Parameter description:

GENOME: Whole genome string concatenated filtered by 'N'
SUFFIX-ARRAY: Suffix array for GENOME in a binary file
PREFIX: File name prefix of the query file
KMER-SIZE: Kmer size to be used to construct the index
NO-OF-QUERIES: How many query k-mers to randomly construct
NO-OF-FILES: How many query files to construct
```

To query the index:
```
./benchmark_index GENOME SUFFIX-ARRAY KMER-SIZE QUERY-FILE INDEX-NAME RUNINFO-FILE INDEX-TYPE QUERY-TYPE

Parameter description:
GENOME: Whole genome string concatenated filtered by 'N'
SUFFIX-ARRAY: Suffix array for GENOME in a binary file
KMER-SIZE: Kmer size of the query k-mers
INDEX-NAME: Index file to use 
RUNINFO-FILE: Where to store the runtime information
INDEX-TYPE: Type of the stored index ("basic-pla" or "repeat-pla")
QUERY-TYPE: Whether to do "search" or "rank" query
```

## Snakemake

To automate the whole process using snakemake, create a folder containing the genome, binary file of the suffix array and query files.
The snakemake file assumes the following naming conventions:
- Genome folder name without any white-spaces(example: GENOME)
- Concatenated genome without N inside the genome folder needs to be named as GENOME/GENOME.ConcatenatedGenome.txt
- Binary suffix array file inside the genome folder as GENOME/GENOME.sa.bin
- Query files with their tag number as GENOME/GENOME.QUERY_TAG.query.txt

For example, let the genome folder name be `ecoli`. 
Then, we will have three files inside `ecoli` folder: `ecoli/ecoli.ConcatenatedGenome.txt`, `ecoli/ecoli.sa.bin`, and `ecoli/ecoli.1.query.txt` (`1` being the tag number of the query)

Then, create a config file using CreateConfigFile.py

Required Parameters:

| Parameter  | Type    | Description    |
|-------------|-------------|-------------|
|--genome_folder | [String] |The name of the folder containing the genome, suffix array and query files|
|--epsilon |  [int]   |Epsilon value to be used to create pla-index|

Optional parameters with required argument:

| Parameter  | Type    | Description    |
|-----------------|-------------|-------------|
|--index_type |[String] | What kind of index type to construct and/or use (default: basic-pla). Possible values: "basic-pla" or "repeat-pla"|
|--query_type |[String] | What kind of query to do (default: search). Possible values: "search" or "rank"|
|--kmer_size |[int] | Kmer size for which pla-index will be calculated (default: 21)|
|--code_path |[String] | Path where the source code is (default: ../src/)|
|--exec_path |[String] | Path where the executables will be stored (default: ../executables/)|
|--l_bits |[int] | How many elements to store in the shortcut array, D. &#124;D&#124; = 2<sup>l</sup> (default: 16)|
|--query_tag |[int] | Which query file to use (default: 1)|
|--sdsl_lib_path |[String] | Path to the SDSL library folder (default: ~/lib)|
|--sdsl_inc_path |[String] | Path to the SDSL include folder (defaul: ~/include)|

Here is an example assuming we have the `ecoli` folder inside the `tests` folder:

```
cd tests
mkdir -p ../executables
python3 ../CreateConfigFile.py --genome_folder ecoli --epsilon 15
snakemake -s ../Snakefile --cores=1 -p
```
