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

## Usage

Now, all the executables are inside the `pla-index/build` folder. 

We follow the following format at the input for both building and querying pla-index:
- The genome is string is filtered by 'N', and concatenated to be a single string
- The suffix array is a binary file of consecutive 64 bit integers (values at the suffix array indices)
- Index type can be either "basic-pla" or "repeat-pla"
- Query type can be either "search" or "rank" query
- Query file contains same length kmers (one per line)

To build the index:
```
./build_index GENOME SUFFIX-ARRAY KMER-SIZE EPS INDEX-NAME L-VALUE INDEX-TYPE QUERY-TYPE

Parameter description:

GENOME: [string], Whole genome string concatenated filtered by 'N'
SUFFIX-ARRAY: [string], Suffix array for GENOME in a binary file
KMER-SIZE: [int], Kmer size to be used to construct the index
EPS: [int], Epsilon value to be used for constructing the pla-index
INDEX-NAME: [string], File name where to save the index
L-VALUE: [int], To determine the size of the shortcut array, D. |D| = 2^l
INDEX-TYPE: [string], Whether to build "basic-pla" or "repeat-pla"
QUERY-TYPE: [string], Whether to do "search" or "rank" query afterwards
```



To query the index:
```
./query_index GENOME SUFFIX-ARRAY KMER-SIZE QUERY-FILE INDEX-NAME RUNINFO-FILE INDEX-TYPE QUERY-TYPE

Parameter description:
GENOME: [string], Whole genome string concatenated filtered by 'N'
SUFFIX-ARRAY: [string], Suffix array for GENOME in a binary file
KMER-SIZE: [int], Kmer size of the query k-mers
QUERY-FILE: [string], Name of the file containing all the queries (one kmer per line)
INDEX-NAME: [string], Index file to use 
RUNINFO-FILE: [string], Where to store the runtime information
INDEX-TYPE: [string], Type of the stored index ("basic-pla" or "repeat-pla")
QUERY-TYPE: [string], Whether to do "search" or "rank" query
```

