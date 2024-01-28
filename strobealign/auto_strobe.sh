#!/bin/bash

echo "Executing for eps: 1023"
sed -i '46s/.*/, rs_pla_index(1023, 16)/' src/index.hpp
make -j -C build
/usr/bin/time -f "%M,%E" --output-file=strobe_sing_mem_wallT.txt  ./build/strobealign  tests/drosophila/ref.fasta tests/drosophila/reads.1.fastq.gz tests/drosophila/reads.2.fastq.gz | samtools view -o droso.bam
cat strobe_sing_mem_wallT.txt >> strobe_all_droso.txt
/usr/bin/time -f "%M,%E" --output-file=strobe_sing_mem_wallT.txt  ./build/strobealign  tests/drosophila/ref.fasta tests/drosophila/reads.1.fastq.gz tests/drosophila/reads.2.fastq.gz | samtools view -o droso.bam
cat strobe_sing_mem_wallT.txt >> strobe_all_droso.txt

