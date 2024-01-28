fn = 'auto_strobe.sh'

strobe_indx_ln = 46

eps_list = [255, 63, 15]

with open(fn, 'w') as f:
    f.write("#!/bin/bash\n\n")

    for eps in eps_list:
        f.write("echo \"Executing for eps: "+str(eps)+"\"\n")
        f.write("sed -i \'"+str(strobe_indx_ln)+"s/.*/, rs_pla_index("+str(eps)+", 16)/\' src/index.hpp\n")
        f.write("make -j -C build\n")

        #f.write("/usr/bin/time -f \"%M,%E\" --output-file=strobe_sing_mem_wallT.txt  ./build/strobealign -t 16 --create-index tests/chm13/chm13v2.0.fa tests/chm13/reads.chm13.1.fq tests/chm13/reads.chm13.2.fq\n")
        for i in range(1, 6):
            f.write("/usr/bin/time -f \"%M,%E\" --output-file=strobe_sing_mem_wallT.txt  ./build/strobealign -t 16 tests/drosophila/ref.fasta tests/drosophila/reads.1.fastq.gz tests/drosophila/reads.2.fastq.gz | samtools view -o droso.bam\n")
            f.write("cat strobe_sing_mem_wallT.txt >> strobe_all_droso.txt\n")
        f.write('\n')

        for i in range(1, 6):
            f.write("/usr/bin/time -f \"%M,%E\" --output-file=strobe_sing_mem_wallT.txt  ./build/strobealign -t 16 tests/chm13/chm13v2.0.fa tests/chm13/reads.chm13.1.fq tests/chm13/reads.chm13.2.fq | samtools view -o chm13.bam\n")
            f.write("cat strobe_sing_mem_wallT.txt >> strobe_all_chm.txt\n")
        f.write('\n')



#        f.write("/usr/bin/time -f \"%M,%E\" --output-file=strobe_sing_mem_wallT.txt  ./build/strobealign -t 16 tests/drosophila/ref.fasta tests/drosophila/reads.1.fastq.gz tests/drosophila/reads.2.fastq.gz | samtools view -o test.bam\n")
#        f.write("cat strobe_sing_mem_wallT.txt >> strobe_all.txt\n")
    # f.write('\n')
