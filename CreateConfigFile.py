import argparse

def CreateConfig(f, args, flags):
    f.write('---\n')
    f.write('CodePath: \"'+args.code_path+'\"\n')
    f.write('ExecPath: \"'+args.exec_path+'\"\n')
    f.write('Kmer_Size: \"'+str(args.kmer_size)+'\"\n')
    
    f.write('Lookup: \"'+str(args.lookup)+'\"\n')
    f.write('Query_Tag: \"'+str(args.query_tag)+'\"\n')
    f.write('Eps: \"'+str(args.epsilon)+'\"\n')
    f.write('SDSL_Inc: \"'+args.sdsl_inc_path+'\"\n')
    f.write('SDSL_Lib: \"'+args.sdsl_lib_path+'\"\n')
    f.write('Compile_Flag: \"'+flags+'\"\n')
    f.write('Genomes: \n')
    f.write('  - '+args.genome_folder+'\n')



def main():
    
    parser = argparse.ArgumentParser(description="PLA-Index config file builder")
    
    parser.add_argument('--genome_folder', type=str, required=True, help="Name of the folder containing the genome, suffix array and query files [Required]")
    parser.add_argument('--epsilon', type=int, required=True, help="Epsilon value to be used to create pla-index [Required]")

    parser.add_argument('--kmer_size', type=int, default=21, help="Kmer size for which pla-index will be calculated (default: 21)")
    parser.add_argument('--code_path', type=str, default="../src/", help="Path where the source code is (default: ../src/)")
    parser.add_argument('--exec_path', type=str, default="../executables/", help="Path where the executables will be stored (default: ../executables/)")  
    parser.add_argument('--lookup', type=int, default=16, help="On average on how many elements the binary search on X array will take place. Used to determine the prefix lookup table. (default: 16)")
    parser.add_argument('--query_tag', type=int, default=1, help="QUERY_TAG to identify which query file to use from GENOME.QUERY_TAG.query.txt (default: 1)")    
    
    parser.add_argument('--sdsl_lib_path', type=str, default='~/lib', help="Path to the SDSL library folder (default: ~/lib)")
    parser.add_argument('--sdsl_inc_path', type=str, default='~/include', help="Path to the SDSL include folder (defaul: ~/include)")

    args = parser.parse_args()
    
    flags = "-O3 -mbmi2 -msse4.2 -DNDEBUG -march=native -pthread"
    
    with open('config.yaml', 'w') as f:
        CreateConfig(f, args, flags)

main()