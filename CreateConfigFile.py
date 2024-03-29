import argparse

# flag: -Ofast -mbmi2 -msse4.2 -DNDEBUG -march=native -I ~/include -L ~/lib 
def CreateConfig(f, args, flags):
    f.write('---\n')
    f.write('CodePath: \"'+args.code_path+'\"\n')
    f.write('ExecPath: \"'+args.exec_path+'\"\n')
    f.write('Kmer_Size: \"'+str(args.kmer_size)+'\"\n')
    f.write('LP_Bits: \"'+str(args.l_bits)+'\"\n')
    f.write('Query_Tag: \"'+str(args.query_tag)+'\"\n')
    f.write('Eps: \"'+str(args.epsilon)+'\"\n')
    f.write('Indx_Type: \"'+args.index_type+'\"\n')
    f.write('Fast_Rank: \"'+args.fast_rank+'\"\n')
    f.write('Query_Type: \"'+args.query_type+'\"\n')
    f.write('SDSL_Inc: \"'+args.sdsl_inc_path+'\"\n')
    f.write('SDSL_Lib: \"'+args.sdsl_lib_path+'\"\n')
    f.write('Compile_Flag: \"'+flags+'\"\n')
    f.write('Genomes: \n')
    f.write('  - '+args.genome_folder+'\n')



def main():
    
    parser = argparse.ArgumentParser(description="PLA-Index config file builder")
    
    parser.add_argument('--genome_folder', type=str, required=True, help="Name of the folder containing the genome, suffix array and query files [Required]")
    parser.add_argument('--epsilon', type=int, required=True, help="Epsilon value to be used to create pla-index [Required]")

    parser.add_argument('--index_type', type=str, default='basic-pla', help="Whether to build \"basic-pla\" or \"repeat-pla\" index. \"repeat-pla\" is generally the better option but can slow down rank queries unless the -r option is also used.(default: basic-pla)")
    parser.add_argument('--fast_rank', type=str, default='n', help="Construct an extra bitvector with the same length as the suffix array. Use this to speed up rank queries. Possible values: \"y\" for yes and \"n\" for no. (default: n)")
    parser.add_argument('--query_type', type=str, default='search', help="Whether to do \"search\" or \"rank\" query. \"rank\" gives you the value of the first position in the suffix array where the query is found. \"search\" on the other hand, also returns the value where the query appears, but it might not be the very first spot. (default: search)")

    parser.add_argument('--kmer_size', type=int, default=21, help="Kmer size for which pla-index will be calculated (default: 21)")
    parser.add_argument('--code_path', type=str, default="../src/", help="Path where the source code is (default: ../src/)")
    parser.add_argument('--exec_path', type=str, default="../executables/", help="Path where the executables will be stored (default: ../executables/)")  
    parser.add_argument('--l_bits', type=int, default=16, help=" How many elements to store in the shortcut array, D. |D| = 2^l (default: 16)")
    parser.add_argument('--query_tag', type=int, default=1, help="QUERY_TAG to identify which query file to use from GENOME.QUERY_TAG.query.txt (default: 1)")        
    parser.add_argument('--sdsl_lib_path', type=str, default='~/lib', help="Path to the SDSL library folder (default: ~/lib)")
    parser.add_argument('--sdsl_inc_path', type=str, default='~/include', help="Path to the SDSL include folder (defaul: ~/include)")
    
    

    args = parser.parse_args()
    
    flags = "-Ofast -mbmi2 -msse4.2 -DNDEBUG -march=native"
    
    with open('config.yaml', 'w') as f:
        CreateConfig(f, args, flags)

main()