import sys
import argparse

# flag: -Ofast -mbmi2 -msse4.2 -DNDEBUG -march=native -I ~/include -L ~/lib 
def CreateConfig(f, args, flags):
    f.write('---\n')
    f.write('CodePath: \"'+args.code_path+'\"\n')
    f.write('ExecPath: \"'+args.exec_path+'\"\n')
    f.write('Kmer_Size: \"'+str(args.kmer_size)+'\"\n')
    f.write('Knot_BS_Thres: \"'+str(args.knot_bs_thres)+'\"\n')
    f.write('LP_Bits: \"'+str(args.l_bits)+'\"\n')
    f.write('Query_Tag: \"'+str(args.query_tag)+'\"\n')
    f.write('Eps: \"'+str(args.epsilon)+'\"\n')
    f.write('Indx_Type: \"'+args.index_type+'\"\n')
    f.write('Query_Type: \"'+args.query_type+'\"\n')
    f.write('SDSL_Inc: \"'+args.sdsl_inc_path+'\"\n')
    f.write('SDSL_Lib: \"'+args.sdsl_lib_path+'\"\n')
    f.write('Compile_Flag: \"'+flags+'\"\n')
    f.write('Genomes: \n')
    f.write('  - '+args.genome_folder+'\n')



def main():
    
    parser = argparse.ArgumentParser(description="PLA-Index config file builder")
    
    parser.add_argument('--genome_folder', type=str, default=None)
    parser.add_argument('--epsilon', type=int, default=None)

    parser.add_argument('--kmer_size', type=int, default=21)
    parser.add_argument('--code_path', type=str, default="../src/")
    parser.add_argument('--exec_path', type=str, default="../executables/")
    parser.add_argument('--knot_bs_thres', type=int, default=64)    
    parser.add_argument('--l_bits', type=int, default=16)
    parser.add_argument('--query_tag', type=int, default=1)    
    parser.add_argument('--index_type', type=str, default='basic-pla')
    parser.add_argument('--query_type', type=str, default='search')
    parser.add_argument('--sdsl_lib_path', type=str, default='~/lib')
    parser.add_argument('--sdsl_inc_path', type=str, default='~/include')

    args = parser.parse_args()
    
    flags = "-Ofast -mbmi2 -msse4.2 -DNDEBUG -march=native"
    # kmer_size = int(sys.argv[1])
    # flags = sys.argv[2]
    # codePath = sys.argv[3]
    # execPath = sys.argv[4]
    # repeat_process = sys.argv[5]
    # topk_thres = sys.argv[6]
    # sa_thres = sys.argv[7]
    # sample_name = sys.argv[8]
    # nBuckets = int(sys.argv[9])
    # lp_bits = int(sys.argv[10])
    # query_tag = sys.argv[11]
    # err_thres = int(sys.argv[12])
    # allowed_rep = int(sys.argv[13])
    
    with open('config.yaml', 'w') as f:
        CreateConfig(f, args, flags)

main()