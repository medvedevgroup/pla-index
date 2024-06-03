#include <stdio.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include "pla_index.hpp"
#include "suffix_array.h"
#include "cmdline.hpp"

int main(int argc, char **argv){
    auto opt = parse_command_line_arguments_build(argc, argv);
    
    string gn_fn = opt.gn_fn; 
    string sa_fn = opt.sa_fn;
    int64_t kmer_size = opt.kmer_size;
    int64_t eps = opt.eps;
    string indx_fn = opt.indx_fn;
    int64_t lookup_count = opt.lookup_count;
    string indx_type = opt.indx_type;

    

    suffix_array<int64_t> sa;
    sa.Load(gn_fn, sa_fn);
    sa.set_kmer_size(kmer_size);

    pla_index pla(eps, lookup_count, sa.get_largest(), sa.size());

    auto s1 = std::chrono::system_clock::now();
    pla.build_index(sa.begin());
    auto s2 = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = s2-s1;

    cout<<"Number of segments: "<<pla.get_num_segments()<<endl;
    cout<<"Index stored in: "<<indx_fn<<endl;
    cout<<"Elapsed seconds in index construction: "<<elapsed_seconds.count()<<" sec"<<endl;

    pla.Save(indx_fn, sa);

    return 0;
}