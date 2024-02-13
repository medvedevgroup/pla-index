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
    int64_t lp_bits = opt.lp_bits;

    // bool isFirstIndexReturned = false;
    INDX_TYPE it = EXACT_PLA;

    suffix_array<int64_t> sa;
    sa.Load(gn_fn, sa_fn);
    sa.set_kmer_size(kmer_size);

    pla_index pla(eps, lp_bits, sa.get_largest(), sa.size(), it);

    auto s1 = std::chrono::system_clock::now();
    pla.build_index(sa.begin());
    auto s2 = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = s2-s1;
    cout<<"Number of segments: "<<pla.get_num_knots()<<endl;
    cout<<"Total size in bytes: "<<pla.get_total_size_in_bytes()<<endl;
    // ofstream indx_time(indx_fn+"_build_time.txt", std::ios_base::app);
    // indx_time<<"Indx file: "<<indx_fn<<endl;
    // indx_time<<"Elapsed seconds in index construction: "<<elapsed_seconds.count()<<" sec"<<endl;
    cout<<"Elapsed seconds in index construction: "<<elapsed_seconds.count()<<" sec"<<endl;

    pla.Save(indx_fn);

    return 0;
}