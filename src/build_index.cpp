#include <stdio.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include "pla_index.hpp"
#include "suffix_array.h"

int main(int argc, char **argv){
    string gn_fn = argv[1]; 
    string sa_fn = argv[2];
    int64_t kmer_size = stoll(argv[3]);
    int64_t eps = stoll(argv[4]);
    string indx_fn = argv[5];
    int64_t lp_bits = stoll(argv[6]);
    bool isFirstIndexReturned = false;
    INDX_TYPE it = EXACT_PLA;

    suffix_array<int64_t> sa;
    sa.Load(gn_fn, sa_fn, kmer_size);

    pla_index pla(eps, lp_bits, sa.get_largest(), isFirstIndexReturned, sa.size(), it);

    auto s1 = std::chrono::system_clock::now();
    pla.build_index(sa.begin());
    auto s2 = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = s2-s1;
    cout<<"Number of segments: "<<pla.get_num_knots()<<endl;
    cout<<"Total size in bytes: "<<pla.get_total_size_in_bytes()<<endl;
    ofstream indx_time(indx_fn+"_build_time.txt", std::ios_base::app);
    indx_time<<"Indx file: "<<indx_fn<<endl;
    indx_time<<"Elapsed seconds in index construction: "<<elapsed_seconds.count()<<" sec"<<endl;
    cout<<"Elapsed seconds in index construction: "<<elapsed_seconds.count()<<" sec"<<endl;

    pla.Save(indx_fn);

    return 0;
}