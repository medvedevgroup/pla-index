#include <stdio.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include "pla_index.hpp"
#include "suffix_array.h"

int main(int argc, char **argv){
    // NFiltered Concatenated single string
    // No new line at the last character
    string gn_fn = argv[1]; 
    string sa_fn = argv[2];
    int64_t kmer_size = stoi(argv[3]);
    int64_t eps = stoi(argv[4]);
    string indx_fn = argv[5];
    int64_t lp_bits = stoi(argv[6]);
    string indx_type = argv[7];
    string query_type = argv[8];
    bool isFirstIndexReturned;

    if (query_type == "search") isFirstIndexReturned = false;
    else if(query_type == "rank") isFirstIndexReturned = true;
    else{
        throw std::logic_error("query type can be either search or rank");
    }
    
    INDX_TYPE it;
    if (indx_type == "basic-pla") it = BASIC_PLA;    
    else if(indx_type == "repeat-pla") it = REPEAT_PLA;
    else{
        throw std::logic_error("indx type can be either basic-pla or repeat-pla");
    }

    suffix_array<int64_t> sa;
    sa.Load(gn_fn, sa_fn, kmer_size);

    pla_index pla(eps, lp_bits, sa.get_largest(), isFirstIndexReturned, sa.size(), it);

    auto s1 = std::chrono::system_clock::now();
    pla.build_index(sa.begin());
    auto s2 = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = s2-s1;

    ofstream indx_time(indx_fn+"_build_time.txt", std::ios_base::app);
    cout<<"Number of segments: "<<pla.get_num_segments()<<endl;
    indx_time<<"Indx file: "<<indx_fn<<endl;
    indx_time<<"Elapsed seconds in index construction: "<<elapsed_seconds.count()<<" sec"<<endl;
    cout<<"Elapsed seconds in index construction: "<<elapsed_seconds.count()<<" sec"<<endl;

    pla.Save(indx_fn);

    return 0;
}