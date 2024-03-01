#include <stdio.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include "pla_index.hpp"

int main(int argc, char **argv){
    std::vector<int64_t> data(1000000);
    std::srand(12345);
    std::generate(data.begin(), data.end(), std::rand);
    data.push_back(42);
    std::sort(data.begin(), data.end());

    int64_t epsilon = 128;
    pla_index pla(epsilon, data);

    auto s1 = std::chrono::system_clock::now();
    pla.build_index(data.begin());
    auto s2 = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = s2-s1;

    cout<<"Number of segments: "<<pla.get_num_segments()<<endl;
    cout<<"Elapsed seconds in index construction: "<<elapsed_seconds.count()<<" sec"<<endl;

    string indx_fn = "pla_index.indx";
    pla.Save(indx_fn, data);

    int64_t q = 42;
    std::cout<<pla.query(q, data)<<endl;


    return 0;
}