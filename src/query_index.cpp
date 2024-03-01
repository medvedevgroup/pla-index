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

    string indx_fn = "pla_index.indx";
    pla_index pla;
    
    pla.Load(indx_fn, data);
    int64_t q = 42;
    std::cout<<pla.query(q, data)<<endl;

    return 0;
}