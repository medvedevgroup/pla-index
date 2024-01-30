#include <stdio.h>
#include <fstream>
#include "pla_index.hpp"
#include "suffix_array.h"

using namespace std::chrono;

void LoadQuery(const std::string &query_file, std::vector<std::string>& query_kmers_vec, 
        std::vector<size_t>& query_kmerVals, suffix_array<int64_t> &sa){
    std::ifstream queryFile(query_file.c_str());
    std::string line, kmer;
    int indx = 0;
    
    while(getline(queryFile, line)){
        // if(indx == 1){
            std::stringstream ss(line);
            ss >> kmer;
            query_kmers_vec.emplace_back(kmer);
            query_kmerVals.emplace_back(sa.GetKmerVal(kmer));
        // }
        // indx++;
        // indx %= 4;
    }
    queryFile.close();
}

int main(int argc, char **argv){
    string gn_fn = argv[1];
    string sa_fn = argv[2];
    int64_t kmer_size = stoll(argv[3]);
    string query_fn = argv[4];
    string indx_fn = argv[5];
    string runInfo_fn = argv[6];
    string indx_type = argv[7];
    string query_type = argv[8];
    // int64_t eps = stoi(argv[10]);
    int64_t knot_bs_thres = 64;
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
        throw std::logic_error("indx type can be either \"basic-pla\" or \"repeat-pla\"");
    }

    suffix_array<int64_t> sa;
    sa.Load(gn_fn, sa_fn, kmer_size);

    pla_index pla(knot_bs_thres, sa.get_largest(), sa.size(), indx_fn, it);

    std::vector<std::string> query_kmers_vec; 
    std::vector<size_t> query_kmerVals;

    LoadQuery(query_fn, query_kmers_vec, query_kmerVals, sa);

    int64_t nQueries = query_kmers_vec.size();
    std::cout<<"Number of queries: "<<nQueries<<std::endl;
    std::cout<<"Starting benchmarking...\n";
    // long long sa_idx;
    int _correct = 0;
    
    std::vector<size_t> str_pos_vec(query_kmers_vec.size());
    auto s1 = std::chrono::system_clock::now();
    for(size_t i=0; i<nQueries; i++){
        str_pos_vec[i] = pla.query(query_kmerVals[i], query_kmers_vec[i], sa);
    }
    auto s2 = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = s2-s1;
    std::ofstream runFile(runInfo_fn.c_str(), std::ios_base::app);
    // runFile<<"Eps: "<<eps<<std::endl;
    runFile<<"Elapsed_seconds in query: "<<elapsed_seconds.count()<<std::endl;    
    std::cout<<"Elapsed_seconds in query: "<<elapsed_seconds.count()<<std::endl;    
    for(size_t i=0; i<query_kmers_vec.size(); i++){
        if(sa.GetKmerAtStrPos(str_pos_vec[i]) == query_kmers_vec[i]){
            _correct++;
        }
        else{
            std::cout<<query_kmers_vec[i]<<std::endl;
        }
    }
    std::cout<<"Correct: "<<_correct<<" Total: "<<query_kmers_vec.size()<<std::endl;
    double accuracy = ((double)_correct/query_kmers_vec.size()) * 100;
    printf("Accuracy: %0.2f %%\n",accuracy);
    runFile<<"Correct: "<<_correct<<" Total: "<<query_kmers_vec.size()<<" Accurancy: "<<accuracy<<"%"<<std::endl;

    // COMMENT LATER
    // pla.save_unpacked(indx_fn);

    return 0;
}


