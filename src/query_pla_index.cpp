#include <stdio.h>
#include <fstream>
#include "pla_index.hpp"
#include "suffix_array.h"
#include "cmdline.hpp"

using namespace std::chrono;

void LoadQuery(const std::string &query_file, std::vector<std::string>& query_kmers_vec, 
        std::vector<size_t>& query_kmerVals, suffix_array<int64_t> &sa){
    std::ifstream queryFile(query_file.c_str());
    std::string line, kmer;
    int indx = 0;
    
    while(getline(queryFile, line)){
        std::stringstream ss(line);
        ss >> kmer;
        query_kmers_vec.emplace_back(kmer);
        query_kmerVals.emplace_back(sa.GetKmerVal(kmer));
    }
    queryFile.close();
}

int main(int argc, char **argv){
    auto opt = parse_command_line_arguments_query(argc, argv);
    string gn_fn = opt.gn_fn; 
    string sa_fn = opt.sa_fn;
    string qfn = opt.query_file;
    string indx_fn = opt.indx_fn;
    string runInfo_fn;
    bool is_stored_on_file = false;
    if(opt.run_info_fn != "-1") is_stored_on_file = true;

    suffix_array<int64_t> sa;
    sa.Load(gn_fn, sa_fn);

    pla_index pla;
    pla.Load(indx_fn, sa);
    
    std::vector<std::string> query_kmers_vec;
    std::vector<size_t> query_kmerVals;

    LoadQuery(qfn, query_kmers_vec, query_kmerVals, sa);
    int64_t nQueries = query_kmers_vec.size();
    std::cout<<"Number of queries: "<<nQueries<<std::endl;
    int _correct = 0;
    
    std::vector<size_t> str_pos_vec(query_kmers_vec.size());
    auto s1 = std::chrono::system_clock::now();
    for(size_t i=0; i<nQueries; i++){
        str_pos_vec[i] = pla.query(query_kmerVals[i], sa);
        
    }
    auto s2 = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = s2-s1;
    
    std::cout<<"Elapsed_seconds in query: "<<elapsed_seconds.count()<<std::endl;
    
    for(size_t i=0; i<query_kmers_vec.size(); i++){
        if(sa.GetKmerAtStrPos(str_pos_vec[i]) == query_kmers_vec[i]){
            _correct++;
        }
        else{
            std::cout<<query_kmers_vec[i]<<std::endl;
            break;
        }
    }
    std::cout<<"Correct: "<<_correct<<" Total: "<<query_kmers_vec.size()<<std::endl;
    double accuracy = ((double)_correct/query_kmers_vec.size()) * 100;
    printf("Accuracy: %0.2f %%\n",accuracy);
    
    if(is_stored_on_file){
        std::ofstream runFile(opt.run_info_fn.c_str());
        runFile<<qfn<<endl;
        runFile<<"Elapsed_seconds in query: "<<elapsed_seconds.count()<<std::endl;
        runFile<<"Correct: "<<_correct<<" Total: "<<query_kmers_vec.size()<<" Accurancy: "<<accuracy<<"%"<<std::endl;
    }
    
    return 0;
}


