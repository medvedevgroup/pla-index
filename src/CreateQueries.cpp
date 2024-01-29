/*
* Suffix array random queries
*/

#include <bits/stdc++.h>
using namespace std;

vector<char> str;
vector<int32_t> sa;
size_t fSize;

string GetKmerAtIndex(int32_t i, int32_t kmer_size){
    int32_t suff_i_start = sa[i];
    int32_t suff_i_end;
    if(str.begin()+suff_i_start+kmer_size <= str.end()){
        string suffix_i(str.begin()+suff_i_start, str.begin()+suff_i_start+kmer_size);
        return suffix_i;
        // suff_i_end = suff_i_start+kmer_size;
    }
    return "-1"; // substring size < k: ignore these
}

void FindUniqueKmers(vector<string> &uniq_kmers,
                    int kmer_size){
    string kmer, prevKmer=""; 
    for(size_t i =0; i<sa.size(); i++){
        kmer = GetKmerAtIndex(i, kmer_size);
        if(kmer == "-1" || kmer == prevKmer) continue;
        uniq_kmers.push_back(kmer);
        prevKmer = kmer;
    }
}

string GetKmerAtStringPos(int64_t idx, int64_t kmer_size){
    string query(str.begin()+idx, str.begin()+idx+kmer_size);
    return query;
}

void WriteQueryKmers(int32_t nQueries, string fn_suff,
        int kmer_size, int nFiles)
{
    int64_t idx;
    vector<string> query_kmers;
    srand(12345);
    for(int64_t i=0; i<nQueries; i++){
        idx = rand() % (fSize - kmer_size);
        query_kmers.push_back(GetKmerAtStringPos(idx, kmer_size));
    }
    
    cout<<"Writing start..\n";

    for(int64_t n=1; n<=nFiles; n++){
        ofstream queryFile(fn_suff+"."+ to_string(n)+".query.txt");
        for(int i=0; i<query_kmers.size(); i++){
            // queryFile<<"@read"<<(i+1)<<endl;
            queryFile<<query_kmers[i]<<endl;
            // queryFile<<"+"<<endl;
            // queryFile<<"999999999999999999999"<<endl;
        }
        queryFile.close();

        random_shuffle(query_kmers.begin(), query_kmers.end());
    }
    
}

void LoadSA(string sa_file){
    sa.resize(fSize);
    FILE *fp = fopen(sa_file.c_str(), "rb");
    size_t err = fread(&sa[0], sizeof(int32_t), fSize, fp);
    fclose(fp);
}

void LoadGenomeString(string genomeFile){
    str.resize(fSize);
    FILE *fp = fopen(genomeFile.c_str(), "rb");
    //Last character is not handled -> should be $ or \0?
    size_t err = fread(&str[0], 1, fSize, fp);
    fclose(fp);
}

std::ifstream::pos_type filesize(const char* filename)
{
    std::ifstream in(filename, std::ifstream::ate | std::ifstream::binary);
    return in.tellg(); 
}


int main(int argc, char **argv){
    string genomeFile = argv[1]; 
    string sa_file = argv[2];
    string fn_suffix = argv[3];
    int kmer_size = stoi(argv[4]);
    int32_t nQueries = stoi(argv[5]);
    int64_t nFiles = stoi(argv[6]);
    fSize = filesize(genomeFile.c_str());
    LoadGenomeString(genomeFile);
    LoadSA(sa_file);    
    WriteQueryKmers(nQueries,fn_suffix, kmer_size, nFiles);
    cout<<nFiles<<" Query File(s) Created"<<endl;
    return 0;
}