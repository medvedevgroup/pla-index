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

void WriteQueryKmers_weighted(int32_t nQueries, string fn_suff,
        int kmer_size)
{
    int64_t idx;
    vector<string> query_kmers;
    srand(12345);
    for(int64_t i=0; i<nQueries; i++){
        idx = rand() % (fSize - kmer_size);
        query_kmers.push_back(GetKmerAtStringPos(idx, kmer_size));
    }
    
    cout<<"Writing start..\n";

    ofstream queryFile("1_QueryKmers_"+fn_suff+".txt");
    for(int i=0; i<query_kmers.size(); i++){
        
            // queryKmers.push_back(uniq_kmers[i]);
        queryFile<<"@read"<<(i+1)<<endl;
        queryFile<<query_kmers[i]<<endl;
        queryFile<<"+"<<endl;
        queryFile<<"999999999999999999999"<<endl;
    }
    queryFile.close();

    random_shuffle(query_kmers.begin(), query_kmers.end());
    ofstream queryFile2("2_QueryKmers_"+fn_suff+".txt");
    for(int i=0; i<query_kmers.size(); i++){
        
            // queryKmers.push_back(uniq_kmers[i]);
        queryFile2<<"@read"<<(i+1)<<endl;
        queryFile2<<query_kmers[i]<<endl;
        queryFile2<<"+"<<endl;
        queryFile2<<"999999999999999999999"<<endl;
    }
    queryFile2.close();

    random_shuffle(query_kmers.begin(), query_kmers.end());
    ofstream queryFile3("3_QueryKmers_"+fn_suff+".txt");
    for(int i=0; i<query_kmers.size(); i++){
        
            // queryKmers.push_back(uniq_kmers[i]);
        queryFile3<<"@read"<<(i+1)<<endl;
        queryFile3<<query_kmers[i]<<endl;
        queryFile3<<"+"<<endl;
        queryFile3<<"999999999999999999999"<<endl;
    }
    queryFile3.close();
}        

void WriteQueryKmers_unique(int32_t nQueries, string fn_suff,
        int kmer_size)
{
    vector<string> uniq_kmers;
    vector<string> query_kmers;
    FindUniqueKmers(uniq_kmers, kmer_size);
    // vector<uint64_t> kmerIndices(uniq_kmers.size(), 0);
    // iota(kmerIndices.begin(), kmerIndices.end(), 0);
    random_shuffle(uniq_kmers.begin(), uniq_kmers.end());
    // vector <uint64_t> queryIndices(kmerIndices.begin(),kmerIndices.begin() + nQueries);
    // sort(queryIndices.begin(), queryIndices.end());

    vector<string> queryKmers;
    int32_t indx = 0;

    cout<<"Writing start..\n";

    ofstream queryFile("1_QueryKmers_"+fn_suff+".txt");
    cout<<"#Unique kmers: "<<uniq_kmers.size()<<endl;
    for(int i=0; i<uniq_kmers.size(); i++){
        
            // queryKmers.push_back(uniq_kmers[i]);
        queryFile<<"@read"<<(indx+1)<<endl;
        queryFile<<uniq_kmers[i]<<endl;
        queryFile<<"+"<<endl;
        queryFile<<"999999999999999999999"<<endl;
        indx++;
        query_kmers.push_back(uniq_kmers[i]);
        if(indx == nQueries) break;
    }
    cout<<"Total number of queries: "<<indx<<endl;
    queryFile.close();

    random_shuffle(query_kmers.begin(), query_kmers.end());
    ofstream queryFile2("2_QueryKmers_"+fn_suff+".txt");
    for(int i=0; i<query_kmers.size(); i++){
        queryFile2<<"@read"<<(indx+1)<<endl;
        queryFile2<<query_kmers[i]<<endl;
        queryFile2<<"+"<<endl;
        queryFile2<<"999999999999999999999"<<endl;
    }
    queryFile2.close();

    random_shuffle(query_kmers.begin(), query_kmers.end());
    ofstream queryFile3("3_QueryKmers_"+fn_suff+".txt");
    for(int i=0; i<query_kmers.size(); i++){
        queryFile3<<"@read"<<(indx+1)<<endl;
        queryFile3<<query_kmers[i]<<endl;
        queryFile3<<"+"<<endl;
        queryFile3<<"999999999999999999999"<<endl;
    }
    queryFile3.close();
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
    fSize = filesize(genomeFile.c_str());
    printf("File size: %ld\n",fSize);
    LoadGenomeString(genomeFile);
    LoadSA(sa_file);
    cout<<"SA: "<<sa.size()<<endl;

    WriteQueryKmers_weighted(nQueries,fn_suffix, kmer_size);
    cout<<"Query done"<<endl;
    return 0;
}
