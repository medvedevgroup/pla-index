#pragma once
#include <vector>
#include <string>

template <typename T>
class suffix_array {
private:
    std::vector<T> sa;
    std::vector<char> str;
    size_t fSize;
    int64_t kmer_size;

public:
    class Iterator {
    private:
        typename std::vector<T>::iterator it;
        suffix_array& sa;

    public:
        Iterator(typename std::vector<T>::iterator iter, suffix_array& sa) 
                : it(iter), sa(sa) {}

        // Overload prefix increment operator (++it)
        Iterator& operator++() {
            ++it;
            return *this;
        }

        // Overload postfix increment operator (it++)
        Iterator operator++(int) {
            Iterator temp = *this;
            ++it;
            return temp;
        }

        // Overload prefix decrement operator (++it)
        Iterator& operator--() {
            --it;
            return *this;
        }

        // Overload postfix increment operator (it++)
        Iterator operator--(int) {
            Iterator temp = *this;
            --it;
            return temp;
        }

        // Overload -= operator (it -= n)
        Iterator& operator-=(std::ptrdiff_t n) {
            it -= n;
            return *this;
        }

        // Overload equality operator (it1 == it2)
        bool operator==(const Iterator& other) const {
            return it == other.it;
        }

        // Overload inequality operator (it1 != it2)
        bool operator!=(const Iterator& other) const {
            return it != other.it;
        }

        // Overload - operator (it1 - it2)
        std::ptrdiff_t operator-(const Iterator& other) const {
            return it - other.it;
        }

        // Optionally, overload dereference operator (*it)
        T operator*() {
            // std::cout<<"iterator val: "<<*it
            // <<" "<<sa.size()<<std::endl;

            return sa.GetValForIterator(*it);
            // return *it;
        }

        T operator[](std::ptrdiff_t index) const{
            return sa.GetValForIterator(*(it + index));
            // return *(it + index);
        }
    };

    // Overloading the [] operator
    T operator[](std::size_t index) {
        if (index < 0 ||  index >= sa.size()) {
            throw std::out_of_range("Index out of range");
        }
        return GetKmerValAtSAIndx(index);
    }

    Iterator begin() {
        return Iterator(sa.begin(), *this);
    }

    Iterator end() {
        return Iterator(sa.end(), *this);
    }

    // Additional member functions for the vector class
    void push_back(const T& value) {
        sa.push_back(value);
    }

    std::size_t size() const {
        return sa.size();
    }

    const int64_t GetKmerValAtSAIndx(const int64_t i){
        return GetKmerVal(GetKmerAtIndex(i));
    }

    const int64_t get_largest(){
        uint64_t last_idx = fSize-1;
        while(GetKmerValAtSAIndx(last_idx) == -1) last_idx--;
        return GetKmerValAtSAIndx(last_idx);
    }

    inline const T get_sa_val(uint64_t idx) const {
        return sa[idx];
    }

    inline const int64_t GetKmerVal(const string &s){
        // if we have a smaller suffix than k-mer size,
        // we want to ignore it so that it never becomes a breakpoint
        if(s.compare("-1") == 0) return -1; 
        int64_t  val=0;
        int i=0;
        while(s[i]){
            val=val<<2;
            if(s[i]=='A'){
                val=val|0;
            }
            else if(s[i]=='C'){
                val=val|1;
            }
            else if(s[i]=='G'){
                val=val|2;
            }
            else{
                val=val|3;
            }
            i++;
        }
        return val;
    }

    const std::string GetKmerAtIndex(const int64_t i){
        int64_t suff_i_start = sa[i];
        if(str.begin()+suff_i_start+kmer_size <= str.end()){
            std::string suffix_i(str.begin()+suff_i_start, str.begin()+suff_i_start+kmer_size);
            return suffix_i;
        }
        return "-1"; // substring size < k: ignore these
    }

    const int64_t GetValForIterator(int64_t i){
        return GetKmerVal(GetKmerAtSAIndex_It(i));
    }
    
    
    const int64_t BinarySearch(const std::string &s, int64_t lo, int64_t hi, bool isFirstIndexReturned) const
    {
        int64_t mid, idx, nLcp, gLcp, loLcp = 0, hiLcp = 0;
        char ch = ' ';
        while(1){
            // nLcp = min(loLcp, hiLcp);  
            // cout<<lo<<" "<<mid<<" "<<hi<<endl;
            // std::cout<<GetKmerValAtSAIndx(lo)<<" "
            // <<GetKmerVal(s)<<" "
            // <<GetKmerValAtSAIndx(hi)<<endl;
            nLcp = loLcp < hiLcp ? loLcp : hiLcp;      
            if(hi - lo <= 1){
                for(int64_t i=lo; i<=hi; i++){
                    // sa_bs_cnt++;
                    idx = sa[i];
                    gLcp = nLcp;
                    for(; gLcp<kmer_size && idx+gLcp<fSize; gLcp++ ) {
                        ch = str[idx+gLcp];
                        if(s[gLcp] != ch) break;
                    }
                    if(gLcp == kmer_size) {
                        if(isFirstIndexReturned) return i;
                        return idx;
                    }
                }
                return -1; // not found
            }
            mid = (lo + hi) >> 1;
            idx = sa[mid]; 
            for(; nLcp<kmer_size && idx+nLcp<fSize; nLcp++ ) {
                ch = str[idx+nLcp];
                if(s[nLcp] != ch) break;
            }
            if(nLcp == kmer_size) {
                if(isFirstIndexReturned) return mid;
                return idx;
            }
            else if(nLcp + idx == fSize || s[nLcp] > ch)
            {
                lo = mid;
                loLcp = nLcp;
            }
            else
            {
                hi = mid;
                hiLcp = nLcp;
            }   
        }
    }

    std::string GetKmerAtSAIndex_It(int64_t suff_i_start){
        // int64_t suff_i_start = sa[i];
        if(str.begin()+suff_i_start+kmer_size <= str.end()){
            std::string suffix_i(str.begin()+suff_i_start, str.begin()+suff_i_start+kmer_size);
            return suffix_i;
        }
        return "-1"; // substring size < k: ignore these
    }

    std::string GetKmerAtStrPos(size_t strPos){
        if(strPos + kmer_size >= str.size()){
            // cout<<strPos<<" " <<saq.str.size()<<endl;
            return "-1";
        }
        if(strPos + kmer_size < str.size()){
            std::string kmer(str.begin()+strPos, str.begin()+strPos+kmer_size);
            return kmer;
        }
        std::string kmer(str.begin()+strPos, str.end());
        return kmer;
    }   

    void LoadGenomeString(std::string genomeFile){
        str.resize(fSize);
        FILE *fp = fopen(genomeFile.c_str(), "rb");
        //Last character is not handled -> should be $ or \0?
        int err = fread(&str[0], 1, fSize, fp);
        // str[fSize-1] = '\0';
        fclose(fp);
    }

    void LoadSA(std::string sa_file){
        sa.resize(fSize);
        FILE *fp = fopen(sa_file.c_str(), "rb");
        int err = fread(&sa[0], sizeof(int64_t), fSize, fp);
        fclose(fp);
    }

    void Load(std::string genomeF, std::string saF, int64_t km_size){
        std::ifstream in(genomeF.c_str(), std::ifstream::ate | std::ifstream::binary);
        fSize = in.tellg();
        printf("#Bases: %ld\n",fSize);
        this->LoadGenomeString(genomeF);        
        this->LoadSA(saF);
        kmer_size = km_size;
        printf("SA loaded\n");
    }
};