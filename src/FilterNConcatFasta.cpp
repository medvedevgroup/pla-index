#include<bits/stdc++.h>
#include <iostream>
#include <fstream>
using namespace std;

// reference: https://stackoverflow.com/a/58090517
void strip(std::string &str)
{   
    if  (str.length() != 0)
    {   
        auto w = std::string(" ") ;
        auto n = std::string("\n") ;
        auto r = std::string("\t") ;
        auto t = std::string("\r") ;
        auto v = std::string(1 ,str.front()); 
        while((v == w) || (v==t) || (v==r) || (v==n))
        {   
            str.erase(str.begin());
            v = std::string(1 ,str.front());
        }
        v = std::string(1 , str.back()); 
        while((v ==w) || (v==t) || (v==r) || (v==n))
        {   
            str.erase(str.end() - 1 );
            v = std::string(1 , str.back());
        }
    }
}

void FilterNFromFile(string inFileName, string outFileName)
{   
    ofstream outFile(outFileName.c_str());
    ifstream inFile(inFileName.c_str());
    string line;
    string filteredLine = "";

    while (getline (inFile, line)){
        strip(line);
        filteredLine = "";
        for(int i=0;i<line.size(); i++){
            if(line[i] == '>'){
                break;
            }
            if(line[i] != 'N' && line[i] != 'n'){
                if(line[i] >= 'a' && line[i] <='z') line[i] -= ('a' - 'A');
                if(line[i] =='A' || line[i] =='C' ||
                    line[i] =='G' ||line[i] =='T'){
                    filteredLine += line[i];
                }
            }
        }
        outFile<<filteredLine;
    }

    outFile.close();
    inFile.close();
}   

int main(int argc, char *argv[]){
    string inFile = argv[1]; //fasta file
    string fn_pref =argv[2]; // same as folder name for snakemake
    strip(inFile);
    string outFile = fn_pref+".ConcatenatedGenome.txt";

    FilterNFromFile(inFile, outFile);

    return 0;
}