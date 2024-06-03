/**
 * inspired from https://github.com/ksahlin/strobealign/blob/main/src/cmdline.cpp
*/

#include "cmdline.hpp"
#include "args.hxx"

CommandLineOptions parse_command_line_arguments_build(int argc, char **argv) {
    args::ArgumentParser parser("build_pla_index" );
    parser.helpParams.showTerminator = false;
    parser.helpParams.helpindent = 20;
    parser.helpParams.width = 90;
    parser.helpParams.programName = "build_pla_index";
    parser.helpParams.shortSeparator = " ";

    args::HelpFlag help(parser, "help", "Print help and exit", {'h', "help"});

    // build_index
    args::ValueFlag<std::string>genome_fasta(parser, "STRING", "Fasta file with one entry and only ACGT characters. [Required]", {'g', "genome_fasta"},args::Options::Required);
    args::ValueFlag<std::string>suffix_array(parser, "STRING", "Suffix array of the genome in a binary file. [Required]", {'s', "suffix_array"},args::Options::Required);
    args::ValueFlag<int64_t> kmer_size(parser, "INT", "Kmer size to be used to construct the index. [default: 21]", {'k', "kmer_size"});
    args::ValueFlag<int64_t> eps(parser, "INT", "Epsilon value to be used for constructing the pla-index. [default: 15]", {'e', "eps"});
    args::ValueFlag<std::string>index_name(parser, "STRING", "File name where to save the index. [default: genome_fasta.index]", {'o', "index"});
    args::ValueFlag<int64_t> lookup_count(parser, "INT", "On average on how many elements the binary search on X array will take place. Used to determine the prefix lookup table. [default: 16]", {'l', "lookup"});
    
    try {
        parser.ParseCLI(argc, argv);
    }
    catch (const args::Completion& e) {
        std::cout << e.what();
        exit(EXIT_SUCCESS);
    }
    catch (const args::Help&) {
        std::cout << parser;
        exit(EXIT_SUCCESS);
    }
    catch (const args::Error& e) {
        std::cerr << parser;
        std::cerr << "Error: " << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }

    CommandLineOptions opt;

    if(genome_fasta) {opt.gn_fn = args::get(genome_fasta);}
    if(suffix_array) {opt.sa_fn = args::get(suffix_array);}
    if(kmer_size) {opt.kmer_size = args::get(kmer_size);}
    if(eps) {opt.eps = args::get(eps);}
    if(index_name) {opt.indx_fn = args::get(index_name);}
    if(lookup_count) {opt.lookup_count = args::get(lookup_count);}
    

    if(opt.indx_fn == "-1"){
        opt.indx_fn = opt.gn_fn+".index";
    }

    return opt;
}

CommandLineOptions parse_command_line_arguments_query(int argc, char **argv) {

    args::ArgumentParser parser("query_pla_index" );
    parser.helpParams.showTerminator = false;
    parser.helpParams.helpindent = 20;
    parser.helpParams.width = 90;
    parser.helpParams.programName = "query_pla_index";
    parser.helpParams.shortSeparator = " ";

    args::HelpFlag help(parser, "help", "Print help and exit", {'h', "help"});

    // build_index
    args::ValueFlag<std::string>genome_fasta(parser, "STRING", "Fasta file with one entry and only ACGT characters [Required]", {'g', "genome_fasta"},args::Options::Required);
    args::ValueFlag<std::string>suffix_array(parser, "STRING", "Suffix array of the genome in a binary file [Required]", {'s', "suffix_array"},args::Options::Required);
    args::ValueFlag<std::string>index_name(parser, "STRING", "Name of the stored index file [Required]", {'i', "index"},args::Options::Required);
    args::ValueFlag<std::string>query_file(parser, "STRING", "Name of a text file containing queries. Each line of the file should contain a single k-mer. The length of the kmer should be the same as the one provided while building the index. [Required]", {'q', "query"},args::Options::Required);
    args::ValueFlag<std::string>run_info(parser, "STRING", "File name to store query times", {"run_info"}, args::Options::Hidden);
        

    try {
        parser.ParseCLI(argc, argv);
    }
    catch (const args::Completion& e) {
        std::cout << e.what();
        exit(EXIT_SUCCESS);
    }
    catch (const args::Help&) {
        std::cout << parser;
        exit(EXIT_SUCCESS);
    }
    catch (const args::Error& e) {
        std::cerr << parser;
        std::cerr << "Error: " << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }

    CommandLineOptions opt;

    if(genome_fasta) {opt.gn_fn = args::get(genome_fasta);}
    if(suffix_array) {opt.sa_fn = args::get(suffix_array);}
    if(index_name) {opt.indx_fn = args::get(index_name);}

    if(query_file) {opt.query_file = args::get(query_file);}
    
    if(run_info) {opt.run_info_fn = args::get(run_info);}

    return opt;
}
