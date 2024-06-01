/**
 * inspired from https://github.com/ksahlin/strobealign/blob/main/src/cmdline.hpp
*/

#pragma once

#include <vector>
#include <string>
#include <utility>

struct CommandLineOptions {
    std::string gn_fn;
    std::string sa_fn;

    int64_t kmer_size {21};
    int64_t eps {15};
    std::string indx_fn {"-1"};
    int64_t lookup_count {16};
    std::string indx_type {"basic-pla"};
    bool is_fast_rank {false};
    std::string query_type {"search"};
    std::string query_file {"-1"};
    std::string run_info_fn {"-1"};
};

CommandLineOptions parse_command_line_arguments_query(int argc, char **argv);
CommandLineOptions parse_command_line_arguments_build(int argc, char **argv);

