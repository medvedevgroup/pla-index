//  Created by Kristoffer Sahlin on 4/21/21.

#ifndef STROBEALIGN_INDEX_HPP
#define STROBEALIGN_INDEX_HPP

#include <chrono>
#include <stdio.h>
#include <string>
#include <vector>
#include <deque>
#include <tuple>
#include <cmath>
#include <iostream>
#include <cassert>
#include "robin_hood.h"
#include "exceptions.hpp"
#include "refs.hpp"
#include "randstrobes.hpp"
#include "indexparameters.hpp"
#include "logger.hpp"
#include "pla_index.hpp"


struct IndexCreationStatistics {
    uint64_t tot_strobemer_count = 0;
    uint64_t tot_occur_once = 0;
    uint64_t tot_high_ab = 0;
    uint64_t tot_mid_ab = 0;
    uint64_t index_cutoff = 0;
    uint64_t filter_cutoff = 0;
    uint64_t distinct_strobemers = 0;

    std::chrono::duration<double> elapsed_hash_index;
    std::chrono::duration<double> elapsed_generating_seeds;
    std::chrono::duration<double> elapsed_counting_hashes;
    std::chrono::duration<double> elapsed_sorting_seeds;
};

struct StrobemerIndex {
    using bucket_index_t = uint64_t;
    StrobemerIndex(const References& references, const IndexParameters& parameters, int bits=-1, int eps)
        : filter_cutoff(0)
        , parameters(parameters)
        , references(references)
        , bits(bits == -1 ? pick_bits(references.total_length()) : bits)  
        , rs_pla_index(eps, 16)
    {
        if (this->bits < 8 || this->bits > 31) {
            throw BadParameter("Bits must be between 8 and 31");
        }
    }
    unsigned int filter_cutoff;
    mutable IndexCreationStatistics stats;

    void write(const std::string& filename) const;
    void read(const std::string& filename);
    void populate(float f, size_t n_threads);

    void check();

    void print_diagnostics(const std::string& logfile_name, int k) const;
    int pick_bits(size_t size) const;
    size_t find(randstrobe_hash_t key) const;
    

    randstrobe_hash_t get_hash(bucket_index_t position) const {
        if (position < randstrobes.size()) {
            return randstrobes[position].hash;
        } else {
            return end();
        }
    }
    
    bool is_filtered(bucket_index_t position) const {
        return get_hash(position) == get_hash(position + filter_cutoff);
    }

    unsigned int get_strobe1_position(bucket_index_t position) const {
        return randstrobes[position].position;
    }

    int strobe2_offset(bucket_index_t position) const {
        return randstrobes[position].strobe2_offset();
    }

    int reference_index(bucket_index_t position) const {
        return randstrobes[position].reference_index();
    }

    RefRandstrobe get_randstrobe(bucket_index_t position) const {
        return randstrobes[position];
    }

    size_t size() const {
        return randstrobes.size();
    }

    unsigned int get_count(bucket_index_t position) const;

    size_t end() const {
        return -1;
    }

    int k() const {
        return parameters.syncmer.k;
    }

    int get_bits() const {
        return bits;
    }

private:
    void assign_all_randstrobes(const std::vector<uint64_t>& randstrobe_counts, size_t n_threads);
    void assign_randstrobes(size_t ref_index, size_t offset);

    const IndexParameters& parameters;
    const References& references;

    /*
     * The randstrobes vector contains all randstrobes sorted by hash.
     *
     * The randstrobe_start_indices vector points to entries in the
     * randstrobes vector. randstrobe_start_indices[x] is the index of the
     * first entry in randstrobes whose top *bits* bits of its hash value are
     * greater than or equal to x.
     *
     * randstrobe_start_indices has one extra guard entry at the end that
     * is always randstrobes.size().
     */

    std::vector<RefRandstrobe> randstrobes;
    // std::vector<bucket_index_t> randstrobe_start_indices;
    int bits; // no. of bits of the hash to use when indexing a randstrobe bucket
    pla_index rs_pla_index;
};

#endif
