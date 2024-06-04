#pragma once

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <string>
#include <cstdlib>
#include <chrono>
#include <cassert>
#include "BitPacking.h"
#include "piecewise_linear_model.hpp"
#include "utils/ef_sequence.hpp"
#include "utils/essentials.hpp"
#include "../../sdsl/dac_vector.hpp"

#include "randstrobes.hpp"
#include "logger.hpp"  

using namespace std;

class pla_index{
private:
    // x value never negative, uses 64 bit hash
    using canonical_segment = typename OptimalPiecewiseLinearModel<uint64_t, uint64_t>::CanonicalSegment;
    int64_t epsilon, lp_bits, knot_bs_thres;
    uint64_t shift_bits, dic_size, total_data_points, max_data_ratio;
    
    vector<uint64_t> indirection_vec, index_bitvec;

    CBitPacking brk_uniform_diff_packed;
    arank::ef_sequence<false> y_range_beg;
    sdsl::dac_vector_dp<> sa_diff_dac_vec;

public:
    pla_index(int64_t eps, int64_t lp_bits);
    
    pla_index(int64_t knot_bs_thres);

    void insert_to_ind_vec(uint64_t brk_beg_kval, int64_t brk_indx);

    void encode_knots(const vector<uint64_t> &brk_kval_vec);    
    
    void SetBit(uint64_t &number, uint64_t pos);
    void Save(std::ostream& os) const;
    void save_unpacked(string dic_fn) const;
    void save_bit_info(string dic_fn, CBitPacking &ind_diff_pack) const;
    const vector<uint64_t> get_processed_ind_vec() const;
    void Load(std::istream& is, const vector<RefRandstrobe> &randstrobes);
    int64_t binary_search(int64_t lo, int64_t hi, const uint64_t kval) const;
    uint64_t BinarySearch_BrkPnt(uint64_t lo, uint64_t hi, const uint64_t query) const;
    
    inline int64_t get_epsilon()const{return epsilon;}

    int64_t query(const vector<RefRandstrobe> &randstrobes, const uint64_t query_val) const;
    uint64_t QueryKmer(const vector<RefRandstrobe> &randstrobes, const uint64_t query_val) const;
    
    int64_t get_pred_diff(const int64_t act_idx, int64_t pred_idx, 
        const vector<RefRandstrobe> &randstrobes) const{
        if(pred_idx > total_data_points) pred_idx = total_data_points-1;
        if(act_idx == pred_idx){
            return 0;
        }
        int64_t lo, hi, mid;
        uint64_t act_data = randstrobes[act_idx].hash;
        if(act_idx < pred_idx){
            lo = act_idx;
            hi = pred_idx;
            while (1){
                // cout<<hi<<" "<<lo<<" "<<mid<<endl;
                if(hi - lo <= 1){
                    if(hi == lo) return hi - pred_idx;
                    if(randstrobes[hi].hash == act_data) return hi - pred_idx;
                    // if(getLcp(sa[hi], act_data, 0) == kmer_size) return hi - pred_idx;
                    // lo is guranteed to be same as act_idx as never updated otherwise
                    return lo - pred_idx; 
                }
                mid = (lo + hi) >> 1;
                if(randstrobes[mid].hash == act_data) lo = mid;
                // if(getLcp(sa[mid], act_data, 0) == kmer_size) lo = mid;
                else hi = mid;
            }
        }
        lo = pred_idx;
        hi = act_idx;
        while(1){
            // cout<<hi<<" "<<lo<<" "<<mid<<endl;
            if(hi - lo <= 1){
                if(hi == lo) return hi - pred_idx;
                if(randstrobes[lo].hash == act_data) return lo - pred_idx;
                // if(getLcp(sa[lo], act_data, 0) == kmer_size) return lo - pred_idx;
                // hi is guranteed to be same as act_idx as never updated otherwise
                return hi - pred_idx; 
            }
            mid = (lo + hi) >> 1;
            if(randstrobes[mid].hash == act_data) hi = mid;
            // if(getLcp(sa[mid], act_data, 0) == kmer_size) hi = mid;
            else lo = mid;
        }
    }
    
    void build_rep_stretch_pla_index_strobe(const vector<RefRandstrobe> &randstrobes);

};