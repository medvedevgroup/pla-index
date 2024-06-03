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
#include "piecewise_linear_model.hpp"
#include "utils/ef_sequence.hpp"
#include "utils/essentials.hpp"
#include "../sdsl/dac_vector.hpp"
#include "utils/compact_vector.hpp"
#include "include/pthash.hpp"


using namespace std;

#define KNOT_BS_THRES 64
#define DEBUG_PRINT_QUERY 0

class pla_index{
private:
    using canonical_segment = typename OptimalPiecewiseLinearModel<int64_t, uint64_t>::CanonicalSegment;
    int64_t epsilon, max_data_ratio, lp_bits;
    uint64_t shift_bits, index_size, total_data_points, diff_bits;
    
    arank::compact_vector xy_cv, prefix_lookup_cv, err_dic_cv;

    arank::ef_sequence<false> y_range_beg;
    sdsl::dac_vector_dp<> sa_diff_dac_vec;
    
    

    typedef pthash::single_phf<pthash::murmurhash2_128,         // base hasher
                       pthash::dictionary_dictionary,  // encoder type
                       true                    // minimal
                       >
        pthash_type;
    pthash_type f;

public:
    // constructor at index build
    pla_index(int64_t eps, int64_t lookup, uint64_t largest,
            uint64_t total_points):
            epsilon(eps),
            lp_bits(lookup),
            total_data_points(total_points)
    {
        // lp_bits = lookup for now
        // later [on prefix calculation] shift_bits -= lp_bits
        if(ceil(log2(largest)) == floor(log2(largest))) shift_bits = ceil(log2(largest)) + 1;
        else shift_bits = ceil(log2(largest));
    }
    
    pla_index(){};
    
    inline uint64_t get_num_segments() const {return index_size;}
    
    template <typename DataType>
    void Save(string index_fn, DataType& data){
        essentials::save(f, index_fn.c_str());
        std::ofstream os(index_fn, std::ios_base::app);
        uint8_t t_lp_bits = lp_bits;
        essentials_arank::save_pod(os, t_lp_bits);
        uint8_t kmer_size = data.get_kmer_size();
        essentials_arank::save_pod(os, kmer_size);
        prefix_lookup_cv.save(os);
        uint32_t nBrkpnts = index_size;
        essentials_arank::save_pod(os, nBrkpnts);
        xy_cv.save(os);
        y_range_beg.save(os);
        uint16_t t_err_th = epsilon;
        essentials_arank::save_pod(os, t_err_th);
        // uint64_t uniq_kmers =err_dic_cv.size();
        // essentials_arank::save_pod(os, uniq_kmers);
        err_dic_cv.save(os);
    }

    template <typename DataType>
    void Load(string index_fn, DataType &data){
        uint64_t total_bytes = essentials::load(f, index_fn.c_str());
        std::ifstream is(index_fn, std::ios::binary);
        is.seekg(total_bytes, std::ios::beg);
        uint8_t lp_bits_t;
        essentials_arank::load_pod(is, lp_bits_t);
        lp_bits = lp_bits_t;
        // cout<<"lp bits: "<<lp_bits<<endl;
       
        uint8_t kmer_size;
        essentials_arank::load_pod(is, kmer_size);
        // cout<<"kmer size: "<<uint64_t(kmer_size)<<endl;
        data.set_kmer_size(kmer_size);
        total_data_points = data.size();
        // cout<<"Total points: "<<total_data_points<<endl;
        uint64_t largest = data.get_largest();
        // cout<<"largest value: "<<largest<<endl;
        uint64_t total_bits;
        if(ceil(log2(largest)) == floor(log2(largest))) total_bits = ceil(log2(largest)) + 1;
        else total_bits = ceil(log2(largest));
        shift_bits = total_bits - lp_bits;
        // cout<<"total: "<<total_bits<<endl;
       
        prefix_lookup_cv.load(is);
        // cout<<"prefix size: "<<prefix_lookup_cv.size()<<endl;
        uint32_t index_size_t;
        essentials_arank::load_pod(is, index_size_t);
        index_size = index_size_t;
        max_data_ratio = ((1ULL << total_bits) - 1)/(index_size-1);
        if(DEBUG_PRINT)
            cout<<"#breakpoints: "<<index_size<<endl;
        xy_cv.load(is);
        y_range_beg.load(is);
        uint16_t max_error_t;
        essentials_arank::load_pod(is, max_error_t);
        epsilon = int64_t(max_error_t);
        diff_bits = (ceil(log2(epsilon*2)) == floor(log2(epsilon*2)))
                ? ceil(log2(epsilon*2)) + 1 + 1
                : ceil(log2(epsilon*2)) + 1;
        err_dic_cv.load(is);

        // Save_temp(index_fn, data);
    }
    
    // suffix array version
    template <typename DataType>
    int64_t query(const int64_t query_val, const DataType &data){
        int64_t query_pred = get_query_pred(query_val, data);
        int64_t true_idx = query_pred + get_signed_value(err_dic_cv[f(query_val)]);
        // if(DEBUG_PRINT_QUERY){
        //     cout<<"\nTrue idx: "<<true_idx<<" pred: "<<query_pred
        //         <<" err: "<<err_dic_cv[f(query_val)]<<endl;
        //     cout<<"idx+1 val: "<<data[true_idx+1]
        //         <<" idx val: "<<data[true_idx]
        //         <<" idx-1 val: "<<data[true_idx-1]<<endl;
        // }
        return data.get_sa_val(true_idx);
    }

    inline int64_t get_signed_value(int64_t value) const{
        return value & 1
            ? -(value>>1)
            : value>>1;
    }

    inline int64_t predict(int64_t q, int64_t x1, int64_t y1, int64_t x2_x1, int64_t y2) const{
        return round(y1 + ((double)(q - x1)/x2_x1 )*(y2 - y1));
    }

    template <typename DataType>
    int64_t get_query_pred(const int64_t query_val, const DataType &data){
        int64_t query_pred;        
        pair<pair<int64_t, int64_t>, pair<int64_t, int64_t> > brkpnt_pair;
        uint64_t lp_bit_value = query_val >> shift_bits;
        const uint64_t diff_mask = (1ULL << diff_bits) - 1;
        // cout<<"lp: "<<lp_bit_value<<" diff_bits: "<<diff_bits
        //     <<" Query: "<<query_val
        //     <<endl;
        int64_t ind_brk_beg_idx = prefix_lookup_cv[lp_bit_value];
        
        int64_t ind_point_val = get_signed_value(xy_cv[ind_brk_beg_idx] >> diff_bits) 
            + max_data_ratio*ind_brk_beg_idx;
        
        // if(DEBUG_PRINT_QUERY){
        //     cout<<"lp: "<<lp_bit_value<<" diff_bits: "<<diff_bits
        //         <<" ind_beg_indx: "<<ind_brk_beg_idx
        //         <<" x[ind]: "<<ind_point_val
        //         <<" query: "<<query_val
        //         <<endl;
        // }

        if(ind_point_val == query_val){
            query_pred = y_range_beg.access(ind_brk_beg_idx) - epsilon;
        }
        else{
            if(ind_point_val < query_val){
                brkpnt_pair = xy_cv.binary_search(ind_brk_beg_idx, 
                    prefix_lookup_cv[lp_bit_value + 1] - 1, query_val, max_data_ratio, diff_bits, diff_mask);
            }
            else{
                auto [x1, y_diff] = xy_cv.get_xy(ind_brk_beg_idx-1, 
                    max_data_ratio,  diff_bits, diff_mask);
                brkpnt_pair = make_pair(make_pair(x1, ind_point_val - x1), 
                    make_pair(y_diff, ind_brk_beg_idx-1));
            }
            auto [val1, val2] = y_range_beg.pair(brkpnt_pair.second.second);
            query_pred = predict(query_val, brkpnt_pair.first.first, val1, 
                    brkpnt_pair.first.second, val2 - get_signed_value(brkpnt_pair.second.first)) - epsilon;            

        }
        return std::clamp(query_pred, int64_t(0), int64_t(total_data_points) - 1);
    }

    template<typename RandomIt>
    int64_t get_pred_diff(const int64_t act_idx, int64_t pred_idx, 
        RandomIt data) const{
        if(pred_idx > total_data_points) pred_idx = total_data_points-1;
        if(act_idx == pred_idx){
            return 0;
        }
        int64_t lo, hi, mid, act_data = data[act_idx];
        if(act_idx < pred_idx){
            lo = act_idx;
            hi = pred_idx;
            while (1){
                // cout<<hi<<" "<<lo<<" "<<mid<<endl;
                if(hi - lo <= 1){
                    if(hi == lo) return hi - pred_idx;
                    if(data[hi] == act_data) return hi - pred_idx;
                    // if(getLcp(sa[hi], act_data, 0) == kmer_size) return hi - pred_idx;
                    // lo is guranteed to be same as act_idx as never updated otherwise
                    return lo - pred_idx; 
                }
                mid = (lo + hi) >> 1;
                if(data[mid] == act_data) lo = mid;
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
                if(data[lo] == act_data) return lo - pred_idx;
                // if(getLcp(sa[lo], act_data, 0) == kmer_size) return lo - pred_idx;
                // hi is guranteed to be same as act_idx as never updated otherwise
                return hi - pred_idx; 
            }
            mid = (lo + hi) >> 1;
            if(data[mid] == act_data) hi = mid;
            // if(getLcp(sa[mid], act_data, 0) == kmer_size) hi = mid;
            else lo = mid;
        }
    }

    // mphf version
    template<typename RandomIt>
    void build_index(const RandomIt &begin){ 
        build_basic_pla_index(begin);
        create_err_table(begin);
    }

    template<typename RandomIt>
    void create_err_table(const RandomIt &data){    
        int64_t indx = 0, curr_kval, prev_kval;
        int64_t query_pred, query_err;
        vector<int64_t> unique_kvals; //keys
        vector<int64_t> true_err; // to first position
        while(1){
            curr_kval = data[indx];
            // cout<<indx<<" "<<GetKmerAtIndex(indx)<<" "<<curr_kval<<endl;
            indx++;        
            if(curr_kval == -1) continue;
            unique_kvals.emplace_back(curr_kval);

            // query_pred_err = get_query_pred_err(indx-1, curr_kval);
            query_pred = get_query_pred(curr_kval, data);
            query_err = add_bit_based_on_sign((indx - 1) - query_pred);
            true_err.emplace_back(query_err);
            prev_kval = curr_kval;
            break;
        }
        
        for(int64_t curr_sa_indx=indx; curr_sa_indx<total_data_points; curr_sa_indx++){     
            // cout<<curr_sa_indx<<" "; 
            curr_kval = data[curr_sa_indx];
            if(curr_kval == prev_kval || curr_kval == -1) {
                // cout<<curr_kval<<endl;
                continue;
            }        
            // cout<<"ok ";
            unique_kvals.emplace_back(curr_kval);
            query_pred = get_query_pred(curr_kval, data);
            query_err = add_bit_based_on_sign(curr_sa_indx - query_pred);
            // cout<<curr_kval<<" ";
            
            true_err.emplace_back(query_err);
            if(DEBUG_PRINT && query_err > 2*epsilon){
                cout<<"query_pred: "<<query_pred
                <<" query_err: "<<query_err
                <<" curr_sa_idx: "<<curr_sa_indx
                <<" curr kval: "<<curr_kval
                <<endl;
                
            }
            prev_kval = curr_kval;
            
        }
        build_mphf(unique_kvals, true_err);   
    }

    void build_mphf(vector<int64_t>& keys, vector<int64_t>& true_err){
        pthash::build_configuration config;
        config.c = 5.0;
        config.alpha = 0.94;
        config.minimal_output = true;  // mphf
        config.verbose_output = false;
        
        auto start = pthash::clock_type::now();
        auto timings = f.build_in_internal_memory(keys.begin(), keys.size(), config);
        double total_seconds = timings.partitioning_seconds + timings.mapping_ordering_seconds +
                            timings.searching_seconds + timings.encoding_seconds;
        
        /* Compute and print the number of bits spent per key. */
        double bits_per_key = static_cast<double>(f.num_bits()) / f.num_keys();
        // mphf_time << "function uses " << bits_per_key << " [bits/key]" << std::endl;
        if(DEBUG_PRINT){
            cout << "function built in " << pthash::seconds(pthash::clock_type::now() - start) << " seconds"
                << std::endl;
            cout << "function uses " << bits_per_key << " [bits/key]" << std::endl;
        }

        // storing the err values in the index from MPHF and Bitpack it
        vector<int64_t>err_dict(keys.size());
        err_dict.resize(keys.size());
        int64_t _max_err = 0;
        for(uint64_t i=0 ;i<keys.size(); i++){
            err_dict[f(keys[i])] = true_err[i];
            if(true_err[i] > _max_err) _max_err = true_err[i];
        }
        
        err_dic_cv.build(err_dict.begin(), err_dict.size());
    }

    void encode_values(const vector<uint64_t> &brk_kval_vec, const vector<uint64_t> &brk_diff_vec){
        vector<uint64_t> encoded_values_vec;
        uint64_t max_bits = shift_bits + lp_bits;
        
        max_data_ratio = ((1ULL << (max_bits)) - 1)/(index_size - 1);
        uint64_t max_value = 0;
        
        for(int64_t i=0; i<index_size; i++){
            
            int64_t diff = int64_t(brk_kval_vec[i]) - i*max_data_ratio;
            encoded_values_vec.emplace_back(add_bit_based_on_sign(diff));
            if(encoded_values_vec[i] > max_value) max_value = encoded_values_vec[i];
        }
        
        auto signed_bit_count = [](int64_t n) -> int64_t {
            return (ceil(log2(n)) == floor(log2(n)))
                ? ceil(log2(n)) + 1 + 1
                : ceil(log2(n)) + 1;
        };

        
        diff_bits = signed_bit_count(2*epsilon);
        const uint64_t diff_mask = (1ULL << diff_bits) - 1;
        
        if(DEBUG_PRINT){
            cout<<"Encoder vec size: "<<encoded_values_vec.size()
                <<" y_diff vec size: "<<brk_diff_vec.size()
                <<endl;
            
        }

        uint64_t max_diff = 0;
        for(int64_t i=0; i<encoded_values_vec.size(); i++){
            uint64_t val = encoded_values_vec[i] << diff_bits;
            val = (i < encoded_values_vec.size()-1)
                ? val | (brk_diff_vec[i] & diff_mask)
                : val;
            if(i < encoded_values_vec.size()-1 && max_diff < brk_diff_vec[i]) 
                max_diff = brk_diff_vec[i];
            encoded_values_vec[i] = val;
        }

        xy_cv.build(encoded_values_vec.begin(), encoded_values_vec.size());
    }

    void insert_to_prefix_lookup(vector<uint64_t>& brk_kval_vec){
        vector<uint64_t> prefix_lookup_vec;
        lp_bits = ceil(log2((brk_kval_vec.size()/lp_bits)));
        shift_bits -= lp_bits; // on constructor, shift_bits = total_bits
        if(DEBUG_PRINT)
            cout<<"Lp bits: "<<lp_bits<<" shift bits: "<<shift_bits<<endl;
        
        
        for(uint64_t i=0; i<brk_kval_vec.size(); i++){
            int64_t lp_bit_value = brk_kval_vec[i] >> shift_bits;
            if(prefix_lookup_vec.size() > lp_bit_value) continue;
            for(int64_t pl=prefix_lookup_vec.size(); pl<lp_bit_value+1; pl++){
                prefix_lookup_vec.emplace_back(i);
            }
        }
        uint64_t prefix_lookup_max = 1ULL << lp_bits;
        
        for(int64_t pl = prefix_lookup_vec.size(); pl < prefix_lookup_max; pl++){
            prefix_lookup_vec.emplace_back(brk_kval_vec.size()-1); 
        }
        // extra element to ease query
        prefix_lookup_vec.emplace_back(brk_kval_vec.size()-1); 
        prefix_lookup_cv.build(prefix_lookup_vec.begin(), prefix_lookup_vec.size());
    }

    /**
     * @brief add a 1 bit at the end of a negative number, and 
     * 0 at the end of a positive number
    */
    inline uint64_t add_bit_based_on_sign(int64_t num){
        return num < 0
            ? ((-num)<<1) | 1
            : num << 1;
    }

    template<typename RandomIt>
    void build_basic_pla_index(const RandomIt &begin){
        vector<uint64_t>  brk_sa_indx_vec, brk_kval_vec, brk_diff_vec;
        bool isFirst = true;
        int64_t prev_knot_ei=0, diff_from_curr_start;
        vector<canonical_segment> out(1);

        auto in_fun = [&](auto i) { return std::pair<int64_t, uint64_t>(begin[i], i); };
        int64_t _max = 0, bit_max = 0;
        auto out_fun = [&](auto cs, const int64_t last_x) { 
            out[0] = cs;
            brk_kval_vec.emplace_back(cs.get_first_x());
            // insert_to_ind_vec(cs.get_first_x(), brk_kval_vec.size()-1);
            
            auto [knot_si, knot_ei] = cs.get_knot_intersection(last_x, epsilon);  
            if(DEBUG_PRINT && cs.get_first_x() == 679337644716){
                cout<<"last x: "<<last_x
                    <<" knot_si: "<<knot_si
                    <<" knot_ei: "<<knot_ei
                    <<endl;
            }
            brk_sa_indx_vec.emplace_back((uint64_t)knot_si); 
            if(isFirst){
                prev_knot_ei = (int64_t)knot_ei;
                isFirst = false;
                return;
            }
            diff_from_curr_start = (int64_t)knot_si - prev_knot_ei;
            if(abs(diff_from_curr_start) > _max) _max = abs(diff_from_curr_start);
            
            brk_diff_vec.emplace_back(add_bit_based_on_sign(diff_from_curr_start));
            if(brk_diff_vec[brk_diff_vec.size()-1] > bit_max) 
                bit_max = brk_diff_vec[brk_diff_vec.size()-1];
            prev_knot_ei = (int64_t)knot_ei;
        };
        auto get_err = [&](int64_t act_idx, int64_t pred_idx){
            return act_idx - pred_idx;
        };

        uint64_t seg_count = make_segmentation_basic_pla(total_data_points, 
                    epsilon, in_fun, out_fun, get_err);

        brk_kval_vec.emplace_back(out[out.size()-1].get_last_x());
        
        insert_to_prefix_lookup(brk_kval_vec); //loopup count = lp_bits
        index_size = brk_kval_vec.size();

        auto [knot_si, knot_ei] = out[0].get_knot_intersection(out[0].get_last_x(), epsilon);
        brk_sa_indx_vec.emplace_back((uint64_t)knot_ei);
        brk_diff_vec.emplace_back(0);
        if(DEBUG_PRINT)
            cout<<"Diff in sa: "<<_max<<" "<<bit_max<<endl;
        
        y_range_beg.encode(brk_sa_indx_vec.data(), brk_sa_indx_vec.size());
        encode_values(brk_kval_vec, brk_diff_vec);
    }
    

};