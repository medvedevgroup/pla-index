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
#include "../sdsl/dac_vector.hpp"
#include "include/pthash.hpp"


using namespace std;

enum INDX_TYPE {BASIC_PLA, REPEAT_PLA, EXACT_PLA};

class pla_index{
private:
    using canonical_segment = typename OptimalPiecewiseLinearModel<int64_t, uint64_t>::CanonicalSegment;
    int64_t epsilon, max_data_ratio, lp_bits, knot_bs_thres;
    uint64_t shift_bits, index_size, total_data_points;
    
    vector<uint64_t> indirection_vec;

    CBitPacking brk_uniform_diff_packed, diff_packd, err_dict_packed;
    arank::ef_sequence<false> y_range_beg;
    sdsl::dac_vector_dp<> sa_diff_dac_vec;
    INDX_TYPE indx_type;

    typedef pthash::single_phf<pthash::murmurhash2_128,         // base hasher
                       pthash::dictionary_dictionary,  // encoder type
                       true                    // minimal
                       >
        pthash_type;
    pthash_type f;

public:
    pla_index(int64_t eps, int64_t lp_bits, uint64_t largest, uint64_t total_points, INDX_TYPE indx_type);
    
    pla_index(int64_t knot_bs_thres);

    void insert_to_ind_vec(int64_t brk_beg_kval, int64_t brk_indx);

    void encode_knots(const vector<int64_t> &brk_kval_vec);    
    
    void SetBit(uint64_t &number, uint64_t pos);
    void Save(string index_fn);
    void save_unpacked(string index_fn);
    void save_bit_info(string index_fn, CBitPacking &ind_diff_pack);
    uint64_t get_total_size_in_bytes();
    const vector<uint64_t> get_processed_ind_vec() const;
    void Load(string index_fn, int64_t total_bits);
    int64_t binary_search(int64_t lo, int64_t hi, const int64_t kval) const;
    inline uint64_t get_num_knots() const {return index_size;}
    
    // suffix array version
    template <typename DataType>
    int64_t query(const int64_t query_val, const string &query, const DataType &data) const{
        int64_t query_pred = get_query_pred(query_val, data);
        int64_t true_idx = query_pred + err_dict_packed.GetValueAt(f(query_val));
        return data.get_sa_val(true_idx);
    }

    // suffix array version
    template <typename DataType>
    int64_t get_query_pred(const int64_t query_val, const DataType &data) const{
        uint64_t lp_bit_value = query_val >> shift_bits;
        int64_t brkpnt, query_pred, direction;
        int64_t ind_brk_beg_idx = indirection_vec[lp_bit_value];
        bool isOnIntersection = false;

        int64_t ind_point_val = brk_uniform_diff_packed.GetValueAt(ind_brk_beg_idx)+max_data_ratio*ind_brk_beg_idx;        

        if(ind_point_val == query_val){
            query_pred = y_range_beg.access(ind_brk_beg_idx) - epsilon;
        }
        else{
            if(ind_point_val < query_val){
                int64_t end_idx = (lp_bit_value == indirection_vec.size() - 1) ? index_size - 1 : indirection_vec[lp_bit_value + 1] - 1;
                brkpnt = binary_search(ind_brk_beg_idx, end_idx, query_val);
            }
            else{
                do{
                    --lp_bit_value;
                    ind_point_val = indirection_vec[lp_bit_value];
                }while(ind_point_val == ind_brk_beg_idx);
                brkpnt = binary_search(ind_point_val, ind_brk_beg_idx, query_val);
            }
            // cout<<query_val<<" "<<brkpnt<<" "<<direction<<endl;
            const int64_t brk_beg_kval = brk_uniform_diff_packed.GetValueAt(brkpnt)+max_data_ratio*brkpnt;
            if(brk_beg_kval == query_val) {
                query_pred = y_range_beg.access(brkpnt) - epsilon;
            }
            else{
                const int64_t brk_end_kval = brk_uniform_diff_packed.GetValueAt(brkpnt+1)+max_data_ratio*(brkpnt+1);
                const int64_t brk_beg_sa_indx = int64_t(y_range_beg.access(brkpnt)) - epsilon;
                const int64_t brk_end_sa_indx = y_range_beg.access(brkpnt + 1) - 
                        diff_packd.GetValueAt(brkpnt) - epsilon;
                    
                // cout<<"Interval: \n"<<brk_beg_kval<<"\n"<<query_val<<"\n"<<brk_end_kval<<endl;
                // cout<<query_val<<endl;
                // cout<<brk_beg_kval<<" "
                //     <<brk_end_kval<<" "
                //     <<brk_beg_sa_indx<<" "
                //     <<brk_end_sa_indx<<endl;
                query_pred = round(brk_beg_sa_indx + 
                        ((double)(query_val - brk_beg_kval)/(brk_end_kval - brk_beg_kval) )*
                        (brk_end_sa_indx - brk_beg_sa_indx));
                // cout<<"Here, Pred: "<<query_pred<<endl;
                // std::cout<<"Pred: "<<query_pred<<" "<<GetKmerAtStrPos(sa[query_pred])<<endl;
            }  
            // cout<<"Second, Pred: "<<query_pred<<endl;          
        }
        // std::cout<<"After bracket Pred:  "<<query_pred<<endl;
        // <<" "<<GetKmerAtStrPos(sa[query_pred])<<endl;
        // query_pred -= epsilon;
        // if(query_pred < int64_t(0)) query_pred = int64_t(0);        
        // if(query_pred > int64_t(total_data_points)-1) query_pred = int64_t(total_data_points)-1;
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
            indx++;        
            if(curr_kval == -1) continue;
            unique_kvals.emplace_back(curr_kval);

            // query_pred_err = get_query_pred_err(indx-1, curr_kval);
            query_pred = get_query_pred(curr_kval, data);
            query_err = (indx - 1) - query_pred;
            true_err.emplace_back(query_err);
            prev_kval = curr_kval;
            break;
        } 
        for(int64_t curr_sa_indx=indx; curr_sa_indx<total_data_points; curr_sa_indx++){     
            // cout<<curr_sa_indx<<endl; 
            curr_kval = data[curr_sa_indx];
            if(curr_kval == prev_kval || curr_kval == -1) {
                // cout<<curr_kval<<endl;
                continue;
            }        
            // cout<<"ok ";
            unique_kvals.emplace_back(curr_kval);
            query_pred = get_query_pred(curr_kval, data);
            query_err = curr_sa_indx - query_pred;
            // cout<<curr_kval<<" ";
            true_err.emplace_back(query_err);
            // cout<<query_pred_err.first<<" "<<fSize<<endl;
            // pred_vec[query_pred_err.first]++;
            // cout<<"ok\n";
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
        // ofstream mphf_time("mphf_time"+to_string(epsilon)+".txt");
        // mphf_time << "function built in " << pthash::seconds(pthash::clock_type::now() - start) << " seconds"
                // << std::endl;
        cout << "mphf function built in " << pthash::seconds(pthash::clock_type::now() - start) << " seconds"
                << std::endl;
        // mphf_time << "computed: " << total_seconds << " seconds" << std::endl;
        /* Compute and print the number of bits spent per key. */
        double bits_per_key = static_cast<double>(f.num_bits()) / f.num_keys();
        // mphf_time << "function uses " << bits_per_key << " [bits/key]" << std::endl;
        cout << "mphf function uses " << bits_per_key << " [bits/key]" << std::endl;

        // storing the err values in the index from MPHF and Bitpack it
        vector<int64_t>err_dict(keys.size());
        err_dict.resize(keys.size());
        for(uint64_t i=0 ;i<keys.size(); i++){
            err_dict[f(keys[i])] = true_err[i];
            // if(keys[i] == 0){
            //     cout<<"mphf: "<<f(keys[i])<<" true err: "<<true_err[i]<<endl;
            // }
        }
        err_dict_packed.BuilidSignedPackedVector(err_dict);
    }

    template<typename RandomIt>
    void build_basic_pla_index(const RandomIt &begin){
        vector<uint64_t>  brk_sa_indx_vec, brk_sa_end_vec;
        vector<int64_t> brk_kval_vec, brk_diff_vec;
        bool isFirst = true;
        int64_t prev_knot_ei=0, diff_from_curr_start;
        vector<canonical_segment> out(1);

        auto in_fun = [&](auto i) { return std::pair<int64_t, uint64_t>(begin[i], i); };
            
        auto out_fun = [&](auto cs, const int64_t last_x) { 
            out[0] = cs;
            brk_kval_vec.emplace_back(cs.get_first_x());
            insert_to_ind_vec(cs.get_first_x(), brk_kval_vec.size()-1);
            
            auto [knot_si, knot_ei] = cs.get_knot_intersection(last_x, epsilon);      
            brk_sa_indx_vec.emplace_back((uint64_t)knot_si); 
            if(isFirst){
                prev_knot_ei = (int64_t)knot_ei;
                isFirst = false;
                return;
            }
            diff_from_curr_start = (int64_t)knot_si - prev_knot_ei;
            brk_diff_vec.emplace_back(diff_from_curr_start);
            prev_knot_ei = (int64_t)knot_ei;
        };
        auto get_err = [&](int64_t act_idx, int64_t pred_idx){
            return get_pred_diff(act_idx, pred_idx, begin);
        };

        uint64_t seg_count = make_segmentation_basic_pla(total_data_points, 
                    epsilon, in_fun, out_fun, get_err);

        brk_kval_vec.emplace_back(out[out.size()-1].get_last_x());
        insert_to_ind_vec(out[out.size()-1].get_last_x(), brk_kval_vec.size()-1);
        index_size = brk_kval_vec.size();
        uint64_t indir_max = 1ULL << lp_bits;
        for(int64_t iv = indirection_vec.size(); iv < indir_max; iv++){
            indirection_vec.emplace_back(brk_kval_vec.size()-1); 
        }

        auto [knot_si, knot_ei] = out[0].get_knot_intersection(out[0].get_last_x(), epsilon);
        brk_sa_indx_vec.emplace_back((uint64_t)knot_ei);
        brk_diff_vec.emplace_back(0);

        diff_packd.BuilidSignedPackedVector(brk_diff_vec);        
        y_range_beg.encode(brk_sa_indx_vec.data(), brk_sa_indx_vec.size());
        encode_knots(brk_kval_vec);
        // cout<<"#Segments: "<<index_size<<endl;
    }

    // suffix array version
    template <typename DataType>
    void Save(string indx_fn, DataType& data){
        essentials::save(f, index_fn.c_str());
        std::ofstream os(indx_fn, std::ios_base::app);
        /**
         * Set flags x 
         * whether indirection settings -> kmer size is stored [nth bit]
         * what kind of dictionary [n-1]
         * is fast rank enabled [n-2]
        */
        uint8_t flag_bits;
        if(data.is_kmer_size_needed()) flag_bits = 1;
        essentials_arank::save_pod(os, flag_bits);

        if(data.is_kmer_size_needed()){
            uint8_t kmer_size = data.get_kmer_size();
            essentials_arank::save_pod(os, kmer_size);
        }

        uint8_t t_lp_bits = lp_bits;
        essentials_arank::save_pod(os, t_lp_bits);
        vector<uint64_t> ind_diff_vec = get_processed_ind_vec();
        CBitPacking ind_diff_pack(false);
        ind_diff_pack.BuilidPackedVector(ind_diff_vec);
        ind_diff_pack.Save_os(os);
        uint32_t nBrkpnts = index_size;
        essentials_arank::save_pod(os, nBrkpnts);
        
        brk_uniform_diff_packed.Save_os(os);
        diff_packd.Save_os(os);

        y_range_beg.save(os);
        uint16_t t_err_th = epsilon;
        essentials_arank::save_pod(os, t_err_th);
        uint64_t uniq_kmers = err_dict_packed.GetNumElements();
        essentials_arank::save_pod(os, uniq_kmers);
        err_dict_packed.Save_os(os);

        // save_unpacked(indx_fn);
        // save_bit_info(indx_fn, ind_diff_pack);
    }

    // suffix array version
    template <typename DataType>
    void Load(string indx_fn, DataType &data){
        uint64_t total_bytes = essentials::load(f, indx_fn.c_str());
        std::ifstream is(indx_fn, std::ios::binary);
        std::streampos current_position = is.tellg();
        is.seekg(total_bytes, std::ios::beg);

        /**
         * Load flags x 
         * whether indirection settings -> kmer size is stored [nth bit]
        */
        uint8_t flag_bits;
        essentials_arank::load_pod(is, flag_bits);
        bool is_kmer_used = false;
        if(flag_bits & 1) is_kmer_used = true;
        indx_type = EXACT_PLA;

        if (is_kmer_used){
            uint8_t t_kmer_size;
            essentials_arank::load_pod(is, t_kmer_size);
            data.set_kmer_size(t_kmer_size);
        }
        uint64_t largest = data.get_largest();
        uint64_t total_bits;
        if(ceil(log2(largest)) == floor(log2(largest))) total_bits = ceil(log2(largest)) + 1;
        else total_bits = ceil(log2(largest));
        
        uint8_t lp_bits_t;
        // is.read(reinterpret_cast<char*>(&lp_bits_t), sizeof(lp_bits_t));
        essentials_arank::load_pod(is, lp_bits_t);
        lp_bits = lp_bits_t;
        shift_bits = total_bits - lp_bits;
        
        uint64_t indir_max = 1ULL << lp_bits;

        CBitPacking ind_diff_pack(false);
        ind_diff_pack.Load_is(is, (indir_max-1));
        // cout<<"ind diff vec isze: "<<indir_max-1<<endl;
        vector<int64_t> indir_diff_vec = ind_diff_pack.GetUnpackedVector();
        // indirection_table = new int64_t[indir_max];
        indirection_vec.resize(indir_max);
        indirection_vec[0] = 0;
        for(int64_t i=0; i<indir_max-1; i++){
            indirection_vec[i+1] = indirection_vec[i] + indir_diff_vec[i];
        }
        // how many breakpoints
        uint32_t index_size_t;
        essentials_arank::load_pod(is, index_size_t);
        // is.read(reinterpret_cast<char*>(&dic_size_t), sizeof(dic_size_t));
        
        index_size = index_size_t;
        // cout<<"#breakpoints: "<<dic_size<<endl;
        // brk_uniform_diff_packed.Load(fp, dic_size);
        brk_uniform_diff_packed.Load_is(is, index_size);
        max_data_ratio = ((1ULL << total_bits) - 1)/(index_size-1);
        diff_packd.Load_is(is, index_size-1);
        
        // cout<<"ef load\n";
        y_range_beg.load(is);
        // cout<<"max err load\n";
        uint16_t max_error_t;
        essentials_arank::load_pod(is, max_error_t);
        // is.read(reinterpret_cast<char*>(&max_error_t), sizeof(max_error_t));
        // err = fread(&max_error_t, sizeof(uint16_t), 1, fp);
        epsilon = (int64_t)max_error_t;

        uint64_t unique_kmers;
        essentials_arank::load_pod(is, unique_kmers);
        err_dict_packed.Load_is(is, unique_kmers);
    }

};