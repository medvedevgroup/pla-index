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


using namespace std;

enum INDX_TYPE {BASIC_PLA, REPEAT_PLA};

class pla_index{
private:
    using canonical_segment = typename OptimalPiecewiseLinearModel<int64_t, uint64_t>::CanonicalSegment;
    int64_t epsilon, max_data_ratio, lp_bits, knot_bs_thres;
    uint64_t shift_bits, dic_size, total_data_points;
    bool isFirstIndexReturned;
    
    vector<uint64_t> indirection_vec, index_bitvec;

    CBitPacking brk_uniform_diff_packed, diff_packd;
    arank::ef_sequence<false> y_range_beg;
    sdsl::dac_vector_dp<> sa_diff_dac_vec;
    INDX_TYPE indx_type;

public:
    pla_index(int64_t eps, int64_t lp_bits, uint64_t largest,
            bool isFirstIndxRet, uint64_t total_points, INDX_TYPE indx_type);
    
    pla_index(int64_t knot_bs_thres, uint64_t largest, 
        uint64_t total_points, string dic_fn, INDX_TYPE indx_type);

    void insert_to_ind_vec(int64_t brk_beg_kval, int64_t brk_indx);

    void encode_knots(const vector<int64_t> &brk_kval_vec);    
    
    void SetBit(uint64_t &number, uint64_t pos);
    void Save(string dic_fn);
    void save_unpacked(string dic_fn);
    void save_bit_info(string dic_fn, CBitPacking &ind_diff_pack);
    uint64_t get_size_in_bytes();
    inline uint64_t get_num_segments() const {return dic_size;}

    const vector<uint64_t> get_processed_ind_vec() const;
    void Load(string dic_fn, int64_t total_bits);
    int64_t binary_search(int64_t lo, int64_t hi, const int64_t kval) const;
    
    template <typename DataType>
    uint64_t get_first_index(uint64_t idx, const DataType &data) const{
        uint64_t bv_idx = idx >> 6; // x/64
        uint64_t bv_pos = idx - (bv_idx<<6) + 1; // x%64 +1
        uint64_t mask = 1ULL << (64-bv_pos); // pos = [1, 64]
        uint64_t offset = 0;
        while (1)
        {
            if(bv_pos == 0){
                bv_pos = 64;
                mask = 1ULL;
                bv_idx--;
            }
            if(index_bitvec[bv_idx] & mask){
                return data.get_sa_val(idx-offset);
            }
            bv_pos--;
            offset++;
            mask <<= 1;
        }        
    }
    
    

    // suffix array version
    template <typename DataType>
    int64_t query(const int64_t query_val, const string &query, const DataType &data) const{
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
                int64_t end_idx = (lp_bit_value == indirection_vec.size() - 1) ? dic_size - 1 : indirection_vec[lp_bit_value + 1] - 1;
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
                int64_t brk_end_sa_indx;
                switch (indx_type){
                case REPEAT_PLA:
                    {
                    int64_t diff_val = sa_diff_dac_vec[brkpnt];
                    brk_end_sa_indx = (diff_val & 1)
                                        ? y_range_beg.access(brkpnt + 1) - (diff_val >> 1) - epsilon
                                        : y_range_beg.access(brkpnt + 1) + (diff_val >> 1) - epsilon;
                    break;
                    }
                case BASIC_PLA:
                    {
                    brk_end_sa_indx = y_range_beg.access(brkpnt + 1) - 
                        diff_packd.GetValueAt(brkpnt) - epsilon;
                    break;
                    }
                }
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
        query_pred = std::clamp(query_pred, int64_t(0), int64_t(total_data_points) - 1);
        
        if (!isFirstIndexReturned){
            return data.BinarySearch(query , max(query_pred - epsilon, int64_t(0)), 
                min(query_pred + epsilon, int64_t(total_data_points)-1), isFirstIndexReturned);
        }
        return get_first_index(data.BinarySearch(query , max(query_pred - epsilon, int64_t(0)), 
                min(query_pred + epsilon, int64_t(total_data_points)-1), isFirstIndexReturned), data);

    }

    template<typename RandomIt>
    void calc_index_bv(RandomIt data){
        index_bitvec.resize(total_data_points/64 + 1);
        int64_t curr_kval, prev_kval = data[0];
        uint64_t position = 1, indx = 0;
        if(prev_kval != -1){
            SetBit(index_bitvec[indx], position);
        }
        position++;
        for(int64_t i=1; i< total_data_points; i++){
            if(position == 65){
                position = 1;
                indx++;
            }
            curr_kval = data[i];
            // -1 check for suffix array type applicaiton
            if(curr_kval != -1 && curr_kval != prev_kval){            
                SetBit(index_bitvec[indx], position);            
            }
            position++;
            prev_kval = curr_kval;
        }
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

    template<typename RandomIt>
    void build_index(const RandomIt &begin){
        if(indx_type == BASIC_PLA)
            build_basic_pla_index(begin);
        else
            build_rep_stretch_pla_index(begin);
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
        dic_size = brk_kval_vec.size();
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
        if(isFirstIndexReturned) calc_index_bv(begin);
        cout<<"#breakpoints: "<<dic_size<<endl;
    }

    template<typename RandomIt>
    void build_rep_stretch_pla_index(const RandomIt &begin){
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
            diff_from_curr_start = prev_knot_ei - (int64_t)knot_si;
            if(diff_from_curr_start < 0){
                diff_from_curr_start *= -2;
                diff_from_curr_start++;
            }
            else{
                diff_from_curr_start *= 2;
            }
            brk_diff_vec.emplace_back(diff_from_curr_start);
            prev_knot_ei = (int64_t)knot_ei;
        };
        auto get_err = [&](int64_t act_idx, int64_t pred_idx){
            return get_pred_diff(act_idx, pred_idx, begin);
        };

        uint64_t seg_count = make_segmentation_rep_pla(total_data_points, epsilon, in_fun, out_fun, get_err);

        brk_kval_vec.emplace_back(out[out.size()-1].get_last_x());
        insert_to_ind_vec(out[out.size()-1].get_last_x(), brk_kval_vec.size()-1);
        dic_size = brk_kval_vec.size();
        uint64_t indir_max = 1ULL << lp_bits;
        for(int64_t iv = indirection_vec.size(); iv < indir_max; iv++){
            indirection_vec.emplace_back(brk_kval_vec.size()-1); 
        }

        auto [knot_si, knot_ei] = out[0].get_knot_intersection(out[0].get_last_x(), epsilon);
        brk_sa_indx_vec.emplace_back((uint64_t)knot_ei);
        brk_diff_vec.emplace_back(0);

        sdsl::dac_vector_dp<> temp(brk_diff_vec);
        sa_diff_dac_vec = temp;
        y_range_beg.encode(brk_sa_indx_vec.data(), brk_sa_indx_vec.size());
        encode_knots(brk_kval_vec);
        if(isFirstIndexReturned) calc_index_bv(begin);
        cout<<"#breakpoints: "<<dic_size<<endl;
    }
};