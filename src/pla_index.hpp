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
    bool is_fast_rank, is_rank_query;
    
    vector<uint64_t> indirection_vec, index_bitvec;

    CBitPacking brk_uniform_diff_packed, diff_packd;
    arank::ef_sequence<false> y_range_beg;
    sdsl::dac_vector_dp<> sa_diff_dac_vec;
    INDX_TYPE indx_type;

public:
    // constructor at index build
    pla_index(int64_t eps, int64_t lp_bits, uint64_t largest,
            bool is_fast_rank, uint64_t total_points, INDX_TYPE indx_type):
            epsilon(eps),
            lp_bits(lp_bits), is_fast_rank(is_fast_rank),
            indx_type(indx_type),
            total_data_points(total_points),
            brk_uniform_diff_packed(true), diff_packd(true)
    {
        uint64_t bits;
        if(ceil(log2(largest)) == floor(log2(largest))) bits = ceil(log2(largest)) + 1;
        else bits = ceil(log2(largest));
        if(bits < lp_bits){
            throw std::logic_error("lp cannot be greater than log2(largest)");
        }
        shift_bits = bits - lp_bits;
    }

    // constructor at index build
    pla_index(int64_t eps, vector<int64_t> &data,
            int64_t lp_bits = 16, 
            bool is_fast_rank = false, INDX_TYPE indx_type = REPEAT_PLA):
            epsilon(eps),
            lp_bits(lp_bits), is_fast_rank(is_fast_rank),
            indx_type(indx_type),
            brk_uniform_diff_packed(true), diff_packd(true)
    {
        total_data_points = data.size();
        int64_t largest = data[total_data_points - 1];
        uint64_t bits;
        if(ceil(log2(largest)) == floor(log2(largest))) bits = ceil(log2(largest)) + 1;
        else bits = ceil(log2(largest));
        if(bits < lp_bits){
            throw std::logic_error("lp cannot be greater than log2(largest)");
        }
        shift_bits = bits - lp_bits;
    }
    
    pla_index(bool is_rank_query = false, int64_t knot_bs_thres = 64):
    knot_bs_thres(knot_bs_thres), 
    is_rank_query(is_rank_query),
    brk_uniform_diff_packed(true), diff_packd(true)
    {}

    void insert_to_ind_vec(int64_t brk_beg_kval, int64_t brk_indx){
        int64_t lp_bit_value = brk_beg_kval >> shift_bits;
        if(indirection_vec.size() > lp_bit_value) return;
        for(int64_t iv=indirection_vec.size(); iv<lp_bit_value+1; iv++){
            indirection_vec.emplace_back(brk_indx);
        }
    }

    void encode_knots(const vector<int64_t> &brk_kval_vec){
        vector<int64_t> brk_uniform_diff_vec;
        uint64_t max_bits = shift_bits + lp_bits;
        
        max_data_ratio = ((1ULL << (max_bits)) - 1)/(dic_size - 1);
        
        for(int64_t i=0; i<dic_size; i++){
            int64_t diff = int64_t(brk_kval_vec[i]) - i*max_data_ratio;
            brk_uniform_diff_vec.emplace_back(diff);
        }
        brk_uniform_diff_packed.BuilidSignedPackedVector(brk_uniform_diff_vec);
    }
    
    void SetBit(uint64_t &number, uint64_t pos){
        uint64_t mask = 1ULL << (64-pos); // pos = [1, 64]
        number |= mask;
    }
    void save_unpacked(string indx_fn){
        indx_fn = indx_fn + ".unpacked";
        std::ofstream os(indx_fn, std::ios::binary);
        uint8_t t_lp_bits = lp_bits;
        essentials_arank::save_pod(os, t_lp_bits);
        os.write(reinterpret_cast<char const*>(indirection_vec.data()),
                    indirection_vec.size() * sizeof(uint64_t));
        uint32_t nBrkpnts = dic_size;
        essentials_arank::save_pod(os, nBrkpnts);
        
        if (indx_type == REPEAT_PLA)
            sa_diff_dac_vec.serialize(os);
        else if(indx_type == BASIC_PLA)
            diff_packd.Save_os(os);
        
        y_range_beg.save(os);
        uint16_t t_err_th = epsilon;
        essentials_arank::save_pod(os, t_err_th);
        essentials_arank::save_pod(os, is_fast_rank);
        if(is_fast_rank){
            os.write(reinterpret_cast<char const*>(index_bitvec.data()),
                    index_bitvec.size() * sizeof(uint64_t));
        }
    }

    void save_bit_info(string indx_fn, CBitPacking &ind_diff_pack){
        int pos = indx_fn.find(".bin");
        string pref = indx_fn.substr(0, pos);
        string suff = indx_fn.substr(pos);
        string bit_fn = pref+"_bitInfo"+suff;
        ofstream bit_info(bit_fn.c_str());
        bit_info<<"#Segments: "<<dic_size<<endl;
        bit_info<<"Indirection vector: each: "<<uint64_t(ind_diff_pack.GetPacketSize())<<endl;
        bit_info<<"knot kmers: each: "<<uint64_t(brk_uniform_diff_packed.GetPacketSize())<<endl;
        bit_info<<"Y beg each: "<<ceil(y_range_beg.num_bits()/dic_size)<<endl;
        if (indx_type == REPEAT_PLA)
            bit_info<<"Y end diff: each: "<<
                ceil((size_in_bytes(sa_diff_dac_vec)*8.0)/sa_diff_dac_vec.size())<<endl;
        else if(indx_type == BASIC_PLA)
            bit_info<<"Y end diff: each: "<<uint64_t(diff_packd.GetPacketSize())<<endl;        
        bit_info.close();
    }

    uint64_t get_size_in_bytes(){
        uint64_t total_size_in_bits = indirection_vec.size()*sizeof(uint64_t) + 
                uint64_t(brk_uniform_diff_packed.GetPacketSize())*brk_uniform_diff_packed.GetNumElements() +
                y_range_beg.num_bits();
        if(indx_type == BASIC_PLA){
            total_size_in_bits += uint64_t(diff_packd.GetPacketSize())*diff_packd.GetNumElements();
        }
        else if(indx_type == REPEAT_PLA){
            total_size_in_bits += size_in_bytes(sa_diff_dac_vec)*8.0;
        }
        if(is_fast_rank) total_size_in_bits += index_bitvec.size()*sizeof(uint64_t);
        return total_size_in_bits/8;
    }

    inline uint64_t get_num_segments() const {return dic_size;}

    const vector<uint64_t> get_processed_ind_vec() const{
        vector<uint64_t> ind_diff_vec;
        int64_t _count = 0;
        for(int64_t i=1; i<indirection_vec.size(); i++){
            uint64_t diff = indirection_vec[i] - indirection_vec[i-1];
            if(diff == 0) _count++;
            ind_diff_vec.push_back(diff);
        }
        return ind_diff_vec;
    }

    int64_t binary_search(int64_t lo, int64_t hi, const int64_t kval) const{
        int64_t mid, midVal;
        while(1){
            if(hi-lo <= knot_bs_thres){
                for (int64_t i = hi; i >= lo; i--)
                {
                    if(brk_uniform_diff_packed.GetValueAt(i)+max_data_ratio*i <= kval) return i;
                }
            } 
            mid = (lo + hi) >> 1;
            midVal = brk_uniform_diff_packed.GetValueAt(mid)+max_data_ratio*mid;
            if (midVal == kval) return mid;
            else if(midVal > kval){ 
                hi = mid;
            }
            else{ 
                lo = mid;
            }
        }  
    }
    /*
    template <typename DataType>
    uint64_t get_first_index(int64_t idx, const uint64_t query_val, 
            const DataType &data) const{
        if(!is_fast_rank){
            while(1){
                if(--idx < 0) return data[0];
                if(data[idx] != query_val) return data.get_sa_val(idx+1);
            }
        }
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
    }*/

    // suffix array version
    template <typename DataType>
    void Save(string indx_fn, DataType& data){
        std::ofstream os(indx_fn, std::ios::binary);
        /**
         * Set flags x 
         * whether indirection settings -> kmer size is stored [nth bit]
         * what kind of dictionary [n-1]
         * is fast rank enabled [n-2]
        */
        uint8_t flag_bits = 0;
        // if(data.is_kmer_size_needed()) flag_bits = 1;
        if(indx_type == BASIC_PLA) flag_bits |= 2;
        if(is_fast_rank) flag_bits |= 4;
        essentials_arank::save_pod(os, flag_bits);

        // if(data.is_kmer_size_needed()){
        //     uint8_t kmer_size = data.get_kmer_size();
        //     essentials_arank::save_pod(os, kmer_size);
        // }

        uint8_t t_lp_bits = lp_bits;
        essentials_arank::save_pod(os, t_lp_bits);
        vector<uint64_t> ind_diff_vec = get_processed_ind_vec();
        CBitPacking ind_diff_pack(false);
        ind_diff_pack.BuilidPackedVector(ind_diff_vec);
        ind_diff_pack.Save_os(os);
        uint32_t nBrkpnts = dic_size;
        
        essentials_arank::save_pod(os, nBrkpnts);
        // os.write(reinterpret_cast<char const*>(&nBrkpnts),
                    //   sizeof(nBrkpnts));
        
        brk_uniform_diff_packed.Save_os(os);

        if (indx_type == REPEAT_PLA)
            sa_diff_dac_vec.serialize(os);
        else if(indx_type == BASIC_PLA)
            diff_packd.Save_os(os);

        y_range_beg.save(os);
        uint16_t t_err_th = epsilon;
        essentials_arank::save_pod(os, t_err_th);
        // os.write(reinterpret_cast<char const*>(&t_err_th),
        //               sizeof(t_err_th));
        essentials_arank::save_pod(os, is_fast_rank);
        if(is_fast_rank){
            os.write(reinterpret_cast<char const*>(index_bitvec.data()),
                    index_bitvec.size() * sizeof(uint64_t));
        }

        // save_unpacked(indx_fn);
        // save_bit_info(indx_fn, ind_diff_pack);
    }
    
    // suffix array version
    template <typename DataType>
    void Load(string indx_fn, DataType &data){
        /**
         * Load flags x 
         * whether indirection settings -> kmer size is stored [nth bit]
         * what kind of dictionary [n-1]
         * is fast rank enabled [n-2]
        */
        std::ifstream is(indx_fn, std::ios::binary);
        uint8_t flag_bits;
        essentials_arank::load_pod(is, flag_bits);
        bool is_kmer_used = false;
        if(flag_bits & 1) is_kmer_used = true;
        if(flag_bits & 2) indx_type = BASIC_PLA;
        else indx_type = REPEAT_PLA;
        if(flag_bits & 4) is_fast_rank = true;
        else is_fast_rank = false;

        // if (is_kmer_used){
        //     uint8_t t_kmer_size;
        //     essentials_arank::load_pod(is, t_kmer_size);
        //     data.set_kmer_size(t_kmer_size);
        // }
        total_data_points = data.size();
        uint64_t largest = data[largest - 1];
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
        uint32_t dic_size_t;
        essentials_arank::load_pod(is, dic_size_t);
        // is.read(reinterpret_cast<char*>(&dic_size_t), sizeof(dic_size_t));
        
        dic_size = dic_size_t;
        // cout<<"#breakpoints: "<<dic_size<<endl;
        // brk_uniform_diff_packed.Load(fp, dic_size);
        brk_uniform_diff_packed.Load_is(is, dic_size);
        max_data_ratio = ((1ULL << total_bits) - 1)/(dic_size-1);
        if (indx_type == REPEAT_PLA)
            sa_diff_dac_vec.load(is);
        else if(indx_type == BASIC_PLA)
            diff_packd.Load_is(is, dic_size-1);
        
        // cout<<"ef load\n";
        y_range_beg.load(is);
        // cout<<"max err load\n";
        uint16_t max_error_t;
        essentials_arank::load_pod(is, max_error_t);
        // is.read(reinterpret_cast<char*>(&max_error_t), sizeof(max_error_t));
        // err = fread(&max_error_t, sizeof(uint16_t), 1, fp);
        epsilon = (int64_t)max_error_t;
        essentials_arank::load_pod(is, is_fast_rank);
        if(is_fast_rank){
            int64_t bitvec_size = total_data_points/64 + 1;
            index_bitvec.resize(bitvec_size);
            is.read(reinterpret_cast<char*>(index_bitvec.data()), 
                static_cast<std::streamsize>(sizeof(uint64_t) * bitvec_size));  
            // err = fread(&index_bitvec[0], sizeof(uint64_t), bitvec_size, fp);
        }
    }

    // suffix array version
    template <typename DataType>
    int64_t query(const int64_t query_val, const DataType &data) const{
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
        
        // TODO: Change to rank and search features?
        auto lo = data.begin() + max(query_pred - epsilon, int64_t(0));
        auto hi = data.begin() + min(query_pred + epsilon, int64_t(total_data_points)-1);
        return *lower_bound(lo, hi, query_val);
        // if (!is_rank_query){
            
        // }
        // return get_first_index(data.BinarySearch(query , max(query_pred - epsilon, int64_t(0)), 
        //         min(query_pred + epsilon, int64_t(total_data_points)-1), is_rank_query), query_val, data);

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
        if(is_fast_rank) calc_index_bv(begin);
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
        if(is_fast_rank) calc_index_bv(begin);
    }
};