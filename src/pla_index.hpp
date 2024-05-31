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

using namespace std;

enum INDX_TYPE {BASIC_PLA, REPEAT_PLA};

#define DEBUG_PRINT 0
#define KNOT_BS_THRES 64

class pla_index{
private:
    using canonical_segment = typename OptimalPiecewiseLinearModel<int64_t, uint64_t>::CanonicalSegment;
    int64_t epsilon, max_data_ratio, lp_bits;
    uint64_t shift_bits, index_size, total_data_points, 
        x_width, ef_width, y_diff_width, pref_width, x_const_width;

    bool is_fast_rank, is_rank_query;
    
    vector<uint64_t> index_bitvec, x_enc_vec, prefix_vec;

    arank::ef_sequence<false> y_range_beg;
    INDX_TYPE indx_type;
    struct brkpnt_values{
        int64_t x1, diff1, ef1, ef2, 
            x2_x1, brkpnt;
    };
    
public:
    // constructor at index build
    pla_index(int64_t eps, bool is_fast_rank, 
            uint64_t total_points, INDX_TYPE it):
            epsilon(eps), is_fast_rank(is_fast_rank),
            total_data_points(total_points),
            indx_type(it)
    {}
    
    // constructor at query
    pla_index(bool is_rank_query):
        is_rank_query(is_rank_query)
    {}


    void SetBit(uint64_t &number, uint64_t pos){
        uint64_t mask = 1ULL << (64-pos); // pos = [1, 64]
        number |= mask;
    }

    void enc_value(uint64_t v, vector<uint64_t> &m_bits, const uint64_t& mask, uint64_t& width,
            uint64_t &m_cur_block, uint64_t &m_cur_shift){
        m_bits[m_cur_block] &= ~(mask << m_cur_shift);
        m_bits[m_cur_block] |= v << m_cur_shift;

        uint64_t res_shift = 64 - m_cur_shift;
        if (res_shift < width) {
            ++m_cur_block;
            m_bits.push_back(0);
            m_bits[m_cur_block] &= ~(mask >> res_shift);
            m_bits[m_cur_block] |= v >> res_shift;
            m_cur_shift = -res_shift;
        }

        m_cur_shift += width;

        if (m_cur_shift == 64) {
            m_cur_shift = 0;
            ++m_cur_block;
            m_bits.push_back(0);
        }
    }

    void encode_values(const vector<uint64_t> &brk_kval_vec, 
        const vector<uint64_t> &brk_diff_vec, int64_t lookup_count){
        
        // lp bits taken input as lookup count
        lp_bits = ceil(log2((brk_kval_vec.size()/lookup_count)));
        uint64_t max_bits = shift_bits;
        shift_bits -= lp_bits;
        max_data_ratio = ((1ULL << (max_bits)) - 1)/(index_size - 1);
        if(DEBUG_PRINT){
            cout<<"lp bits: "<<lp_bits
                <<" shift bits: "<<shift_bits
                <<" total bits: "<<max_bits
                <<endl;
        }
        uint64_t largest = 0;

        for(int64_t i=0; i<index_size; i++){
            int64_t diff = int64_t(brk_kval_vec[i]) - i*max_data_ratio;
            uint64_t val = add_bit_based_on_sign(diff);
            if(val > largest) largest = val;
        }
        
        auto bits_needed = [](int64_t n) -> int64_t {
            return (ceil(log2(n)) == floor(log2(n)))
                ? ceil(log2(n)) + 1
                : ceil(log2(n));
        };

        
        // +1 because includes negative number
        y_diff_width = bits_needed(epsilon*2) + 1;
        x_width = bits_needed(largest);
        pref_width = bits_needed(index_size - 1);
        ef_width = y_range_beg.cv_width();

        auto get_mask = [] (int64_t n) -> uint64_t {
            return (1ULL << n) - 1;
        };

        const uint64_t y_diff_mask = get_mask(y_diff_width);
        const uint64_t x_mask = get_mask(x_width);
        const uint64_t pref_mask = get_mask(pref_width);
        const uint64_t ef_mask = get_mask(ef_width);

        x_const_width = x_width + ef_width + y_diff_width;

        if(DEBUG_PRINT){
            cout<<"width\n"
                <<"x: "<<x_width
                <<" y_diff: "<<y_diff_width
                <<" pref: "<<pref_width
                <<" ef: "<<ef_width
                <<" x_const_width: "<<x_const_width
                <<" max data ratio: "<<max_data_ratio
                <<" largest: "<<largest
                <<endl;
        }

        

        uint64_t prefix_count=0, x_enc_block = 0, x_enc_shift = 0, 
            prefix_enc_block=0, prefix_enc_shift=0;
        
        x_enc_vec.push_back(0);
        prefix_vec.push_back(0);

        for(int64_t i=0; i<index_size; i++){
            int64_t diff = int64_t(brk_kval_vec[i]) - i*max_data_ratio;
            uint64_t val = add_bit_based_on_sign(diff);

            enc_value(val, x_enc_vec, x_mask, 
                x_width, x_enc_block, x_enc_shift);

            enc_value(y_range_beg.cv_access(i), x_enc_vec, 
                ef_mask, ef_width, x_enc_block, x_enc_shift);
            
            enc_value(brk_diff_vec[i], x_enc_vec, y_diff_mask, 
                y_diff_width, x_enc_block, x_enc_shift);

            if(DEBUG_PRINT){
                if(i > 1635555 && i < 1635558){
                    cout<<i<<" x_val: "<<brk_kval_vec[i]
                        <<" cv val: "<<y_range_beg.cv_access(i)
                        <<" y_start[]: "<<y_range_beg.access(i, 
                                    y_range_beg.cv_access(i), ef_width)
                        <<" diff: "<<brk_diff_vec[i]
                        <<endl;
                }
            }
            
            int64_t lp_bit_value = brk_kval_vec[i] >> shift_bits;
            if(prefix_count <= lp_bit_value){
                for(int64_t pl=prefix_count; pl<lp_bit_value+1; pl++){
                    prefix_count++;
                    enc_value(i, prefix_vec, pref_mask, pref_width, 
                        prefix_enc_block, prefix_enc_shift);
                }
            }
        }

        // this ensures upto index_size-1 will be checked for the last entry
        enc_value(index_size-1, prefix_vec, pref_mask, pref_width, 
            prefix_enc_block, prefix_enc_shift);
    }
    
    
    // uint64_t get_size_in_bytes();
    inline uint64_t get_num_segments() const {return index_size;}
    
    template <typename DataType>
    uint64_t get_first_index(int64_t idx, const uint64_t query_val, 
            const DataType &data) const{
        if(!is_fast_rank){
            while(1){
                if(--idx < 0) return data.get_sa_val(0);
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
    }

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
        if(data.is_kmer_size_needed()) flag_bits = 1;
        if(indx_type == BASIC_PLA) flag_bits |= 2;
        if(is_fast_rank) flag_bits |= 4;
        essentials_arank::save_pod(os, flag_bits);

        if(data.is_kmer_size_needed()){
            uint8_t kmer_size = data.get_kmer_size();
            essentials_arank::save_pod(os, kmer_size);
        }

        uint8_t t_lp_bits = lp_bits;
        essentials_arank::save_pod(os, t_lp_bits);

        uint8_t pref_width_t = pref_width;
        essentials_arank::save_pod(os, pref_width_t);
        essentials_arank::save_vec(os, prefix_vec);

        uint32_t index_size_t = index_size;        
        essentials_arank::save_pod(os, index_size_t);
        
        uint8_t x_width_t = x_width;
        essentials_arank::save_pod(os, x_width_t);
        uint8_t ef_width_t = ef_width;
        essentials_arank::save_pod(os, ef_width_t);

        essentials_arank::save_vec(os, x_enc_vec);

        y_range_beg.save(os);

        uint16_t epsilon_t = epsilon;
        essentials_arank::save_pod(os, epsilon_t);

        if(is_fast_rank){
            os.write(reinterpret_cast<char const*>(index_bitvec.data()),
                    index_bitvec.size() * sizeof(uint64_t));
        }
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

        if (is_kmer_used){
            uint8_t t_kmer_size;
            essentials_arank::load_pod(is, t_kmer_size);
            data.set_kmer_size(t_kmer_size);
        }

        total_data_points = data.size();
        // uint64_t largest = data.get_largest();
        uint64_t total_bits = 2*data.get_kmer_size();
        // if(ceil(log2(largest)) == floor(log2(largest))) total_bits = ceil(log2(largest)) + 1;
        // else total_bits = ceil(log2(largest));
        
        uint8_t lp_bits_t;
        essentials_arank::load_pod(is, lp_bits_t);
        lp_bits = lp_bits_t;
        shift_bits = total_bits - lp_bits;

        uint8_t pref_width_t;
        essentials_arank::load_pod(is, pref_width_t);
        pref_width = pref_width_t;
        essentials_arank::load_vec(is, prefix_vec);

        // if(DEBUG_PRINT){
            // cout<<"Size: Prefix table: \n"
            // <<"width: "<<pref_width
            // <<" size: "<<prefix_vec.size()
            // <<endl;
        // }
        
        uint32_t index_size_t;
        essentials_arank::load_pod(is, index_size_t);
        index_size = index_size_t;
        max_data_ratio = ((1ULL << total_bits) - 1)/(index_size-1);

        // if(DEBUG_PRINT)
            cout<<"#breakpoints: "<<index_size<<endl;
            cout<<"lp bits: "<<lp_bits<<endl;
        
        uint8_t x_width_t, ef_width_t;
        essentials_arank::load_pod(is, x_width_t);
        essentials_arank::load_pod(is, ef_width_t);
        x_width = x_width_t;
        ef_width = ef_width_t;

        essentials_arank::load_vec(is, x_enc_vec);

        y_range_beg.load(is);
        
        uint16_t epsilon_t;
        essentials_arank::load_pod(is, epsilon_t);
        epsilon = epsilon_t;
        
        y_diff_width = (ceil(log2(epsilon*2)) == floor(log2(epsilon*2)))
                ? ceil(log2(epsilon*2)) + 1 + 1
                : ceil(log2(epsilon*2)) + 1;
        
        x_const_width = x_width + ef_width + y_diff_width;

        // if(DEBUG_PRINT){
            // cout<<"width\n"
            //     <<"x: "<<x_width
            //     <<" y_diff: "<<y_diff_width
            //     <<" pref: "<<pref_width
            //     <<" ef: "<<ef_width
            //     <<" x_const_width: "<<x_const_width
            //     <<endl;
        // }
        
        cout<<"PLA index loaded"<<endl;
        
        if(is_fast_rank){
            int64_t bitvec_size = total_data_points/64 + 1;
            index_bitvec.resize(bitvec_size);
            is.read(reinterpret_cast<char*>(index_bitvec.data()), 
                static_cast<std::streamsize>(sizeof(uint64_t) * bitvec_size));
        }
    }

    inline int64_t get_signed_value(int64_t value) const{
        return value & 1
            ? -(value>>1)
            : value>>1;
    }


    inline int64_t predict(int64_t q, int64_t x1, int64_t y1, int64_t x2_x1, int64_t y2) const{
        return round(y1 + ((double)(q - x1)/x2_x1 )*(y2 - y1));
    }

    inline uint64_t get_val(uint64_t& pos, const uint64_t& mask, const uint64_t& width) const{
        uint64_t block = pos >> 6;
        uint64_t shift = pos & 63;
        pos += width;
        return shift + width <= 64
                    ? x_enc_vec[block] >> shift & mask
                    : (x_enc_vec[block] >> shift) | (x_enc_vec[block + 1] << (64 - shift) & mask);
    }

    void binary_search(uint64_t lo, uint64_t hi, const int64_t query_val, 
            brkpnt_values &bv) const
    {
        int64_t mid, midVal, pos, rest_val,
            block, shift, prev_x = -1, curr_x, rest_bits = x_const_width - x_width, mdr;
        const uint64_t x_total_width = x_width + rest_bits;
        const uint64_t rest_mask = (1ULL << rest_bits) - 1;
        const uint64_t ef_mask = (1ULL << ef_width) - 1;
        const uint64_t x_mask = (1ULL << x_width) - 1;

        while(1){
            if(hi-lo <= KNOT_BS_THRES){
                mdr = max_data_ratio * lo;
                pos = lo*x_total_width;
                for (int64_t i = lo; i <= hi; i++)
                {
                    block = pos >> 6;
                    shift = pos & 63;
                    curr_x = shift + x_width <= 64
                        ? x_enc_vec[block] >> shift & x_mask
                        : (x_enc_vec[block] >> shift) | (x_enc_vec[block + 1] << (64 - shift) & x_mask);
                    curr_x = get_signed_value(curr_x) + mdr;
                    // if(DEBUG_PRINT){
                    //     cout<<i<<" "<<lo<<" "<<hi<<" "<<x_total_width<<" "<<pos<<" "
                    //     <<curr_x<<" "<<query_val<<endl;
                    // }
                    
                    if(curr_x >= query_val) {
                        pos += x_width;
                        block = pos >> 6;
                        shift = pos & 63;
                        rest_val = shift + rest_bits <= 64
                            ? x_enc_vec[block] >> shift & rest_mask
                            : (x_enc_vec[block] >> shift) | 
                                (x_enc_vec[block + 1] << (64 - shift) & rest_mask);
                        bv.ef2 = rest_val & ef_mask;

                        if(prev_x != -1){
                            bv.x1 = prev_x;
                            pos -= x_total_width;
                            block = pos >> 6;
                            shift = pos & 63;
                            rest_val = shift + rest_bits <= 64
                                ? x_enc_vec[block] >> shift & rest_mask
                                : (x_enc_vec[block] >> shift) | 
                                    (x_enc_vec[block + 1] << (64 - shift) & rest_mask);                        
                            
                            bv.ef1 = rest_val & ef_mask;
                            bv.diff1 = rest_val >> ef_width;
                            // bv.x2_x1 = curr_x - prev_x;
                            // if(DEBUG_PRINT){
                            //     cout<<"In bs x:\n";
                            //     cout<<"diff1: "<<bv.diff1<<" ef1 "<<bv.ef1
                            //         <<" ef2: "<<bv.ef2<<endl;
                            // }                            
                        }

                        // if prev_x == -1, this means bv already has all the value it needs
                        bv.x2_x1 = curr_x - bv.x1;
                        bv.brkpnt = i-1;
                        return;
                    }
                    prev_x = curr_x;
                    pos += x_total_width;
                    mdr += max_data_ratio;
                }
                cout<<"Not found\n";
                cout<<query_val<<endl;
                exit(-1);
            } 
            
            mid = (lo + hi) >> 1;
            // prefetch

            pos = mid * x_total_width;
            mdr = max_data_ratio*mid;
            
            block = pos >> 6;
            shift = pos & 63;
            midVal = shift + x_width <= 64
                ? x_enc_vec[block] >> shift & x_mask
                : (x_enc_vec[block] >> shift) | (x_enc_vec[block + 1] << (64 - shift) & x_mask);
            midVal = get_signed_value(midVal) + mdr;

            // if(DEBUG_PRINT)
            //     cout<<mid<<" "<<midVal<<" "<<query_val<<endl;
            
            if (midVal == query_val) {
                pos += x_width;
                block = pos >> 6;
                shift = pos & 63;
                rest_val = shift + rest_bits <= 64
                    ? x_enc_vec[block] >> shift & rest_mask
                    : (x_enc_vec[block] >> shift) | 
                        (x_enc_vec[block + 1] << (64 - shift) & rest_mask);
                
                bv.x1 = midVal;
                bv.ef1 = rest_val & ef_mask;
                bv.diff1 = rest_val >> ef_width;
                bv.brkpnt = mid;
                
                pos += rest_bits;
                block = pos >> 6;
                shift = pos & 63;
                bv.x2_x1 = shift + x_width <= 64
                    ? x_enc_vec[block] >> shift & x_mask
                    : (x_enc_vec[block] >> shift) | (x_enc_vec[block + 1] << (64 - shift) & x_mask);
                bv.x2_x1  = get_signed_value(bv.x2_x1) + mdr+ max_data_ratio - midVal;
                
                pos += x_width;
                block = pos >> 6;
                shift = pos & 63;
                rest_val = shift + rest_bits <= 64
                    ? x_enc_vec[block] >> shift & rest_mask
                    : (x_enc_vec[block] >> shift) | 
                        (x_enc_vec[block + 1] << (64 - shift) & rest_mask);
                
                bv.ef2 = rest_val & ef_mask;
                // if(DEBUG_PRINT)
                //     cout<<"ef "<<bv.ef1<<" "<<rest_val<<" "<<bv.diff1<<" "<<bv.x2_x1<<endl;
                return;
            }
            else if(midVal > query_val){ 
                hi = mid;
            }
            else{ 
                lo = mid;
            }
        }

    }

    template <typename DataType>
    inline int64_t query(const int64_t query_val, 
        const string &query, const DataType &data) const
    {
        if(indx_type == BASIC_PLA){
            return query_concise_basic_pla(query_val, query, data);
        }
        else{
            ;
        }
    }

    template <typename DataType>
    int64_t query_concise_basic_pla(const int64_t query_val, 
        const string &query, const DataType &data) const
    {
        int64_t query_pred, lo, brkpnt;
        brkpnt_values bv;
        uint64_t lp_bit_value = query_val >> shift_bits;
        const uint64_t ef_mask = (1ULL << ef_width) - 1;
        const uint64_t x_mask = (1ULL << x_width) - 1;
        const uint64_t y_diff_mask = (1ULL << y_diff_width) - 1;
        const uint64_t pref_mask = (1ULL << pref_width) - 1;

        uint64_t pos = lp_bit_value*pref_width;
        uint64_t block = pos >> 6;
        uint64_t shift = pos & 63;
        uint64_t indx1 = (shift + pref_width <= 64) 
            ?   prefix_vec[block]>>shift & pref_mask
            :   (prefix_vec[block]>>shift) | (prefix_vec[block+1] << (64-shift) & pref_mask);
        
        pos += pref_width;
        block = pos >> 6;
        shift = pos & 63;
        uint64_t indx2 = (shift + pref_width <= 64) 
            ?   prefix_vec[block]>>shift & pref_mask
            :   (prefix_vec[block]>>shift) | (prefix_vec[block+1] << (64-shift) & pref_mask);
        
        int64_t x2, ef2, diff2, mdr = max_data_ratio*indx1;
        uint64_t rest_bits = x_const_width - x_width;
        const uint64_t rest_mask = (1ULL << rest_bits)-1;
        uint64_t rest_val;

        if(indx1 > 0){
            pos = (indx1-1) * x_const_width;
            x2 = get_signed_value(get_val(pos, x_mask, x_width)) + mdr - max_data_ratio;
            rest_val = get_val(pos, rest_mask, rest_bits);

            ef2 = rest_val & ef_mask;
            diff2 = rest_val >> ef_width;
        }
        else{
            pos = indx1 * x_const_width;
        }

        bv.x1 = get_signed_value(get_val(pos, x_mask, x_width)) + mdr;
        rest_val = get_val(pos, rest_mask, rest_bits);

        bv.ef1 = rest_val & ef_mask;
        bv.diff1 = rest_val >> ef_width;

        // if(DEBUG_PRINT){
        //     cout<<"Before match\n"
        //         <<"prefix index 1: "<<indx1<<" index 2: "<<indx2<<endl
        //         <<"bv x1: "<<bv.x1<<" ef1: "<<bv.ef1<<" diff1: "<<bv.diff1<<endl;
        //     if(indx1 > 0){
        //         cout<<"x2: "<<x2<<" ef2: "<<ef2<<" diff2: "<<diff2<<endl;
        //     }
        // }

        if(bv.x1 == query_val){
            query_pred = y_range_beg.access(indx1, bv.ef1, ef_width) - epsilon;
        }
        else{
            if(bv.x1 < query_val){
                // indx2 - 1 should also work fine
                binary_search(indx1+1, indx2, query_val, bv);
            }
            else{
                bv.x2_x1 = bv.x1 - x2;
                bv.x1 = x2; // x2 is the previous element
                bv.ef2 = bv.ef1;
                bv.ef1 = ef2;
                bv.diff1 = diff2;
                bv.brkpnt = indx1 - 1;
            }

            auto [knot_si, val2] = y_range_beg.pair(bv.brkpnt, bv.ef1, bv.ef2, ef_width);
            int64_t knot_ei = val2 -  get_signed_value(bv.diff1);

            query_pred = predict(query_val, bv.x1, knot_si, bv.x2_x1,
                knot_ei) - epsilon;

            // if(DEBUG_PRINT){
            //     cout<<"After query_pred\n"
            //         <<"x1: "<<bv.x1<<" q: "<<query_val<<" x2: "<<bv.x1 + bv.x2_x1
            //         <<"\ny1: "<<knot_si<<" y2: "<<knot_ei
            //         <<" next knot si: "<<val2
            //         <<" pred: "<<query_pred
            //         <<endl;
                    
            // }
        }
        
        query_pred = std::clamp(query_pred, int64_t(0), int64_t(total_data_points) - 1);
        
        return !is_rank_query
            ? data.BinarySearch(query , max(query_pred - epsilon, int64_t(0)), 
                min(query_pred + epsilon, int64_t(total_data_points)-1), false)
            : get_first_index(data.BinarySearch(query , max(query_pred - epsilon, int64_t(0)), 
                min(query_pred + epsilon, int64_t(total_data_points)-1), true), query_val, data);        
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
        RandomIt& data) const{
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
    void build_index(const RandomIt &begin, const uint64_t lookup_count,
        int64_t kmer_size){
        if(indx_type == BASIC_PLA)
            build_basic_pla_index(begin, lookup_count, kmer_size);
        else
            // build_rep_stretch_pla_index(begin, lookup_count);
            ;
    }

    /**
     * @brief add a 1 bit at the end of a negative number, and 
     * 0 at the end of a positive number
    */
    inline uint64_t add_bit_based_on_sign(int64_t &num){
        return num < 0
            ? ((-num)<<1) | 1
            : num << 1;
    }

    template<typename RandomIt>
    void build_basic_pla_index(const RandomIt &begin, 
        const uint64_t lookup_count,
        int64_t kmer_size){
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
            if(DEBUG_PRINT){
                if(cs.get_first_x() == 4398046511102){
                    cout<<"Inside out func\n"
                        <<"knot start: "<<cs.get_first_x()
                        <<" knot end: "<<last_x
                        <<" knot si: "<<knot_si
                        <<" knot ei: "<<knot_ei
                        <<endl;
                }
            }
            
        };
        auto get_err = [&](int64_t act_idx, int64_t pred_idx){
            return get_pred_diff(act_idx, pred_idx, begin);
        };

        uint64_t seg_count = make_segmentation_basic_pla(total_data_points, 
                    epsilon, in_fun, out_fun, get_err);

        if(brk_kval_vec[brk_kval_vec.size()-1] != out[0].get_last_x()){
            brk_kval_vec.emplace_back(out[0].get_last_x());
            uint64_t knot_si = begin.rank_of_last_entry() + epsilon;
            brk_sa_indx_vec.emplace_back(knot_si);
            diff_from_curr_start = (int64_t)knot_si - prev_knot_ei;
            brk_diff_vec.emplace_back(add_bit_based_on_sign(diff_from_curr_start));
            if(DEBUG_PRINT){
                cout<<"Last entry rank: "<<knot_si - epsilon<<endl;
                cout<<"prev knot ei: "<<prev_knot_ei
                    <<" diff_from_curr: "<<diff_from_curr_start<<endl;
            }
        }
        index_size = brk_kval_vec.size();

        brk_diff_vec.emplace_back(0);// dummy/senitel
                
        y_range_beg.encode(brk_sa_indx_vec.data(), brk_sa_indx_vec.size());
        
        shift_bits = 2* kmer_size;
        
        encode_values(brk_kval_vec, brk_diff_vec, lookup_count);

        if(is_fast_rank) calc_index_bv(begin);
    }

};