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
#include "utils/util.hpp"
#include "utils/compact_vector.hpp"
#include "../../sdsl/dac_vector.hpp"
#include "../../sdsl/sd_vector.hpp"

#include "randstrobes.hpp"
#include "logger.hpp"  

#define DEBUG_PRINT 0
#define KNOT_BS_THRES 64

using namespace std;
using namespace std::chrono;

enum INDX_TYPE {BASIC_PLA, REPEAT_PLA};
static Logger& logger_pla = Logger::get();
class pla_index{
private:
    using canonical_segment = typename OptimalPiecewiseLinearModel<uint64_t, uint64_t>::CanonicalSegment;
    int64_t epsilon, max_data_ratio, lp_bits;
    uint64_t shift_bits, dic_size, total_data_points;
    uint64_t ef_width, x_width, num_dac_levels,
        dac_width, pref_width, x_const_width;
    
    vector<uint64_t> x_enc_vec, prefix_vec;
    vector<uint64_t> dac_level_width_vec;
    vector< vector<uint64_t>> dac_level_vecs;
    arank::ef_sequence<false> y_range_beg;
    sdsl::dac_vector_dp<> sa_diff_dac_vec;
    
    vector<sdsl::sd_vector<> > dac_sd_vecs;

    vector<sdsl::sd_vector<>::rank_1_type> rank_support_dac_vec;
    
    struct brkpnt_values{
        uint64_t x1, ef1, x2_x1,
            brkpnt, is_lvl_sec;
        int64_t diff1;
    };
public:
   // constructor at index build
    pla_index(int64_t eps):
            epsilon(eps)
    {
        x_enc_vec.push_back(0);
        prefix_vec.push_back(0);
    }
    
    // constructor at index query
    pla_index()
        
    {   
    }

    int64_t get_epsilon() const {return epsilon;}

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

    void process_level_vec(){
        
        sa_diff_dac_vec.get_bit_per_level(dac_level_width_vec);
        
        if(DEBUG_PRINT){
            logger_pla.debug()<<"dac_width\n";
            for(size_t i=0; i<dac_level_width_vec.size(); i++) 
                logger_pla.debug()<<dac_level_width_vec[i]<<" ";
            logger_pla.debug()<<endl;
        }
        
        for(size_t i=dac_level_width_vec.size()-1; i>=1; i--){
            dac_level_width_vec[i] -= dac_level_width_vec[i-1];
        }

        if(DEBUG_PRINT){
            logger_pla.debug()<<"dac_width after update\n";
            for(size_t i=0; i<dac_level_width_vec.size(); i++) logger_pla.debug()<<dac_level_width_vec[i]<<" ";
            logger_pla.debug()<<endl;
        }
    }

    // TODO: Check for bugs
    void encode_knots(const vector<uint64_t> &brk_kval_vec, int64_t lookup_count){
        vector<uint64_t> prefix_lookup_vec, level_counter;
        sdsl::bit_vector dac_bv(dic_size);
        // dac_bv.resize(dic_size);

        if(lookup_count > dic_size) lookup_count = 1;

        lp_bits = ceil(log2(dic_size / double(lookup_count))) ;
        shift_bits -= lp_bits; // on build_index, shift_bits = total_bits
        
        uint64_t max_bits = shift_bits + lp_bits;
        ef_width = y_range_beg.cv_width();
        const uint64_t ef_mask = (1ULL << ef_width) - 1;
        
        dac_width = 1; // 1st level
        const uint64_t dac_mask = 1;

        if(DEBUG_PRINT){
            logger_pla.debug()<<"Shift bits: "<<shift_bits<<" lp bits: "<<lp_bits<<endl;
            logger_pla.debug()<<"ef width: "<<ef_width<<endl;
        }
        
        // TODO: Get rid of the x vec encoding, if it needs 64 bits anyway
        if(max_bits == 64) max_data_ratio = ~UINT64_C(0)/(dic_size - 1);
        else max_data_ratio = ((1ULL << (max_bits)) - 1)/(dic_size - 1);
        // max_data_ratio = ((1ULL << (max_bits)) - 1)/(dic_size - 1);
        // TODO: SIMD [lowest priority]
        uint64_t largest = 0;
        for(int64_t i=0; i<dic_size; i++){
            int64_t diff = __int128_t(brk_kval_vec[i]) - __int128_t(i*__int128_t(max_data_ratio));
            uint64_t val = add_bit_based_on_sign(diff);
            // x_diff_vec.emplace_back(val);
            if(val > largest) {
                largest = val;
                // logger_pla.debug()<<i<<" "<<(brk_kval_vec[i]>>1)<<" "<<i*max_data_ratio<<" "<<largest<<endl;
            }
        }
        if(ceil(log2(largest)) == floor(log2(largest))) x_width = ceil(log2(largest)) + 1;
        else x_width = ceil(log2(largest));
        if(DEBUG_PRINT){
            logger_pla.debug()<<"Max data raio: "<<max_data_ratio
                <<" Largest: "<<largest
                <<" x_width: "<<x_width<<endl;
        }
        
        // TODO: Correct this in the main code!!
        const uint64_t x_mask = (1ULL << x_width) -1;


        if(ceil(log2(dic_size-1)) == floor(log2(dic_size-1))) pref_width = ceil(log2(dic_size-1)) + 1;
        else pref_width = ceil(log2(dic_size-1));

        // x_mask = (1ULL << x_width) - 1;
        const uint64_t pref_mask = (1ULL << pref_width) - 1;

        // gets the bit width at each level, and the corresponding mask 
        process_level_vec();
        vector<sdsl::bit_vector> dac_bvs(dac_level_width_vec.size()); // TODO: size() - 1
        dac_sd_vecs.resize(sa_diff_dac_vec.levels()-1);
        
        // no need for last level, but included for now
        uint64_t curr_size = 0;
        for(int64_t i=sa_diff_dac_vec.levels()-1; i>=0; i--){
            curr_size += sa_diff_dac_vec.get_size_at_level(i);
            if(i<num_dac_levels - 1)
                dac_bvs[i].resize(curr_size);
        }
        
        vector<uint64_t> dac_level_block_vec(sa_diff_dac_vec.levels()), 
                    dac_level_shift_vec(sa_diff_dac_vec.levels()),
                    dac_level_counter(sa_diff_dac_vec.levels(), 0);
        
        dac_level_vecs.resize(sa_diff_dac_vec.levels());
        dac_level_block_vec[0] = 0;
        dac_level_shift_vec[0] = 0;
        dac_level_vecs[0].push_back(0);
        uint64_t total_bits = 0;
        for(size_t i=1; i<sa_diff_dac_vec.levels(); i++){
            // total_bits += sa_diff_dac_vec.get_bit_at_level(i-1)*sa_diff_dac_vec.get_size_at_level(i-1);
            // dac_level_block_vec[i] = total_bits >> 6;
            // dac_level_shift_vec[i] = total_bits & 63;
            dac_level_block_vec[i] = 0;
            dac_level_shift_vec[i] = 0;
            dac_level_vecs[i].push_back(0);
        }

        // vector<uint64_t> level_max_bit(dac_level_width_vec.size()); // can also be done with a select call I think
        
        // x + bit_for_2nd_level + ef_width + 1st level dac
        x_const_width = x_width + 1 + ef_width + dac_level_width_vec[0];
        
        uint64_t curr_kval_position, prev_level = 0, prefix_count=0;
        uint64_t x_enc_block = 0, x_enc_shift = 0, prefix_enc_block=0, prefix_enc_shift=0,
            dac_enc_block = 0, dac_enc_shift = 0;
        uint64_t slope_width = 1, slope_mask = 1;
        const uint64_t dac_0_mask = (1ULL << dac_level_width_vec[0]) - 1;
        if(DEBUG_PRINT)
            logger_pla.debug()<<"DAC mask: "<<dac_mask<<" width: "<<dac_width<<endl;
        for(size_t i=0; i<dic_size; i++){
            // enc x value
            int64_t diff = __int128_t(brk_kval_vec[i]) - __int128_t(i*__int128_t(max_data_ratio));
            uint64_t val = add_bit_based_on_sign(diff);
            uint64_t curr_level = 0;
            
            enc_value(val, x_enc_vec, x_mask, x_width, x_enc_block, x_enc_shift);
            
            // enc slope
            // enc_value(brk_kval_vec[i]&1, x_enc_vec, slope_mask, slope_width, x_enc_block, x_enc_shift);
            
            //enc whether to go to 2nd level
            auto [level, element] = sa_diff_dac_vec.access_element(i);
            uint64_t enc_level = level > curr_level ? 1:0;
            
            enc_value(enc_level, x_enc_vec, dac_mask, dac_width, x_enc_block, x_enc_shift);
            
            dac_bvs[curr_level][i] = enc_level; // either 0 or 1

            // enc ef
            enc_value(y_range_beg.cv_access(i), x_enc_vec, ef_mask, ef_width, x_enc_block, x_enc_shift);
            
            // encode lsbs of all level dac
            enc_value(sa_diff_dac_vec[i]&dac_0_mask, x_enc_vec, 
                dac_0_mask, dac_level_width_vec[0], x_enc_block, x_enc_shift);
            
            if(level > 0){
                int64_t val = sa_diff_dac_vec[i] >> dac_level_width_vec[curr_level];
                while(1){
                    curr_level++;
                    enc_level = level > curr_level ? 1:0;
                    
                    uint64_t dac_mask = (1ULL<<dac_level_width_vec[curr_level])-1;
                    enc_value(val & dac_mask, dac_level_vecs[curr_level-1], dac_mask,
                        dac_level_width_vec[curr_level], dac_level_block_vec[curr_level-1], dac_level_shift_vec[curr_level-1]);
                    
                    if(curr_level < num_dac_levels-1)
                        dac_bvs[curr_level][dac_level_counter[curr_level]++] = enc_level;

                    if(!enc_level) break;
                    val >>= dac_level_width_vec[curr_level];
                }

                // enc_value(sa_diff_dac_vec[i] >> dac_level_width_vec[0], dac_lvl2_vec, dac_level_mask_vec[1],
                //     dac_level_width_vec[1], dac_enc_block, dac_enc_shift);
            }
            if(DEBUG_PRINT && i>=8139467 && i<=8139470){
                logger_pla.debug()<<i<<" x: "<<(brk_kval_vec[i])<<" y_start: "
                    <<y_range_beg.access(i, y_range_beg.cv_access(i), ef_width)<<" "
                    <<" cv: "<<y_range_beg.cv_access(i)<<" lvl: "
                    <<level<<" diff1: "
                    <<sa_diff_dac_vec[i]<<" y_end: "
                    <<y_range_beg.access(i, y_range_beg.cv_access(i), ef_width) + sa_diff_dac_vec[i]
                    <<endl;
            }
            int64_t lp_bit_value = (brk_kval_vec[i]) >> shift_bits;
            if(prefix_count <= lp_bit_value){
                for(int64_t pl=prefix_count; pl<lp_bit_value+1; pl++){
                    prefix_count++;
                    enc_value(i, prefix_vec, pref_mask, pref_width, prefix_enc_block, prefix_enc_shift);
                }
            }
        }
        // this ensures upto dic_size-1 will be checked for the last entry
        enc_value(dic_size-1, prefix_vec, pref_mask, pref_width, prefix_enc_block, prefix_enc_shift);

        rank_support_dac_vec.resize(dac_sd_vecs.size());
        for(size_t i=0; i<dac_sd_vecs.size(); i++){
            sdsl::sd_vector<> sd(dac_bvs[i]);
            dac_sd_vecs[i] = sd;
            rank_support_dac_vec[i].set_vector(&dac_sd_vecs[i]);
        }
        x_const_width = x_width + 1 + dac_level_width_vec[0] + ef_width;
    }
    
    
    
    // void save_bit_info(string indx_fn, CBitPacking &ind_diff_pack);
    // uint64_t get_size_in_bytes();
    inline uint64_t get_num_segments() const {return dic_size;}

    
    void Save(std::ofstream &os) const{
        uint8_t t_lp_bits = lp_bits;
        essentials_arank::save_pod(os, t_lp_bits);
        
        uint8_t pref_width_t = pref_width;
        essentials_arank::save_pod(os, pref_width_t);
        essentials_arank::save_vec(os, prefix_vec);        
        logger_pla.debug()<<"prefix vec size: "<<prefix_vec.size()<<endl;
        
        uint32_t nBrkpnts = dic_size;
        essentials_arank::save_pod(os, nBrkpnts);

        uint8_t x_width_t = x_width;
        essentials_arank::save_pod(os, x_width_t);
        essentials_arank::save_vec(os, x_enc_vec);
        
        essentials_arank::save_pod(os, num_dac_levels);
        essentials_arank::save_vec(os, dac_level_width_vec);
        for(size_t i=0; i<num_dac_levels; i++){
            essentials_arank::save_vec(os, dac_level_vecs[i]);
            if(i < num_dac_levels-1)
                dac_sd_vecs[i].serialize(os);
        }
        
        // y_diff
        uint8_t ef_width_t = ef_width;
        essentials_arank::save_pod(os, ef_width_t);
        y_range_beg.save(os);
        
        // meta info
        uint16_t t_err_th = epsilon;
        essentials_arank::save_pod(os, t_err_th);
    }
    
    // suffix array version
    
    void Load(std::ifstream& is, const vector<RefRandstrobe> &randstrobes){
        total_data_points = randstrobes.size();
        uint64_t largest = randstrobes[total_data_points-1].hash;
        uint64_t total_bits;
        if(ceil(log2(largest)) == floor(log2(largest))) total_bits = ceil(log2(largest)) + 1;
        else total_bits = ceil(log2(largest));

        if(total_bits > 64) total_bits = 64;
                
        // prefix lookup table info
        uint8_t lp_bits_t;
        essentials_arank::load_pod(is, lp_bits_t);
        lp_bits = lp_bits_t;
        
        shift_bits = total_bits - lp_bits;

        uint8_t pref_width_t;
        essentials_arank::load_pod(is, pref_width_t);
        pref_width = pref_width_t;
        // pref_mask = (1ULL << pref_width) - 1;
        
        essentials_arank::load_vec(is, prefix_vec);
        
        // x_vec
        uint32_t dic_size_t;
        essentials_arank::load_pod(is, dic_size_t);
        dic_size = dic_size_t;
        
        uint8_t x_width_t;
        essentials_arank::load_pod(is, x_width_t);
        x_width = x_width_t;
        // x_mask = (1ULL << x_width) - 1;
        // logger_pla.debug()<<"x width: "<<x_width<<" "<<x_mask<<endl;
        essentials_arank::load_vec(is, x_enc_vec);
        
        
        essentials_arank::load_pod(is, num_dac_levels);
        essentials_arank::load_vec(is, dac_level_width_vec);
        
        
        dac_level_vecs.resize(num_dac_levels);
        dac_sd_vecs.resize(num_dac_levels);
        rank_support_dac_vec.resize(num_dac_levels);

        for(int64_t i=0; i<num_dac_levels; i++){
            essentials_arank::load_vec(is, dac_level_vecs[i]);
            if(i < num_dac_levels-1){
                dac_sd_vecs[i].load(is);
                rank_support_dac_vec[i].set_vector(&dac_sd_vecs[i]);
            }
            
        }
        
        
        max_data_ratio = ((1ULL << total_bits) - 1)/(dic_size-1);
        
        uint8_t ef_width_t;
        essentials_arank::load_pod(is, ef_width_t);
        ef_width = ef_width_t;
        
        y_range_beg.load(is);
        
        x_const_width = x_width + 1 + dac_level_width_vec[0] + ef_width;

        uint16_t max_error_t;
        essentials_arank::load_pod(is, max_error_t);
        epsilon = (int64_t)max_error_t;
    }

    inline int64_t get_signed_value(int64_t value) const{
        return value & 1
            ? -(value>>1)
            : value>>1;
    }

    

    inline int64_t predict(uint64_t q, uint64_t x1, int64_t y1, double x2_x1, int64_t y2) const{
        // if(DEBUG_PRINT)
        //     logger_pla.debug()<<"Predict: y1: "<<y1<<" q-x1: "<<(q - x1)<<
        //         " x2_x1: "<<x2_x1<<" y2-y1: "<<(y2 - y1)<<" slope: "<<((double)(q - x1)/x2_x1 )*(y2 - y1)<<endl;
        return round(y1 + ((q - x1)/x2_x1 )*(y2 - y1));
    }

    // TODO: Feels like there is something that can be done to make this efficient
    /**
     * Gets the value from width bits at pos, also increments pos to width value
     * 
     * Check out the bextr instruction
    */
    inline uint64_t get_val(uint64_t& pos, const uint64_t& mask, const uint64_t& width) const{
        // pos += width;
        // return arank::util::readInt(x_enc_vec.data(), pos-width, width);
        uint64_t block = pos >> 6;
        uint64_t shift = pos & 63;
        pos += width;
        return shift + width <= 64
                    ? x_enc_vec[block] >> shift & mask
                    : (x_enc_vec[block] >> shift) | (x_enc_vec[block + 1] << (64 - shift) & mask);
    }
    
    
    // TODO: Do timing without using inline
    inline void binary_search(uint64_t lo, uint64_t hi, const uint64_t query_val, 
            brkpnt_values &bv) const
    {
        int64_t mid, pos,
            block, shift, rest_bits = x_const_width - x_width;
        uint64_t curr_x=0, prev_x = 0, midVal, rest_val, mdr;
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
                    //     logger_pla.debug()<<i<<" "<<lo<<" "<<hi<<" "<<x_total_width<<" "<<pos<<" "
                    //     <<curr_x<<" "<<query_val<<endl;
                    // }
                    
                    if(curr_x > query_val) {
                        if(prev_x != 0){
                            bv.x1 = prev_x;
                            pos -= rest_bits;
                            block = pos >> 6;
                            shift = pos & 63;
                            rest_val = shift + rest_bits <= 64
                                ? x_enc_vec[block] >> shift & rest_mask
                                : (x_enc_vec[block] >> shift) | (x_enc_vec[block + 1] << (64 - shift) & rest_mask);
                        
                            // bv.is_next_slope = rest_val & 1;
                            bv.is_lvl_sec = rest_val & 1;
                            
                            bv.ef1 = (rest_val>>1) & ef_mask;
                            bv.diff1 = rest_val >> (ef_width + 1);
                            // bv.x2_x1 = curr_x - prev_x;
                            // if(DEBUG_PRINT){
                            //     logger_pla.debug()<<"In bs x: sec lvl: "<<bv.is_lvl_sec<<endl;
                            //     logger_pla.debug()<<"rest val: "<<rest_val<<" "<<rest_bits<<" "<<(1ULL << (rest_bits-2))<<endl;
                            //     logger_pla.debug()<<"diff: "<<bv.diff1<<" ef "<<bv.ef1<<endl;
                            //     logger_pla.debug()<<"x2_x1: "<<curr_x - bv.x1
                            //     <<" brkpnt: "<<i-1
                            //     <<endl;
                            //     // logger_pla.debug()<<"block "<<block<<" "<<x_enc_vec[block]<<" "<<shift<<endl;
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
                // this case only happens for the last k-mer
                // match or no match, takes to the last entry
                // knot_ei will take to the correct position
                // if(curr_x == query_val){
                bv.x1 = curr_x;
                bv.x2_x1 = 1; // actually 0, set 1 for not zero division (dummy)
                pos -= rest_bits;
                block = pos >> 6;
                shift = pos & 63;
                // curr_x = shift + x_width <= 64
                //     ? x_enc_vec[block] >> shift & x_mask
                //     : (x_enc_vec[block] >> shift) | (x_enc_vec[block + 1] << (64 - shift) & x_mask);
                // curr_x = get_signed_value(curr_x) + mdr-(max_data_ratio<<1);
                
                // block = pos >> 6;
                // shift = pos & 63;
                rest_val = shift + rest_bits <= 64
                    ? x_enc_vec[block] >> shift & rest_mask
                    : (x_enc_vec[block] >> shift) | (x_enc_vec[block + 1] << (64 - shift) & rest_mask);
            
                // bv.is_next_slope = rest_val & 1;
                bv.is_lvl_sec = 0;
                
                bv.ef1 = (rest_val>>1) & ef_mask;
                bv.diff1 = 0;
                bv.brkpnt = hi;
                // if(DEBUG_PRINT)
                //     logger_pla.debug()<<"match: "<<curr_x<<" "<<bv.x2_x1<<" "<<rest_val<<endl;
                return;
                // }
                // logger_pla.debug()<<"Not found\n";
                // logger_pla.debug()<<query_val<<endl;
                // exit(-1);
            } 
            mid = (lo + hi) >> 1;
            pos = mid * x_total_width;
            mdr = max_data_ratio*mid;
            block = pos >> 6;
            shift = pos & 63;
            midVal = shift + x_width <= 64
                ? x_enc_vec[block] >> shift & x_mask
                : (x_enc_vec[block] >> shift) | (x_enc_vec[block + 1] << (64 - shift) & x_mask);
            midVal = get_signed_value(midVal) + mdr;
            // if(DEBUG_PRINT)
            //     logger_pla.debug()<<mid<<" "<<midVal<<" "<<query_val<<endl;
            if (midVal == query_val) {
                pos += x_width;
                block = pos >> 6;
                shift = pos & 63;
                rest_val = shift + rest_bits <= 64
                    ? x_enc_vec[block] >> shift & rest_mask
                    : (x_enc_vec[block] >> shift) | (x_enc_vec[block + 1] << (64 - shift) & rest_mask);
                bv.x1 = midVal;
                // bv.is_next_slope = rest_val & 1;
                bv.is_lvl_sec = rest_val & 1;

                bv.ef1 = (rest_val>>1) & ef_mask;
                bv.diff1 = rest_val >> (ef_width+1);
                bv.brkpnt = mid;
                
                pos += rest_bits;
                block = pos >> 6;
                shift = pos & 63;
                bv.x2_x1 = shift + x_width <= 64
                    ? x_enc_vec[block] >> shift & x_mask
                    : (x_enc_vec[block] >> shift) | (x_enc_vec[block + 1] << (64 - shift) & x_mask);
                bv.x2_x1  = get_signed_value(bv.x2_x1) + mdr+ max_data_ratio - midVal;
                // TODO: Update to direct check
                
                // if(DEBUG_PRINT)
                //     logger_pla.debug()<<"ef "<<bv.ef1<<" "<<rest_val<<" "<<bv.diff1<<" "<<bv.x2_x1<<endl;
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

    
    int64_t query(const vector<RefRandstrobe> &randstrobes,
        const uint64_t query_val) const{
        int64_t query_pred, lo,  brkpnt;
        brkpnt_values bv;
        uint64_t lp_bit_value = query_val >> shift_bits;
        const uint64_t ef_mask = (1ULL << ef_width) - 1;
        const uint64_t x_mask = (1ULL << x_width) - 1;
        const uint64_t pref_mask = (1ULL << pref_width) - 1;
        // if(DEBUG_PRINT){
        //     logger_pla.debug()<<"shift: "<<shift_bits<<" "<<lp_bits<<" "<<lp_bit_value<<" "<<query_val<<endl;
        //     logger_pla.debug()<<x_width<<" "<<x_const_width<<" "<<dac_level_width_vec[0]
        //         <<" "<<ef_width<<" "<<pref_width<<endl;
        // }
        
        uint64_t pos = lp_bit_value*pref_width;
        uint64_t block = pos >> 6;
        uint64_t shift = pos & 63;
        uint64_t indx = (shift + pref_width <= 64) 
            ?   prefix_vec[block]>>shift & pref_mask
            :   (prefix_vec[block]>>shift) | (prefix_vec[block+1] << (64-shift) & pref_mask);
        pos += pref_width;
        block = pos >> 6;
        shift = pos & 63;
        uint64_t indx2 = (shift + pref_width <= 64) 
            ?   prefix_vec[block]>>shift & pref_mask
            :   (prefix_vec[block]>>shift) | (prefix_vec[block+1] << (64-shift) & pref_mask);
    
        
        // logger_pla.debug()<<"pref pair: "<<indx<<" "<<indx2<<endl;
        
        
        // TODO: mdr may need 128 bits to compute
        int64_t ef2=0, diff2=0,
            is_sec_lvl_2=0;
        uint64_t x2=0, mdr = max_data_ratio*indx;
        
        uint64_t rest_bits = x_const_width - x_width;
        const uint64_t rest_mask = (1ULL << rest_bits)-1;
        uint64_t rest_val;
        
        // s1 = std::chrono::system_clock::now();
        if(indx > 0){
            pos = (indx-1) * x_const_width;
            x2 = get_signed_value(get_val(pos, x_mask, x_width)) + mdr - max_data_ratio;
            rest_val = get_val(pos, rest_mask, rest_bits);
            is_sec_lvl_2 = rest_val & 1;

            ef2 = (rest_val>>1) & ef_mask;
            diff2 = rest_val >> (ef_width + 1);
            bv.x1 = get_signed_value(get_val(pos, x_mask, x_width)) + mdr;
        }
        else{
            pos = indx * x_const_width;
            bv.x1 = get_signed_value(get_val(pos, x_mask, x_width)) + mdr;
            if(query_val < bv.x1) return -1;
        }
        
        // TODO: there can be intger bit issue here
        
        rest_val = get_val(pos, rest_mask, rest_bits);
        bv.is_lvl_sec = rest_val & 1;
        bv.ef1 = (rest_val>>1) & ef_mask;
        bv.diff1 = rest_val >> (ef_width + 1);

        // s2 = std::chrono::system_clock::now();
        // portion_sec[1] += (s2-s1);
        // if(DEBUG_PRINT){
        //     logger_pla.debug()<<"from prefix: x1: "<<bv.x1<<" query: "<<query_val
        //         <<" x2: (<x1): "<<x2<<" indx1: "<<indx<<" indx2: "<<indx2
        //         <<" bv ef1: "<<bv.ef1<<" diff1: "<<bv.diff1<<" islvlsec: "<<bv.is_lvl_sec
        //         <<endl;
        // }
        
        // s1 = std::chrono::system_clock::now();
        if(bv.x1 == query_val){
            query_pred = y_range_beg.access(indx, bv.ef1, ef_width) - epsilon;
        }
        else{
            if(bv.x1 < query_val){
                // int64_t end_idx = (lp_bit_value == prefix_lookup_cv.size() - 1) ? dic_size - 1 : prefix_lookup_cv[lp_bit_value + 1] - 1;
                // brkpnt_pair = x_cv.binary_search(ind_pair.first, ind_pair.second, query_val, brkpnt);
                // TODO: Can be made indx2-1 and lo+1 w/ careful thinking
                // auto s3 = std::chrono::system_clock::now();
                binary_search(indx+1, indx2, query_val, bv);
                // auto s4 = std::chrono::system_clock::now();
                // portion_sec[3] += (s4-s3);
            }
            else{
                bv.x2_x1 = bv.x1 - x2;
                bv.x1 = x2; // x2 is the previous element
                bv.ef1 = ef2;
                bv.diff1 = diff2;
                bv.is_lvl_sec = is_sec_lvl_2;
                bv.brkpnt = indx - 1;
            }
            int64_t knot_si = y_range_beg.access(bv.brkpnt, bv.ef1, ef_width);
            // if(DEBUG_PRINT)
            //     logger_pla.debug()<<"knot si: "<<knot_si<<endl;
            int64_t knot_ei = bv.diff1;
            
            if(bv.is_lvl_sec){
                int64_t curr_level = 0;
                int64_t total_width = dac_level_width_vec[0];
                int64_t dac_indx = bv.brkpnt;
                while(1){
                    
                    dac_indx = rank_support_dac_vec[curr_level](dac_indx);
                    const int64_t dac_width = dac_level_width_vec[curr_level + 1];
                    const uint64_t dac_mask = (1ULL << dac_width) - 1;
                    
                    pos = dac_indx*dac_width;
                    block = pos >> 6;
                    shift = pos & 63;
                    int64_t temp = (shift + dac_width <= 64) 
                    ?   dac_level_vecs[curr_level][block]>>shift & dac_mask
                    :   (dac_level_vecs[curr_level][block]>>shift) 
                        | (dac_level_vecs[curr_level][block+1] << (64-shift) & dac_mask);
                    
                    knot_ei |= (temp << total_width);
                    
                    total_width += dac_width;
                    // if(DEBUG_PRINT){
                    //     logger_pla.debug()<<query_val<<" level: "<<curr_level<<" knot_ei: "<<knot_ei<<endl;
                    //     logger_pla.debug()<<"dac indx: "<<dac_indx<<" width: "<<dac_width<<" mask: "<<dac_mask<<endl;
                    //     logger_pla.debug()<<"dac temp: "<<temp<<endl;
                    //     logger_pla.debug()<<"total width: "<<total_width<<endl;
                    //     logger_pla.debug()<<"knot ei: "<<knot_ei<<endl;
                    // }
                    if(curr_level == num_dac_levels-2 || !dac_sd_vecs[curr_level+1][dac_indx]) break;
                    curr_level++;
                }
            }
            knot_ei += knot_si;
            // if(DEBUG_PRINT){
            //     logger_pla.debug()<<"After breakpoint: \n"
            //     <<bv.x1<<" "<<query_val<<" "<<bv.x1+bv.x2_x1<<" \n"
            //     <<"si: "<<knot_si<<" ei: "<<knot_ei<<" \n"
            //     <<"brkpnt: "<<bv.brkpnt<<" val_si "<<randstrobes[knot_si-epsilon].hash<<
            //     " val_ei "<<randstrobes[knot_ei - epsilon].hash<<" sec_lvl: "
            //     <<bv.is_lvl_sec<<" "
            //     <<endl<<endl;
            // }
            

            query_pred = predict(query_val, bv.x1, knot_si, bv.x2_x1, knot_ei) - epsilon;
            
        }
        query_pred = std::clamp(query_pred, int64_t(0), int64_t(total_data_points) - 1);
        // s2 = std::chrono::system_clock::now();
        // portion_sec[2] += (s2-s1);
        // if(DEBUG_PRINT){
        //     logger_pla.debug()<<"pred: "<<query_pred<<endl;
        // }
        uint64_t hash = randstrobes[query_pred].hash;
        if(hash == query_val) {
            while(1){
                if(!query_pred) return 0;
                if(randstrobes[--query_pred].hash != query_val) return query_pred+1;
            }
        }
        if(hash < query_val){
            int64_t last_idx = min(int64_t(query_pred + epsilon), int64_t(total_data_points - 1));
            for(int64_t idx = query_pred+1; idx<=last_idx; idx++){
                if(randstrobes[idx].hash == query_val) return idx; // guranteed to be the first index
            }
            return -1;
            
        }
        else{
            int64_t last_idx = max(query_pred - epsilon, int64_t(0));
            for(int64_t idx = query_pred-1; idx>=last_idx; idx--){
                if(randstrobes[idx].hash == query_val) {                
                    while(1){
                        if(!idx) return 0;
                        if(randstrobes[--idx].hash != query_val) return idx+1;
                    }
                }
            }
            return -1;
        }
        
    }

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
                // logger_pla.debug()<<hi<<" "<<lo<<" "<<mid<<endl;
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
            // logger_pla.debug()<<hi<<" "<<lo<<" "<<mid<<endl;
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

    template<typename RandomIt>
    void build_index(const RandomIt &begin, const uint64_t lookup_count){
        build_rep_stretch_pla_index(begin, lookup_count);
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
    
    void build_rep_stretch_pla_index(const vector<RefRandstrobe> &randstrobes, const uint64_t lookup_count){
        vector<uint64_t>  brk_sa_indx_vec, brk_kval_vec;
        vector<int64_t> y_diff_vec;
        bool isFirst = true;
        int64_t prev_knot_ei=0, diff_from_curr_start;
        vector<canonical_segment> out(1);

        uint64_t largest = randstrobes[randstrobes.size()-1].hash;
        if(ceil(log2(largest)) == floor(log2(largest))) shift_bits = ceil(log2(largest)) + 1;
        else shift_bits = ceil(log2(largest));
        total_data_points = randstrobes.size();

        auto in_fun = [&](auto i) { return std::pair<uint64_t, uint64_t>(randstrobes[i].hash, i); };
            
        auto out_fun = [&](auto cs, const int64_t last_x) { 
            out[0] = cs;
            uint64_t curr_x = cs.get_first_x();
            brk_kval_vec.emplace_back(curr_x);
            // insert_to_ind_vec(cs.get_first_x(), brk_kval_vec.size()-1);
            
            auto [knot_si, knot_ei] = cs.get_knot_intersection(last_x, epsilon);      
            brk_sa_indx_vec.emplace_back((uint64_t)knot_si); 
            // for this version, lets keep this, and see what does the space look like
            // calculate the dac storage only, and compare w/ the previous version
            diff_from_curr_start = int64_t(knot_ei) - int64_t(knot_si);
            if(diff_from_curr_start < 0 ){
                logger_pla.debug()<<"neg value: "<<diff_from_curr_start<<endl;
            }
            y_diff_vec.emplace_back(diff_from_curr_start);
            
        };
        auto get_err = [&](int64_t act_idx, int64_t pred_idx){
            return get_pred_diff(act_idx, pred_idx, randstrobes);
        };

        uint64_t seg_count = make_segmentation_rep_pla(total_data_points, epsilon, in_fun, out_fun, get_err);

        // TODO: Correct this in the main code!! brk[] >> 1 -> incorrect version
        if(out[0].get_last_x() != (brk_kval_vec[brk_kval_vec.size()-1])){
            brk_kval_vec.emplace_back(out[0].get_last_x());
            // auto [knot_si, knot_ei] = out[0].get_knot_intersection(out[0].get_last_x(), epsilon);   
            uint64_t last_rank = randstrobes.size()-1;
            while(randstrobes[--last_rank].hash == randstrobes[randstrobes.size()-1].hash);
            uint64_t knot_si = last_rank + 1 + epsilon;
            
            brk_sa_indx_vec.emplace_back(knot_si);
            y_diff_vec.emplace_back(0);
        }
        
        dic_size = brk_kval_vec.size();
        logger_pla.debug()<<"#breakpoints: "<<dic_size<<endl;
        sdsl::dac_vector_dp<> temp(y_diff_vec);
        sa_diff_dac_vec = temp;
        // ofstream os("Dac_storage.bin", std::ios::binary);
        // sa_diff_dac_vec.serialize(os);
        num_dac_levels = sa_diff_dac_vec.levels();
        logger_pla.debug()<<"Number of dac levels: "<<num_dac_levels<<endl;
        
        y_range_beg.encode(brk_sa_indx_vec.data(), brk_sa_indx_vec.size());
        // temp_storage(brk_kval_vec, lookup_count);
        encode_knots(brk_kval_vec, lookup_count);
        logger_pla.debug()<<"Done encoding\n";
    }
};