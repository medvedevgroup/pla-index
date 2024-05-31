/**
 * Source: https://github.com/jermp/pthash/
*/

#pragma once

#include <vector>
#include <cmath>

#include "essentials.hpp"

using namespace std;

namespace arank {

struct compact_vector {
    template <typename Data>
    struct enumerator {
        enumerator() {}

        enumerator(Data const* data, uint64_t i = 0)
            : m_i(i)
            , m_cur_val(0)
            , m_cur_block((i * data->m_width) >> 6)
            , m_cur_shift((i * data->m_width) & 63)
            , m_data(data) {}

        uint64_t operator*() {
            read();
            return m_cur_val;
        }

        enumerator& operator++() {
            ++m_i;
            return *this;
        }

        inline uint64_t value() {
            read();
            return m_cur_val;
        }

        inline void next() {
            ++m_i;
        }

        bool operator==(enumerator const& other) const {
            return m_i == other.m_i;
        }

        bool operator!=(enumerator const& other) const {
            return !(*this == other);
        }

    private:
        uint64_t m_i;
        uint64_t m_cur_val;
        uint64_t m_cur_block;
        int64_t m_cur_shift;
        Data const* m_data;

        void read() {
            if (m_cur_shift + m_data->m_width <= 64) {
                m_cur_val = m_data->m_bits[m_cur_block] >> m_cur_shift & m_data->m_mask;
            } else {
                uint64_t res_shift = 64 - m_cur_shift;
                m_cur_val = (m_data->m_bits[m_cur_block] >> m_cur_shift) |
                            (m_data->m_bits[m_cur_block + 1] << res_shift & m_data->m_mask);
                ++m_cur_block;
                m_cur_shift = -res_shift;
            }

            m_cur_shift += m_data->m_width;

            if (m_cur_shift == 64) {
                m_cur_shift = 0;
                ++m_cur_block;
            }
        }
    };

    struct builder {
        builder() : m_size(0), m_width(0), m_mask(0), m_back(0), m_cur_block(0), m_cur_shift(0) {}

        builder(uint64_t n, uint64_t w) {
            resize(n, w);
        }

        void resize(size_t n, uint64_t w) {
            m_size = n;
            m_width = w;
            m_mask = -(w == 64) | ((uint64_t(1) << w) - 1);
            m_back = 0;
            m_cur_block = 0;
            m_cur_shift = 0;
            m_bits.resize(
                /* use 1 word more for safe access() */
                essentials_arank::words_for(m_size * m_width) + 1, 0);
        }

        template <typename Iterator>
        builder(Iterator begin, uint64_t n, uint64_t w) : builder(n, w) {
            fill(begin, n);
        }

        template <typename Iterator>
        void fill(Iterator begin, uint64_t n) {
            if (!m_width) throw std::runtime_error("width must be greater than 0");
            for (uint64_t i = 0; i != n; ++i, ++begin) push_back(*begin);
        }

        void set(uint64_t i, uint64_t v) {
            assert(m_width);
            assert(i < m_size);
            if (i == m_size - 1) m_back = v;

            uint64_t pos = i * m_width;
            uint64_t block = pos >> 6;
            uint64_t shift = pos & 63;

            m_bits[block] &= ~(m_mask << shift);
            m_bits[block] |= v << shift;

            uint64_t res_shift = 64 - shift;
            if (res_shift < m_width) {
                m_bits[block + 1] &= ~(m_mask >> res_shift);
                m_bits[block + 1] |= v >> res_shift;
            }
        }

        void push_back(uint64_t v) {
            assert(m_width);
            m_back = v;
            m_bits[m_cur_block] &= ~(m_mask << m_cur_shift);
            m_bits[m_cur_block] |= v << m_cur_shift;

            uint64_t res_shift = 64 - m_cur_shift;
            if (res_shift < m_width) {
                ++m_cur_block;
                m_bits[m_cur_block] &= ~(m_mask >> res_shift);
                m_bits[m_cur_block] |= v >> res_shift;
                m_cur_shift = -res_shift;
            }

            m_cur_shift += m_width;

            if (m_cur_shift == 64) {
                m_cur_shift = 0;
                ++m_cur_block;
            }
        }

        friend struct enumerator<builder>;

        typedef enumerator<builder> iterator;

        iterator begin() const {
            return iterator(this);
        }

        iterator end() const {
            return iterator(this, size());
        }

        void build(compact_vector& cv) {
            cv.m_size = m_size;
            cv.m_width = m_width;
            cv.m_mask = m_mask;
            cv.m_bits.swap(m_bits);
            builder().swap(*this);
        }

        void swap(compact_vector::builder& other) {
            std::swap(m_size, other.m_size);
            std::swap(m_width, other.m_width);
            std::swap(m_mask, other.m_mask);
            std::swap(m_cur_block, other.m_cur_block);
            std::swap(m_cur_shift, other.m_cur_shift);
            m_bits.swap(other.m_bits);
        }

        uint64_t back() const {
            return m_back;
        }

        uint64_t size() const {
            return m_size;
        }

        uint64_t width() const {
            return m_width;
        }

        std::vector<uint64_t>& bits() {
            return m_bits;
        }

    private:
        uint64_t m_size;
        uint64_t m_width;
        uint64_t m_mask;
        uint64_t m_back;
        uint64_t m_cur_block;
        int64_t m_cur_shift;
        std::vector<uint64_t> m_bits;
    };

    compact_vector() : m_size(0), m_width(0), m_mask(0) {}

    template <typename Iterator>
    void build(Iterator begin, uint64_t n) {
        assert(n > 0);
        uint64_t max = *std::max_element(begin, begin + n);
        uint64_t width = max == 0 ? 1 : std::ceil(std::log2(max + 1));
        build(begin, n, width);
    }

    template <typename Iterator>
    void build(Iterator begin, uint64_t n, uint64_t w) {
        compact_vector::builder builder(begin, n, w);
        builder.build(*this);
    }

    inline uint64_t operator[](uint64_t i) const {
        assert(i < size());
        uint64_t pos = i * m_width;
        uint64_t block = pos >> 6;
        // __builtin_prefetch(m_bits.data()+block);
        // __builtin_prefetch(m_bits.data()+block+1);
        uint64_t shift = pos & 63;
        return shift + m_width <= 64
                   ? m_bits[block] >> shift & m_mask
                   : (m_bits[block] >> shift) | (m_bits[block + 1] << (64 - shift) & m_mask);
    }

    // get ((x[i], y[i]), (x[i-1], y[i-1]))
    inline pair<pair<int64_t, uint64_t>, pair<int64_t, uint64_t> > get_x_pair(const int64_t &i) const{
        uint64_t pos = i * m_width;
        int64_t ratio = mdr * i;
        // cout<<"mdr: "<<mdr<<endl;
        uint64_t block = pos >> 6;
        uint64_t shift = pos & 63;
        uint64_t val1 = shift + m_width <= 64
                   ? m_bits[block] >> shift & m_mask
                   : (m_bits[block] >> shift) | (m_bits[block + 1] << (64 - shift) & m_mask);
        // int64_t x_val1 = get_signed_value(val1 << ef_bits) + ratio;
        // uint64_t ef_val1 = val1 & ef_mask;
        
        // int64_t x_val2 = 0; // these are dummy values; if i=0, first condition will be checked in query
        // uint64_t ef_val2 = 0;
        uint64_t val2;
        if(i > 0){
            pos -= m_width;
            // ratio -= mdr;
            block = pos >> 6;
            shift = pos & 63;
            val2 = shift + m_width <= 64
                   ? m_bits[block] >> shift & m_mask
                   : (m_bits[block] >> shift) | (m_bits[block + 1] << (64 - shift) & m_mask);
            // x_val2 = get_signed_value(val2 << ef_bits) + ratio;
            // ef_val2 = val2 & ef_mask;
        }
        return make_pair(make_pair(get_signed_value(val1 >> ef_bits) + ratio, val1 & ef_mask),
                make_pair(get_signed_value(val2 >> ef_bits) + ratio - mdr, val2 & ef_mask) );
    }

    // get (prefix[i], prefix[i+1])
    inline pair<int64_t, int64_t> get_ind_pair(const int64_t &i, uint64_t &indx_size) const{
        uint64_t pos = i * m_width;
        uint64_t block = pos >> 6;
        uint64_t shift = pos & 63;
        int64_t val1 = shift + m_width <= 64
                   ? m_bits[block] >> shift & m_mask
                   : (m_bits[block] >> shift) | (m_bits[block + 1] << (64 - shift) & m_mask);
        int64_t val2 = indx_size;
        if(i < size()-1){
            pos += m_width;
            block = pos >> 6;
            shift = pos & 63;
            val2 = shift + m_width <= 64
                   ? m_bits[block] >> shift & m_mask
                   : (m_bits[block] >> shift) | (m_bits[block + 1] << (64 - shift) & m_mask);
        }
        return make_pair(val1, val2-1);
    }

    inline int64_t get_signed_value(int64_t value) const{
        return value & 1
            ? -(value>>1)
            : value>>1;
    }

    int64_t get_x_val(int64_t i) const{
        return get_signed_value((*this)[i] >> ef_bits) + mdr*i;
    }
    
    // returns ((x[i], ef_cv[i]), (x[i+1] - x[i], ef_cv[i+1])), also sets the brkpnt at parameter value
    pair<pair<int64_t, uint64_t>, pair<int64_t, uint64_t> > binary_search(int64_t lo, int64_t hi, 
                const int64_t kval, int64_t &brkpnt) const{
        int64_t mid, ratio, prev_x_val = -1, x_val;
        uint64_t pos, block, shift, value, prev_value;
        // cout<<lo<<" "<<ef_bits<<" "<<hi<<endl;
        // cout<<get_x_val(lo)<<" "<<kval<<" "<<get_x_val(hi)<<endl;
        // int _count=0;
        while(1){
            if(hi-lo <= 64){
                // cout<<lo<<" "<<ef_bits<<" "<<hi<<endl;
                // cout<<get_x_val(lo)<<" "<<kval<<" "<<get_x_val(hi)<<endl;
                ratio = mdr * hi;
                pos = hi * m_width;
                for (int64_t i = hi; i >= lo; i--)
                {
                    block = pos >> 6;
                    shift = pos & 63;
                    
                    value = (shift + m_width <= 64) ?
                        m_bits[block] >> shift & m_mask
                    :   (m_bits[block] >> shift) | (m_bits[block + 1] << (64 - shift) & m_mask);
                    
                    x_val = get_signed_value(value >> ef_bits) + ratio;
                    // cout<<i<<" "<<_count<<" "<<x_val<<" "<<kval<<endl;
                    if(x_val <= kval){
                        brkpnt = i;
                        if(prev_x_val != -1){
                            return make_pair(make_pair(x_val, value & ef_mask), 
                                    make_pair(prev_x_val - x_val, prev_value & ef_mask));
                        }
                        else{
                            ratio += mdr;
                            pos += m_width;
                            block = pos >> 6;
                            shift = pos & 63;
                            prev_value = (shift + m_width <= 64) ?
                                    m_bits[block] >> shift & m_mask
                                :   (m_bits[block] >> shift) | (m_bits[block + 1] << (64 - shift) & m_mask);
                            prev_x_val = get_signed_value(prev_value >> ef_bits) + ratio;
                            return make_pair(make_pair(x_val, value & ef_mask), 
                                    make_pair(prev_x_val - x_val, prev_value & ef_mask));

                        }
                    }
                    ratio -= mdr;
                    prev_value = value;
                    prev_x_val = x_val;
                    pos -= m_width;
                    // _count++;
                    // if(_count > 64) break;
                }
                // cout<<_count<<endl;
                // exit(1);
                
            } 
            mid = (lo + hi) >> 1;
            value = (*this)[mid];
            x_val = get_signed_value(value >> ef_bits) + mdr*mid;
            if (x_val == kval) {
                brkpnt = mid;
                uint64_t prev_value = (*this)[mid+1];
                return make_pair(make_pair(x_val, value & ef_mask),  
                        make_pair(get_signed_value(prev_value >> ef_bits)+ mdr * mid + mdr - x_val, prev_value & ef_mask));
            }
            else if(x_val > kval){ 
                hi = mid;
            }
            else{ 
                lo = mid;
            }
        }  
    }


    /**
     * From pos value we get the ith bits position where the data is at
     * Then, i>>3 -> gives the 64 bit alignment where the data is located
     * For example width is 4, pos is 30, so data is at 120.
     * So, so the 64 bit data retrieved will have value from [120, 123] bits
     * ptr + i>>3  is char pointer operations, so ptr + 15 byte approach
     * 15 *8 = 120
     * so, that 
     * but we dont need 124-127 bits
     * 120 % 8 = 0
    */

    // it retrieves at least 57 bits
    inline uint64_t access(uint64_t pos) const {
        assert(pos < size());
        uint64_t i = pos * m_width;
        const char* ptr = reinterpret_cast<const char*>(m_bits.data());
        return (*(reinterpret_cast<uint64_t const*>(ptr + (i >> 3))) >> (i & 7)) & m_mask;
    }

    inline pair<uint64_t, uint64_t> get_two_consecutive_access(uint64_t pos) const{
        uint64_t i = pos * m_width;
        uint64_t j = i + m_width;
        const char* ptr = reinterpret_cast<const char*>(m_bits.data());
        return make_pair((*(reinterpret_cast<uint64_t const*>(ptr + (i >> 3))) >> (i & 7)) & m_mask,
            (*(reinterpret_cast<uint64_t const*>(ptr + (j >> 3))) >> (j & 7)) & m_mask);
    }

    inline uint64_t get_prev_unmatched(int64_t &curr_val, uint64_t &idx){
        uint64_t pos, block, shift, value;
        pos = idx * m_width;
        do{
            block = pos >> 6;
            shift = pos & 63;
            value = (shift + m_width <= 64)
                   ? m_bits[block] >> shift & m_mask
                   : (m_bits[block] >> shift) | (m_bits[block + 1] << (64 - shift) & m_mask);
            pos -= m_width;
        }while(value == curr_val);
        return value;        
    }

    uint64_t back() const {
        return operator[](size() - 1);
    }

    inline uint64_t size() const {
        return m_size;
    }

    inline uint64_t width() const {
        return m_width;
    }

    typedef enumerator<compact_vector> iterator;

    iterator begin() const {
        return iterator(this);
    }

    iterator end() const {
        return iterator(this, size());
    }

    iterator at(uint64_t pos) const {
        return iterator(this, pos);
    }

    std::vector<uint64_t> const& bits() const {
        return m_bits;
    }

    size_t bytes() const {
        return sizeof(m_size) + sizeof(m_width) + sizeof(m_mask) + essentials_arank::vec_bytes(m_bits);
    }

    void swap(compact_vector& other) {
        std::swap(m_size, other.m_size);
        std::swap(m_width, other.m_width);
        std::swap(m_mask, other.m_mask);
        m_bits.swap(other.m_bits);
    }

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visitor.visit(m_size);
        visitor.visit(m_width);
        visitor.visit(m_mask);
        visitor.visit(m_bits);
    }

    void save(std::ostream &os){
        essentials_arank::save_pod(os, m_size);
        essentials_arank::save_pod(os, m_width);
        essentials_arank::save_pod(os, m_mask);
        essentials_arank::save_vec(os, m_bits);
    }

    void load(std::istream &is){
        essentials_arank::load_pod(is, m_size);
        essentials_arank::load_pod(is, m_width);
        essentials_arank::load_pod(is, m_mask);
        essentials_arank::load_vec(is, m_bits);
    }

    void save_fp(FILE* fp){
        fwrite(&m_size, sizeof(m_size), 1, fp);
        fwrite(&m_width, sizeof(m_width), 1, fp);
        fwrite(&m_mask, sizeof(m_mask), 1, fp);
        essentials_arank::save_vec_fp(fp, m_bits);
    }

    void load_fp(FILE* fp){
        uint64_t err = fread(&m_size, sizeof(m_size), 1, fp);
        err = fread(&m_width, sizeof(m_width), 1, fp);
        err = fread(&m_mask, sizeof(m_mask), 1, fp);
        essentials_arank::load_vec_fp(fp, m_bits);
        if(!err) std::cout<<"Read 0 bytes in compact vec: load_fp. should not happen\n";
    }

    int64_t mdr;
    uint64_t ef_bits, ef_mask;

private:
    uint64_t m_size;
    uint64_t m_width;
    uint64_t m_mask;
    std::vector<uint64_t> m_bits;
};

}  // namespace arank