/**
 * Source: https://github.com/jermp/pthash/
*/

#pragma once

#include "bit_vector.hpp"
#include "darray.hpp"
#include "compact_vector.hpp"

namespace arank {

template <bool encode_prefix_sum = false>
struct ef_sequence{
    ef_sequence() {}

    template <typename Iterator>
    void encode(Iterator begin, uint64_t n) {
        if (n == 0) return;

        uint64_t u = *(begin + n - 1);
        if constexpr (encode_prefix_sum) {
            u = std::accumulate(begin, begin + n, static_cast<uint64_t>(0));
            n = n + 1;  // because I will add a zero at the beginning
        }

        uint64_t l = uint64_t((n && u / n) ? util::msb(u / n) : 0);
        // std::cout<<"u value: "<<u<<" n: "<<n<<endl;
        // std::cout<<"l value: "<<l<<endl;
        bit_vector_builder bvb_high_bits(n + (u >> l) + 1);
        compact_vector::builder cv_builder_low_bits(n, l);

        uint64_t low_mask = (uint64_t(1) << l) - 1;
        uint64_t last = 0;
        // I add a zero at the beginning
        if constexpr (encode_prefix_sum) {
            if (l) cv_builder_low_bits.push_back(0);
            bvb_high_bits.set(0, 1);
            n = n - 1;  // restore n
        }
        for (size_t i = 0; i < n; ++i, ++begin) {
            auto v = *begin;
            if constexpr (encode_prefix_sum) {
                v = v + last;             // prefix sum
            } else if (i and v < last) {  // check the order
                std::cerr << "error at " << i << "/" << n << ":\n";
                std::cerr << "last " << last << "\n";
                std::cerr << "current " << v << "\n";
                throw std::runtime_error("ef_sequence is not sorted");
            }
            if (l) cv_builder_low_bits.push_back(v & low_mask);
            bvb_high_bits.set((v >> l) + i + encode_prefix_sum, 1);
            last = v;
        }

        bit_vector(&bvb_high_bits).swap(m_high_bits);
        cv_builder_low_bits.build(m_low_bits);
        m_high_bits_d1.build(m_high_bits);
    }

    

    inline uint64_t access(uint64_t i, uint64_t ef_cv_val, uint64_t ef_width) const {
        assert(i < size());
        return ((m_high_bits_d1.select(m_high_bits, i) - i) << ef_width) |
                ef_cv_val;
        //       m_low_bits.access(i);
    }

    inline uint64_t cv_access(uint64_t i) const{
        // TODO: Remove this condition
        // if(m_low_bits.access(i) != m_low_bits[i]){
        //     std::cerr<<"Access and [] does not return same value\n";
        //     exit(EXIT_FAILURE);
        // }
        return m_low_bits.access(i);
    }

    inline uint64_t cv_width() const{
        return m_low_bits.width();
    }

    inline std::pair<int64_t, int64_t> pair(uint64_t i, uint64_t low1, 
        uint64_t low2, uint64_t ef_width) const {
        assert(i < size());  // and encode_prefix_sum);
        // uint64_t low1 = m_low_bits.access(i);
        // uint64_t low2 = m_low_bits.access(i + 1);
        // auto [low1, low2] = m_low_bits.get_two_consecutive_access(i);
        // uint64_t l = m_low_bits.width();
        uint64_t pos = m_high_bits_d1.select(m_high_bits, i);
        uint64_t h1 = pos - i;
        uint64_t h2 = bit_vector::unary_iterator(m_high_bits, pos + 1).next() - i - 1;
        uint64_t val1 = (h1 << ef_width) | low1;
        uint64_t val2 = (h2 << ef_width) | low2;
        return {val1, val2};
    }

    inline uint64_t diff(uint64_t i) const {
        auto [val1, val2] = pair(i);
        return val2 - val1;
    }

    // MODIFIED -> will not work, cv stored separately
    inline uint64_t size() const {
        return m_low_bits.size();
    }

    // MODIFIED -> removed low bit info
    uint64_t num_bits() const {
        return 8 * (m_high_bits.bytes() + m_high_bits_d1.bytes());
    }

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visitor.visit(m_high_bits);
        visitor.visit(m_high_bits_d1);
        visitor.visit(m_low_bits);
    }

    void save(std::ostream &os) {
        m_high_bits.save(os);
        m_high_bits_d1.save(os);
    }

    void load(std::istream &is){
        m_high_bits.load(is);
        m_high_bits_d1.load(is);
    }

    void save_fp(FILE* fp){
        m_high_bits.save_fp(fp);
        m_high_bits_d1.save_fp(fp);
        m_low_bits.save_fp(fp);
    }

    void load_fp(FILE* fp){
        m_high_bits.load_fp(fp);
        m_high_bits_d1.load_fp(fp);
        m_low_bits.load_fp(fp);
    }

private:
    bit_vector m_high_bits;
    darray1 m_high_bits_d1;
    compact_vector m_low_bits;
};

}  // namespace arank