// This file is part of PGM-index <https://github.com/gvinciguerra/PGM-index>.
// Copyright (c) 2018 Giorgio Vinciguerra.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#pragma once

#ifdef _OPENMP
#include <omp.h>
#else
// #warning Compilation with -fopenmp is recommended
typedef int omp_int_t;
inline omp_int_t omp_get_max_threads() { return 1; }
#endif

#include <cmath>
#include <limits>
#include <string>
#include <vector>
#include <stdexcept>
#include <type_traits>

#define DEBUG_PRINT_PLA 0
#define DEBUG_PRINT_FIRST 1


template<typename T>
using LargeSigned = typename std::conditional_t<std::is_floating_point_v<T>,
                                                long double,
                                                std::conditional_t<(sizeof(T) < 8), int64_t, __int128>>;

template<typename X, typename Y>
class OptimalPiecewiseLinearModel {
private:
    using SX = LargeSigned<X>;
    using SY = LargeSigned<Y>;

    struct Slope {
        SX dx{};
        SY dy{};

        bool operator<(const Slope &p) const {
            return dy * p.dx < dx * p.dy;
        }

        bool operator>(const Slope &p) const {
            return dy * p.dx > dx * p.dy;
        }

        bool operator==(const Slope &p) const {
            return dy * p.dx == dx * p.dy;
        }

        bool operator!=(const Slope &p) const {
            return dy * p.dx != dx * p.dy;
        }

        explicit operator long double() const {
            return dy / (long double) dx;
        }
    };

    struct StoredPoint {
        X x;
        Y y;
    };

    struct Point {
        X x{};
        SY y{};

        Slope operator-(const Point &p) const {
            return {SX(x) - p.x, y - p.y};
        }
    };

    template<bool Upper>
    struct Hull : private std::vector<StoredPoint> {
        const SY epsilon;
        
        explicit Hull(SY epsilon) : std::vector<StoredPoint>(), epsilon(Upper ? epsilon : -epsilon) {}

        Point operator[](size_t i) const {
            auto &p = std::vector<StoredPoint>::operator[](i);
            return {p.x, int64_t(p.y) + epsilon};
        }

        void clear() { std::vector<StoredPoint>::clear(); }
        void resize(size_t n) { std::vector<StoredPoint>::resize(n); }
        void reserve(size_t n) { std::vector<StoredPoint>::reserve(n); }
        size_t size() const { return std::vector<StoredPoint>::size(); }
        void push(X x, Y y) { std::vector<StoredPoint>::emplace_back(StoredPoint{x, y}); };
    };

    const Y epsilon;
    Hull<false> lower;
    Hull<true> upper;
    X first_x = 0;
    X last_x = 0;    
    size_t lower_start = 0;
    size_t upper_start = 0;
    size_t points_in_hull = 0;
    Point rectangle[4];

    

    auto cross(const Point &O, const Point &A, const Point &B) const {
        auto OA = A - O;
        auto OB = B - O;
        return (OA.dx * OB.dy) - (OA.dy * OB.dx);
    }

public:
    uint64_t seg_start_indx;
    class CanonicalSegment;
    /**
     * The explicit keyword in C++ is primarily used in constructors to indicate that the constructor 
     * should not be implicitly used for implicit type conversions. When you mark a constructor with 
     * explicit, you're telling the compiler not to perform automatic type conversions using that constructor.
    */
    explicit OptimalPiecewiseLinearModel(Y epsilon) : epsilon(epsilon), lower(epsilon), upper(epsilon) {
        if (epsilon < 0)
            throw std::invalid_argument("epsilon cannot be negative");

        upper.reserve(1u << 8);
        lower.reserve(1u << 8);
    }

    bool add_point(const X &x, const Y &y,  const Y &occurenceCount) {
        if (points_in_hull > 0 && x <= last_x){
            throw std::logic_error("Points must be increasing by x.: "
                +std::to_string(uint64_t(x))+" "+std::to_string(uint64_t(last_x)));
        }

        Point p1{x, SY(y) + epsilon + occurenceCount -1};
        Point p2{x, int64_t(y) - epsilon};
        
        if (points_in_hull == 0) {
            first_x = x;
            rectangle[0] = p1;
            rectangle[1] = p2;
            upper.clear();
            lower.clear();
            upper.push(x, y + occurenceCount -1);
            lower.push(x, y);
            upper_start = lower_start = 0;
            ++points_in_hull;
            last_x = x;
            return true;
        }

        if (points_in_hull == 1) {
            rectangle[2] = p2;
            rectangle[3] = p1;
            upper.push(x, y+ occurenceCount -1);
            lower.push(x, y);
            ++points_in_hull;
            last_x = x;
            return true;
        }

        if (epsilon == 0) {
            auto p1_on_line1 = p1 - rectangle[0] == rectangle[2] - rectangle[0];
            points_in_hull = p1_on_line1 ? points_in_hull + 1 : 0;
            last_x = x;
            return p1_on_line1;
        }

        auto slope1 = rectangle[2] - rectangle[0];
        auto slope2 = rectangle[3] - rectangle[1];
        bool outside_line1 = p1 - rectangle[2] < slope1; // minus has higher precedence
        bool outside_line2 = p2 - rectangle[3] > slope2;

        if (outside_line1 || outside_line2) {
            // points_in_hull = 0;
            return false;
        }
        last_x = x;
        if (p1 - rectangle[1] < slope2) {
            // Find extreme slope
            auto min = lower[lower_start] - p1;
            auto min_i = lower_start;
            for (auto i = lower_start + 1; i < lower.size(); i++) {
                auto val = (lower[i] - p1);
                if (val > min)
                    break;
                min = val;
                min_i = i;
            }

            rectangle[1] = lower[min_i];
            rectangle[3] = p1;
            lower_start = min_i;

            // Hull update
            auto end = upper.size();
            for (; end >= upper_start + 2 && cross(upper[end - 2], upper[end - 1], p1) <= 0; --end);
            upper.resize(end);
            upper.push(x, y+ occurenceCount -1);
        }

        if (p2 - rectangle[0] > slope1) {
            // Find extreme slope
            auto max = upper[upper_start] - p2;
            auto max_i = upper_start;
            for (auto i = upper_start + 1; i < upper.size(); i++) {
                auto val = (upper[i] - p2);
                if (val < max)
                    break;
                max = val;
                max_i = i;
            }

            rectangle[0] = upper[max_i];
            rectangle[2] = p2;
            upper_start = max_i;

            // Hull update
            auto end = lower.size();
            for (; end >= lower_start + 2 && cross(lower[end - 2], lower[end - 1], p2) >= 0; --end);
            lower.resize(end);
            lower.push(x, y);
        }

        ++points_in_hull;

        return true;
    }

    CanonicalSegment get_segment() const {
        if (points_in_hull == 1)
            return CanonicalSegment(rectangle[0], rectangle[1], first_x, last_x, points_in_hull);
        return CanonicalSegment(rectangle, first_x, last_x, points_in_hull);
    }

    void reset() {
        points_in_hull = 0;
        lower.clear();
        upper.clear();
    }
};

template<typename X, typename Y>
class OptimalPiecewiseLinearModel<X, Y>::CanonicalSegment {
    friend class OptimalPiecewiseLinearModel;

    Point rectangle[4];
    X first, last;
    int64_t points_in_hull;

    CanonicalSegment(const Point &p0, const Point &p1, X first, X last, int64_t pih) 
        : rectangle{p0, p1, p0, p1}, first(first), last(last), points_in_hull(pih) {};

    CanonicalSegment(const Point (&rectangle)[4], X first, X last, int64_t pih)
        : rectangle{rectangle[0], rectangle[1], rectangle[2], rectangle[3]}, 
            first(first), last(last), points_in_hull(pih)  {};

    bool one_point() const {
        return rectangle[0].x == rectangle[2].x && rectangle[0].y == rectangle[2].y
            && rectangle[1].x == rectangle[3].x && rectangle[1].y == rectangle[3].y;
    }

    bool two_point() const{
        return points_in_hull == 2;
    }

public:

    CanonicalSegment() = default;

    explicit CanonicalSegment(X first) : CanonicalSegment({first, 0}, {first, 0}, first) {};

    X get_first_x() const { return first; }
    X get_last_x() const { return last; }

    CanonicalSegment copy(X x) const {
        auto c(*this);
        c.first = x;
        return c;
    }

   
    std::string ConvertToString(__int128 num) {
        std::string str;
        do {
            int digit = num % 10;
            str = std::to_string(digit) + str;
            num = (num - digit) / 10;
        } while (num != 0);
        return str;
    }

    
    std::tuple<int64_t, int64_t> get_knot_intersection(X knot_end, int64_t eps){
        if (one_point()){
            int64_t knot_si  = rectangle[1].y;// TODO: 2*eps
            
            return {knot_si, knot_si};
        }
        X knot_start = first;
        // TODO: Can adjust the values so that at query, mid is the first match
        // slope does not matter, but knot_si has to be increasing
        if(two_point()){
            auto p1 = Point{rectangle[0].x,  rectangle[0].y - 2*eps};
            auto p2 = Point{rectangle[2].x,  rectangle[2].y + eps};
            auto slope = p2-p1;
            
            int64_t knot_si = rectangle[0].y - 2*eps;
            int64_t knot_ei = round(int64_t(knot_si)+
                    (uint64_t(slope.dy))*((__int128_t(knot_end)-__int128_t(rectangle[1].x))
                    /(double)slope.dx));
            
            return {knot_si, knot_ei};
        }

        
        auto max_slope = rectangle[3] - rectangle[1];
        
        
        int64_t knot_si = round(int64_t(rectangle[1].y) +
                    (max_slope.dy*(__int128_t(knot_start)-__int128_t(rectangle[1].x)))
                    /(double)max_slope.dx);
        int64_t knot_ei = round(int64_t(rectangle[1].y)+
                    (uint64_t(max_slope.dy))*((__int128_t(knot_end)-__int128_t(rectangle[1].x))
                    /(double)max_slope.dx));
        

        return {knot_si, knot_ei};
    }

    bool isSlopeAtHalfPoint(X knot_end, int64_t eps){
        if (one_point() || two_point()){
            return false;
        }                
        X knot_start = first;
        auto max_slope = rectangle[3] - rectangle[1];
        
        double knot_si = int64_t(rectangle[1].y) +
                    (max_slope.dy*(__int128_t(knot_start)-__int128_t(rectangle[1].x)))
                    /(double)max_slope.dx;
        double knot_ei = int64_t(rectangle[1].y)+
                    (uint64_t(max_slope.dy))*((__int128_t(knot_end)-__int128_t(rectangle[1].x))
                    /(double)max_slope.dx);
        return ((floor(knot_si-0.5)== floor(knot_si))&& (floor(knot_ei - 0.5)==floor(knot_ei)));
        
    }

    int64_t Pred_idx(int64_t brk_beg_sa_indx, int64_t brk_end_sa_indx,
                uint64_t query_val, uint64_t brk_beg_kval, 
                uint64_t brk_end_kval, int64_t eps){
        return round(brk_beg_sa_indx + 
                    (brk_end_sa_indx - brk_beg_sa_indx)*
                    ((uint64_t(query_val) - uint64_t(brk_beg_kval)) /
                    (double)(brk_end_kval - brk_beg_kval)) ) - eps;
    }
    
};



template<typename Fin, typename Fout, typename Ferr>
size_t make_segmentation_rep_pla(int64_t n, int64_t epsilon, Fin in, Fout out, Ferr GetErr) {
    if (n == 0)
        return 0;

    using X = typename std::invoke_result_t<Fin, size_t>::first_type;
    using Y = typename std::invoke_result_t<Fin, size_t>::second_type;
    size_t c = 0;// count variable
    size_t start = 0;
    
    
    bool isIncluded = false;
    OptimalPiecewiseLinearModel<X, Y> opt(epsilon);
    
    auto p = in(start);
    opt.seg_start_indx = start;
    // counting repeats to increase range
    int64_t curr_sa_indx = start, occurenceCount = 1;
    while(1){
        curr_sa_indx++;
        auto temp = in(curr_sa_indx);
        if(temp.first == p.first) occurenceCount++;
        else break;
    }
    
    start = curr_sa_indx; // as at curr indx new kmer has started        
    opt.add_point(p.first, p.second+epsilon, occurenceCount);
    
    bool isForcedKnot = false, isSlopeIssuePresent = false;
    uint64_t issueKmer = 0;
    bool isSlopeOk = false;
    while(!isSlopeOk){
        for (int64_t i = start; i < n; ++i) {
           
            auto next_p = in(i);
            
            if ((i != start && next_p.first == p.first))
                continue;
            p = next_p;
            curr_sa_indx = i;
            occurenceCount = 1;
            while(1){
                curr_sa_indx++;
                auto temp = in(curr_sa_indx);
                if(temp.first == p.first) occurenceCount++;
                else break;
            }
            i+=(occurenceCount-1);
            if(isSlopeIssuePresent && p.first == issueKmer){
                isSlopeIssuePresent = false;
                isSlopeOk = true;
                OptimalPiecewiseLinearModel<uint64_t, uint64_t>::CanonicalSegment cs = opt.get_segment();
                
                if(cs.isSlopeAtHalfPoint(issueKmer, epsilon)){
                    uint64_t chk_idx= opt.seg_start_indx;
                    
                    auto [knot_si, knot_ei] = cs.get_knot_intersection(issueKmer, epsilon);
                    
                    auto prev = in(chk_idx);
                    uint64_t prev_seed = prev.first;
                    for(chk_idx = opt.seg_start_indx+1; chk_idx<i-occurenceCount+1; chk_idx++){
                        auto point = in(chk_idx);
                        if(point.first == prev_seed) continue;
                        int64_t pred;
                        pred = cs.Pred_idx(knot_si, knot_ei, point.first, 
                                cs.get_first_x(), issueKmer, epsilon);
                        
                        if(point.first == DEBUG_PRINT_PLA){
                            int64_t err = GetErr(chk_idx, pred);
                        }
                        
                        if(pred < 0) pred = 0;
                        if(pred > chk_idx){
                            int64_t err = GetErr(chk_idx, pred);
                            if(err<0) err*= -1;
                            if(err > epsilon){
                                isSlopeOk = false;  
                                issueKmer = point.first;
                                break;
                            }
                        }
                        prev_seed = point.first;
                    }
                    if(!isSlopeOk){                        
                        isSlopeIssuePresent = true;
                        opt.reset();
                        i = opt.seg_start_indx - 1;
                        start = i+1;
                        continue;
                    }
                }
                out(opt.get_segment(), issueKmer);
                ++c;
                opt.reset();
                opt.seg_start_indx = i-occurenceCount+1;
            }
            isIncluded = opt.add_point(p.first, p.second+epsilon, occurenceCount);
            if (!isIncluded) {
                isSlopeOk = true;
                OptimalPiecewiseLinearModel<uint64_t, uint64_t>::CanonicalSegment cs = opt.get_segment();
                
                if(cs.isSlopeAtHalfPoint(p.first, epsilon)){
                
                    uint64_t chk_idx= opt.seg_start_indx;
                    
                    auto [knot_si, knot_ei] = cs.get_knot_intersection(p.first, epsilon);
                    
                    auto prev = in(chk_idx);
                    
                    uint64_t prev_seed = prev.first;
                    for(chk_idx = opt.seg_start_indx+1; chk_idx<i-occurenceCount+1; chk_idx++){
                        auto point = in(chk_idx);
                        if(point.first == prev_seed) continue;
                        int64_t pred;
                        pred = cs.Pred_idx(knot_si, knot_ei, point.first, 
                                cs.get_first_x(), p.first, epsilon);
                        
                        if(pred < 0) pred = 0;
                        
                        if(pred > chk_idx){
                            int64_t err = GetErr(chk_idx, pred);
                            if(err<0) err*= -1;
                            if(err > epsilon){
                                isSlopeOk = false;
                                issueKmer = point.first;
                                break;
                            }
                        }
                        prev_seed = point.first;
                    }
                     
                    if(!isSlopeOk){                        
                        isSlopeIssuePresent = true;
                        opt.reset();
                        i = opt.seg_start_indx - 1;
                        start = i+1;
                    }
                }
                if(isSlopeOk){
                    out(opt.get_segment(), p.first); 
                    opt.reset();
                    start = i - occurenceCount + 1;
                    i = start - 1;
                    opt.seg_start_indx = start;                    
                    ++c;
                }
            }
        }
        // have to check for the last segment whether slope is ok
        OptimalPiecewiseLinearModel<uint64_t, uint64_t>::CanonicalSegment cs = opt.get_segment();
        isSlopeOk = true;
        if(cs.isSlopeAtHalfPoint(p.first, epsilon)){
            uint64_t chk_idx= opt.seg_start_indx;
           
            auto [knot_si, knot_ei] = cs.get_knot_intersection(p.first, epsilon);
            
            auto prev = in(chk_idx);
            uint64_t prev_seed = prev.first;
            for(chk_idx = opt.seg_start_indx+1; chk_idx<n-occurenceCount; chk_idx++){
                auto point = in(chk_idx);
                if(point.first == prev_seed) continue;
                int64_t pred;
                pred = cs.Pred_idx(knot_si, knot_ei, point.first, 
                        cs.get_first_x(), p.first, epsilon);


                if(pred < 0) pred = 0;            
                if(pred > chk_idx){
                    int64_t err = GetErr(chk_idx, pred);
                    if(err<0) err*= -1;
                    if(err > epsilon){
                        isSlopeOk = false;
                        issueKmer = point.first;
                        break;
                    }
                }
                prev_seed = point.first;
            }
            
            if(!isSlopeOk){            
                isSlopeIssuePresent = true;
                opt.reset();
                start = opt.seg_start_indx;
            
            }
        }
        if(isSlopeOk){
            out(opt.get_segment(), p.first); // this segment does not include curr_kval
            opt.reset();
            ++c;
        }

    }
    // have to check whether there is slope issue for the last segment as well
    return c;
}
