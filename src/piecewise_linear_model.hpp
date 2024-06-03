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

#define DEBUG_PRINT 0

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
        // std::cout<<OA.dx <<" "<< OB.dy<<" "<<OA.dy <<" "<< OB.dx<<std::endl;
        // std::cout<<OA.dx * OB.dy<<" "<<OA.dy * OB.dx<<std::endl;
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
            // std::std::cout<<x<< " "<<last_x<<std::endl;
            throw std::logic_error("Points must be increasing by x.: "
                +std::to_string(uint64_t(x))+" "+std::to_string(uint64_t(last_x)));
        }
            

        // occurenceCount = (int64_t)1;
        // std::cout<<"curr point: "<<x<<std::endl;
        
        Point p1{x, SY(y) + epsilon + occurenceCount -1};
        Point p2{x, int64_t(y) - epsilon};
        // if(x >= 2188 && x<= 2334){
        //     std::cout<<"["<<(int64_t)SY(y) + epsilon + occurenceCount -1<<", ";
        //     std::cout<<(int64_t)SY(y) - epsilon<<"]"<<std::endl;
        // }
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
            points_in_hull = 0;
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
            return CanonicalSegment(rectangle[0], rectangle[1], first_x, last_x);
        return CanonicalSegment(rectangle, first_x, last_x);
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

    CanonicalSegment(const Point &p0, const Point &p1, X first, X last) : rectangle{p0, p1, p0, p1}, first(first), last(last) {};

    CanonicalSegment(const Point (&rectangle)[4], X first, X last)
        : rectangle{rectangle[0], rectangle[1], rectangle[2], rectangle[3]}, first(first), last(last) {};

    bool one_point() const {
        return rectangle[0].x == rectangle[2].x && rectangle[0].y == rectangle[2].y
            && rectangle[1].x == rectangle[3].x && rectangle[1].y == rectangle[3].y;
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

    std::tuple<int64_t, uint8_t, SY> get_fixed_point_segment(X origin, X max_input) const {
        if (one_point())
            return {0, 0, (rectangle[0].y + rectangle[1].y) / 2};
        // std::cout<<"in get_fixed_ps: "<<max_input<<std::endl;
        auto &p1 = rectangle[1];
        auto max_slope = rectangle[3] - rectangle[1];

        auto is_slope_integral = max_slope.dy % max_slope.dx == 0;
        auto slope_exponent = is_slope_integral ? 0 : (uint8_t) std::ceil(std::log2(max_input)) + 1;
        auto slope_significand = (max_slope.dy << slope_exponent) / max_slope.dx;

        auto intercept_n = max_slope.dy * (SX(origin) - p1.x);
        auto intercept_d = max_slope.dx;
        auto rounding_term = ((intercept_n < 0) ^ (intercept_d < 0) ? -1 : +1) * intercept_d / 2;
        auto intercept = (intercept_n + rounding_term) / intercept_d + p1.y;

        return {slope_significand, slope_exponent, intercept};
    }
    
    SY get_prediction(X kmer) const{
        auto &p1 = rectangle[1];
        auto max_slope = rectangle[3] - rectangle[1];
        auto origin = first;
        auto max_input = last-first+1;

        auto is_slope_integral = max_slope.dy % max_slope.dx == 0;
        auto slope_exponent = is_slope_integral ? 0 : (uint8_t) std::ceil(std::log2(max_input)) + 1;
        auto slope_significand = (max_slope.dy << slope_exponent) / max_slope.dx;

        auto intercept_n = max_slope.dy * (SX(origin) - p1.x);
        auto intercept_d = max_slope.dx;
        auto rounding_term = ((intercept_n < 0) ^ (intercept_d < 0) ? -1 : +1) * intercept_d / 2;
        auto intercept = (intercept_n + rounding_term) / intercept_d + p1.y;

        auto prediction = (((slope_significand) * (kmer - first)) >> slope_exponent) + intercept;
        return prediction;
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
            int64_t knot_si  = rectangle[1].y + eps;
            // int64_t knot_ei = rectangle[2].y;
            return {knot_si, knot_si};
        }
        X knot_start = first;
        auto max_slope = rectangle[3] - rectangle[1];
        int64_t knot_si = round(int64_t(rectangle[1].y) +
                    (max_slope.dy*(__int128_t(knot_start)-__int128_t(rectangle[1].x)))
                    /(double)max_slope.dx);
        int64_t knot_ei = round(int64_t(rectangle[1].y)+
                    (uint64_t(max_slope.dy))*((__int128_t(knot_end)-__int128_t(rectangle[1].x))
                    /(double)max_slope.dx));
        
        if(knot_start == 199799769495){
            std::cout<<"in get knot intersection: knot start: "<<knot_start
                <<" knot_end: "<<knot_end
                <<" knot si: "<<knot_si
                <<" knot ei: "<<knot_ei
                <<" slope: "<<double(max_slope.dy)/double(max_slope.dx)
                <<std::endl;
            std::cout<<" rectangle: \n";
            for(int64_t i=0; i<4; i++){
                std::cout<<int64_t(rectangle[i].x)<<" "<<int64_t(rectangle[i].y)<<"\n";
            }
        }

        return {knot_si, knot_ei};
    }

    bool isSlopeAtHalfPoint(X knot_end){
        if (one_point()){
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
size_t make_segmentation_basic_pla(int64_t n, int64_t epsilon, Fin in, Fout out, Ferr GetErr) {
    // std::cout<<"start of segment: "<<n<<" "<<epsilon<<std::endl;
    if (n == 0)
        return 0;

    using X = typename std::invoke_result_t<Fin, size_t>::first_type;
    using Y = typename std::invoke_result_t<Fin, size_t>::second_type;
    size_t c = 0;// count variable
    size_t start = 0;
    
    bool isIncluded = false;
    OptimalPiecewiseLinearModel<X, Y> opt(epsilon);
    opt.seg_start_indx = start;
    auto p = in(start);
    while(1){
        p = in(start); // pair of (x,y)
        if(p.first != -1){
            break;
        }
        start++;
    }
    // counting repeats to increase range
    int64_t curr_sa_indx = start, occurenceCount = 1;
    while(1){
        curr_sa_indx++;
        auto temp = in(curr_sa_indx);
        if(temp.first == p.first) occurenceCount++;
        else break;
    }
    
    start = curr_sa_indx; // as at curr indx new kmer has started        
    opt.add_point(p.first, p.second+epsilon, 1);
    
    // knot_sx = p.first;
    // std::cout<<p.first<<" "<<p.second<<" "<<std::endl;
    // int64_t unique_kmers = 1;
    
    // vector<int64_t> kval_vec;
    // vector<int64_t> range_start_vec;
    // vector<int64_t> range_end_vec;
    bool isForcedKnot = false, isSlopeIssuePresent = false;
    int64_t issueKmer = 0;
    bool isSlopeOk = false;
    while(!isSlopeOk){
        for (int64_t i = start; i < n; ++i) {
            if(isForcedKnot){
                isForcedKnot = false;
                out(opt.get_segment(), issueKmer);
                opt.reset();
                // std::cout<<"Forced index: "<<i<<std::endl;
                start = i-1;
                while(1){
                    p = in(start);
                    if(p.first != -1){
                        int64_t temp = p.first;
                        while (1)
                        {
                            start--;
                            auto t = in(start);
                            if(t.first != temp) break;
                        }
                        
                        i = start++;
                        break;
                    }
                    start--;
                }
                // std::cout<<"New i: "<<i<<std::endl;
                opt.seg_start_indx = i+1;
                
                ++c;
                continue;
            }
            auto next_p = in(i);
            
            if ((i != start && next_p.first == p.first) || next_p.first == -1)
                continue;
    //        std::cout<<i<<std::endl;
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
    //        std::cout<<occurenceCount<<std::endl;
            
            isIncluded = opt.add_point(p.first, p.second+epsilon, 1);
            if (!isIncluded) {
                /**
                 * TODO: Check to see if the seg slope is *.5.
                 * If it is, check whether any point in this segment has the error 
                 * value beyond epsilon
                 * If there is some point, put a forced knot at that point and 
                 * recalculate (?) starting from the first point at this segment
                 * */ 
                isSlopeOk = true;
                OptimalPiecewiseLinearModel<int64_t, uint64_t>::CanonicalSegment cs = opt.get_segment();
                // if(0){
                if(cs.get_first_x() == 199799769495){
                    auto [knot_si, knot_ei] = cs.get_knot_intersection(cs.get_last_x(), epsilon);
                    int64_t pred = cs.Pred_idx(knot_si, knot_ei, 215752959, 
                                cs.get_first_x(), cs.get_last_x(), epsilon);
                    std::cout<<"in add point: \n"
                        <<" knot start: "<<cs.get_first_x()
                        <<" knot_end: "<<cs.get_last_x()
                        <<" knot si: "<<knot_si
                        <<" knot ei: "<<knot_ei
                        <<" pred: "<<pred
                        <<" slope half point: "<<cs.isSlopeAtHalfPoint(cs.get_last_x())
                        <<std::endl;
                    uint64_t chk_idx= opt.seg_start_indx;
                    auto prev = in(chk_idx);
                    uint64_t prev_seed = prev.first;
                    for(chk_idx = opt.seg_start_indx+1; chk_idx<i-occurenceCount; chk_idx++){
                        auto point = in(chk_idx);
                        if(point.first == prev_seed || point.first == -1) continue;
                        int64_t pred = cs.Pred_idx(knot_si, knot_ei, point.first, 
                                cs.get_first_x(), cs.get_last_x(), epsilon);
                        
                        if(pred < 0) pred = 0;
                        
                        int32_t err = chk_idx - pred;;
                        // std::cout<<chk_idx<<" pred: "<<pred<<" err: "<<err<<std::endl;
                        if(err<0) err*= -1;
                        std::cout<<"query: "<<int64_t(point.first)
                            <<"chk idx "<<chk_idx<<
                            " pred: "<<pred
                            <<" err: "<<err<<std::endl;
                        
                        prev_seed = point.first;
                        if(prev_seed > 199799769501) break;
                    }
                }
                if(cs.isSlopeAtHalfPoint(cs.get_last_x())){
                    
                    uint64_t chk_idx= opt.seg_start_indx;
                    auto [knot_si, knot_ei] = cs.get_knot_intersection(cs.get_last_x(), epsilon);
                    // std::cout<<"SlopeHalfPoint, seg start idx: "<<chk_idx<<" not working idx: "<<i
                    //     <<" occur: "<<occurenceCount
                    //     <<" seg: "<<c
                    //     <<" knot_si: "<<knot_si<<" knot_ei: "<<knot_ei
                    //     <<std::endl;
                    /**
                     * curr_kval is not included in the curr segment, so looping only in the last segment
                     * first and last kmer/hash is guaranteed to be valid
                     * */ 
                    auto prev = in(chk_idx);
                    uint64_t prev_seed = prev.first;
                    for(chk_idx = opt.seg_start_indx+1; chk_idx<i-occurenceCount; chk_idx++){
                        auto point = in(chk_idx);
                        if(point.first == prev_seed || point.first == -1) continue;
                        int64_t pred = cs.Pred_idx(knot_si, knot_ei, point.first, 
                                cs.get_first_x(), cs.get_last_x(), epsilon);
                        
                        // std::cout<<"chk idx "<<chk_idx<<" pred: "<<pred<<std::endl;
                        
                        if(pred < 0) pred = 0;
                        
                        if(pred > chk_idx){
                            int32_t err = chk_idx - pred;
                            // std::cout<<chk_idx<<" pred: "<<pred<<" err: "<<err<<std::endl;
                            if(err<0) err*= -1;
                            if(err > epsilon){
                                isSlopeOk = false;
                                issueKmer = point.first;
                                break;
                            }
                        }
                        prev_seed = point.first;
                    }
                    // std::cout<<"isSlopeOk: "<<isSlopeOk<<" chkIdx: "<<chk_idx<<" kmer: "<<issueKmer<<std::endl;
                    if(!isSlopeOk){                        
                        isSlopeIssuePresent = true;
                        opt.reset();
                        i = opt.seg_start_indx - 1;
                        start = i+1;
                        // std::cout<<"slope issue: "<<c<<std::endl;
                    }
                }
                if(isSlopeOk){
                    out(opt.get_segment(), cs.get_last_x()); // this segment does not include curr_kval
                    opt.reset();
                    // current i points to the last occurence of the curr_kval
                    // need to point i to the (first occurence of the)
                    //last kmer at the last written segment
                    start = i - occurenceCount;
                    while(1){
                        p = in(start);
                        if(p.first != -1){
                            int64_t temp = p.first;
                            while (1)
                            {
                                start--;
                                //     if(start < 0) break; // should not happen
                                auto t = in(start);
                                if(t.first != temp) break;
                            }
                            
                            i = start++;
                            break;
                        }
                        start--;
                    }
                    opt.seg_start_indx = start;                    
                    ++c;
                }

                
            }
            else{
                if(isSlopeIssuePresent){
                    if(issueKmer == p.first){
                        isForcedKnot = true;
                        isSlopeIssuePresent = false;
                        // std::cout<<"slope issue found\n";
                    }
                }
            }
        }
        // have to check for the last segment whether slope is ok
        OptimalPiecewiseLinearModel<int64_t, uint64_t>::CanonicalSegment cs = opt.get_segment();
        isSlopeOk = true;
        if(cs.isSlopeAtHalfPoint(cs.get_last_x())){
            
            uint64_t chk_idx= opt.seg_start_indx;
            // std::cout<<"SlopeHalfPoint, seg start idx: "<<chk_idx<<" not working idx: "<<n
            //     <<" occur: "<<occurenceCount
            //     <<" seg: "<<c<<std::endl;
            auto [knot_si, knot_ei] = cs.get_knot_intersection(cs.get_last_x(), epsilon);
            /**
             * curr_kval is not included in the curr segment, so looping only in the last segment
             * first and last kmer/hash is guaranteed to be valid
             * */ 
            auto prev = in(chk_idx);
            uint64_t prev_seed = prev.first;
            for(chk_idx = opt.seg_start_indx+1; chk_idx<n-occurenceCount; chk_idx++){
                auto point = in(chk_idx);
                if(point.first == prev_seed || point.first == -1) continue;
                int64_t pred = cs.Pred_idx(knot_si, knot_ei, point.first, 
                        cs.get_first_x(), cs.get_last_x(), epsilon);
                if(pred < 0) pred = 0;            
                if(pred > chk_idx){
                    int32_t err = chk_idx - pred;;
                    if(err<0) err*= -1;
                    if(err > epsilon){
                        isSlopeOk = false;
                        issueKmer = point.first;
                        break;
                    }
                }
                prev_seed = point.first;
            }
            // std::cout<<"isSlopeOk: "<<isSlopeOk<<" chkIdx: "<<chk_idx<<" kmer: "<<issueKmer<<std::endl;
            if(!isSlopeOk){            
                isSlopeIssuePresent = true;
                opt.reset();
                start = opt.seg_start_indx;
                // std::cout<<"slope issue: "<<c<<std::endl;
            }
        }
        if(isSlopeOk){
            out(opt.get_segment(), cs.get_last_x()); // this segment does not include curr_kval
            opt.reset();
            ++c;
        }

    }
    // have to check whether there is slope issue for the last segment as well
    return c;
}


