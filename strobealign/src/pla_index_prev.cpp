#include "pla_index.hpp"

static Logger& logger = Logger::get();

// constructor at index build
pla_index::pla_index(int64_t eps, int64_t lp_bits):
            epsilon(eps),
            lp_bits(lp_bits),
            brk_uniform_diff_packed(true)
{
    knot_bs_thres = 64;
}

// constructor at index query
pla_index::pla_index(int64_t knot_bs_thres):
    knot_bs_thres(knot_bs_thres), 
    brk_uniform_diff_packed(true)
{}

void pla_index::insert_to_ind_vec(uint64_t brk_beg_kval, int64_t brk_indx){
    int64_t lp_bit_value = brk_beg_kval >> shift_bits;
    if(indirection_vec.size() > lp_bit_value) return;
    for(int64_t iv=indirection_vec.size(); iv<lp_bit_value+1; iv++){
        indirection_vec.emplace_back(brk_indx);
    }
}

void pla_index::encode_knots(const vector<uint64_t> &brk_kval_vec){
    vector<int64_t> brk_uniform_diff_vec;
    uint64_t max_bits = shift_bits + lp_bits;
    if(max_bits == 64) max_data_ratio = ~UINT64_C(0)/(dic_size - 1);
    else max_data_ratio = ((1ULL << (max_bits)) - 1)/(dic_size - 1);
    logger.debug()<<"Max Data Ratio: "<<max_data_ratio<<endl;
    // ofstream seed_knot("knot_diff.txt");
    // seed_knot<<max_data_ratio<<endl;

    for(uint64_t i=0; i<dic_size; i++){
        int64_t diff = __int128_t(brk_kval_vec[i]) - __int128_t(i*max_data_ratio);
        // seed_knot<<brk_kval_vec[i] << " "<<diff<<" "<<(static_cast<uint64_t>(diff) + max_data_ratio*i)<<endl;
        brk_uniform_diff_vec.emplace_back(diff);
    }
    
    brk_uniform_diff_packed.BuilidSignedPackedVector(brk_uniform_diff_vec);
    logger.debug()<<"Max Data Ratio: "<<brk_uniform_diff_vec[0]<<" "<<
        brk_uniform_diff_vec[dic_size-1]<<endl;
    logger.debug()<<"bits needed: "<<brk_uniform_diff_packed.GetPacketSize()<<endl;
}

void pla_index::SetBit(uint64_t &number, uint64_t pos){
    uint64_t mask = 1ULL << (64-pos); // pos = [1, 64]
    number |= mask;
}


uint64_t pla_index::BinarySearch_BrkPnt(uint64_t lo, uint64_t hi, const uint64_t query) const{
    uint64_t mid, midVal;
    // cout<<query<<endl;
    while(1){
        // if(hi - lo <= 1){
        //     if(dic_vals[hi] <= query) return hi;
        //     else return lo;
        // }
        // cout<<lo<<" "<<mid<<" "<<hi<<endl;
        // cout<<dic_vals[lo]<<" "<<dic_vals[mid]<<" "<<dic_vals[hi]<<endl;
        if(hi-lo <= this->knot_bs_thres){
            for (uint64_t i = hi; i >= lo; i--)
            {
                if(brk_uniform_diff_packed.GetValueAt(i)+max_data_ratio*i  <= query) return i;
                // if(knot_packd.GetValueAt(i) <= query) return i;
            }
            logger.debug()<<"Problem in topk base case\n";
            exit(0);
        }
        mid = (lo + hi) >> 1;
        // cout<<lo<<" "<<mid<<" "<<hi<<" "<<dic_vals[hi]<<" "<<query<<endl;
        midVal = brk_uniform_diff_packed.GetValueAt(mid)+max_data_ratio*mid;
        // midVal = knot_packd.GetValueAt(mid);
        if (midVal == query) return mid;
        else if(midVal > query){ //search left half
            hi = mid;
        }
        else{ // search right half
            lo = mid;
        }
    }
}

int64_t pla_index::binary_search(int64_t lo, int64_t hi, const uint64_t kval) const{
    uint64_t mid, midVal;
    while(1){
        if(hi-lo <= knot_bs_thres){
            for (uint64_t i = hi; i >= lo; i--)
            {
                if(static_cast<uint64_t>(brk_uniform_diff_packed.GetValueAt(i))+max_data_ratio*i <= kval) return i;
            }
            logger.debug()<<"Problem in topk base case\n"
                <<"kval: "<<kval
                <<endl;
            exit(0);
        } 
        mid = (lo + hi) >> 1;
        midVal = static_cast<uint64_t>(brk_uniform_diff_packed.GetValueAt(mid))+max_data_ratio*mid;
        if (midVal == kval) return mid;
        else if(midVal > kval){ 
            hi = mid;
        }
        else{ 
            lo = mid;
        }
    }  
}

uint64_t pla_index::QueryKmer(const vector<RefRandstrobe> &randstrobes, 
            const uint64_t query_val) const{
    // logger.debug()<<"\n\n";
    // ofstream all_query("query_new_all.txt", ios_base::app);
    int64_t lp_bit_value = query_val >> shift_bits;
    int64_t brkpnt, query_pred, direction;
    int64_t ind_brk_beg_idx = indirection_vec[lp_bit_value];
    // all_query<<endl<<query_val<<" "<<ind_brk_beg_idx<<" "
    //     << brk_uniform_diff_packed.GetValueAt(ind_brk_beg_idx)+
    //         max_data_ratio*ind_brk_beg_idx<<" "
    //     << lp_bit_value<<" ";
    bool isOnIntersection = false;
    // logger.debug()<<"q: "<<query_val<<" lp: "<< lp_bit_value<<" brk[ind_beg]: "<<brk_kval_vec[ind_brk_beg_idx] <<endl;
    if(brk_uniform_diff_packed.GetValueAt(ind_brk_beg_idx)+
            max_data_ratio*ind_brk_beg_idx > query_val) direction = 0;
    else if(brk_uniform_diff_packed.GetValueAt(ind_brk_beg_idx)+
        max_data_ratio*ind_brk_beg_idx < query_val) direction = 1;
    else{
        isOnIntersection = true;
        query_pred = y_range_beg.access(ind_brk_beg_idx) - epsilon;
    }
    
    if(!isOnIntersection){
        // logger.debug()<<"dir: "<<direction<<endl;
        if(direction == 1){
            if(lp_bit_value == indirection_vec.size()-1){
                brkpnt =  binary_search(ind_brk_beg_idx, dic_size-1, query_val);    
                // all_query<<dic_size-1<<" ";
            }
            else{
                brkpnt =  binary_search(ind_brk_beg_idx, indirection_vec[lp_bit_value+1]-1, query_val);
                // all_query<<indirection_vec[lp_bit_value+1]-1<<" ";
            }
        }
        else{
            lp_bit_value--;
            if(lp_bit_value<0) return -1;
            int64_t ind_brk_end_idx;
            while(1){
                ind_brk_end_idx = indirection_vec[lp_bit_value];
                if(ind_brk_end_idx != ind_brk_beg_idx) break;
                lp_bit_value--;
            }
            brkpnt = binary_search(ind_brk_end_idx, ind_brk_beg_idx, query_val);
            // all_query<<ind_brk_end_idx<<" ";
        }
        // all_query<<brkpnt<<" "<<brk_uniform_diff_packed.GetValueAt(brkpnt)+max_data_ratio*brkpnt<<" ";
        // logger.debug()<<"brk: "<<brkpnt<<" val at brk: "<<brk_kval_vec[brkpnt]<<endl;
        if(brk_uniform_diff_packed.GetValueAt(brkpnt)+max_data_ratio*brkpnt  == query_val) {
            query_pred = y_range_beg.access(brkpnt) - epsilon;
        }
        else{
            uint64_t brk_beg_kval, brk_end_kval;
            int64_t brk_beg_sa_indx, brk_end_sa_indx;
            brk_beg_kval = brk_uniform_diff_packed.GetValueAt(brkpnt)+max_data_ratio*brkpnt;
            brk_end_kval = brk_uniform_diff_packed.GetValueAt(brkpnt+1)+max_data_ratio*(brkpnt+1);
            brk_beg_sa_indx = y_range_beg.access(brkpnt) - epsilon;
            int64_t diff_val = sa_diff_dac_vec[brkpnt];
            if(diff_val & 1){
                brk_end_sa_indx = y_range_beg.access(brkpnt+1) - (diff_val>>1)
                                     - epsilon;
            }
            else{
                brk_end_sa_indx = y_range_beg.access(brkpnt+1) + (diff_val>>1)
                                     - epsilon;
            }
            // all_query<<brk_beg_kval<<" "
            //     <<brk_end_kval<<" "
            //     <<brk_beg_sa_indx<<" "
            //     <<brk_end_sa_indx<<" ";
            // logger.debug()<<"Interval: \n"<<brk_beg_kval<<"\n"<<query_val<<"\n"<<brk_end_kval<<endl;
            // logger.debug()<<query_val<<endl;
            // logger.debug()<<brk_beg_kval<<" "
            //     <<brk_end_kval<<" "
            //     <<brk_beg_sa_indx<<" "
            //     <<brk_end_sa_indx<<endl;
            // query_pred = round(brk_beg_sa_indx + 
            //         ((double)(query_val - brk_beg_kval)/(brk_end_kval - brk_beg_kval) )*
            //         (brk_end_sa_indx - brk_beg_sa_indx));
            
            query_pred = round(brk_beg_sa_indx + 
                (brk_end_sa_indx - brk_beg_sa_indx)*
                ((uint64_t(query_val) - uint64_t(brk_beg_kval)) /
                (double)(brk_end_kval - brk_beg_kval)) );
            // logger.debug()<<"Here, Pred: "<<query_pred<<endl;
            // logger.debug()<<"Pred: "<<query_pred<<" "<<GetKmerAtStrPos(sa[query_pred])<<endl;
        }  
        // logger.debug()<<"Second, Pred: "<<query_pred<<endl;          
    }
    // logger.debug()<<"After bracket Pred: "<<query_pred<<endl;
    if(query_pred < int64_t(0)) query_pred = int64_t(0);
    // logger.debug()<<"Pred: "<<query_pred<<" "<<GetKmerAtStrPos(sa[query_pred])<<endl;
    if(query_pred > int64_t(total_data_points)-1) query_pred = int64_t(total_data_points)-1;
    // logger.debug()<<"after mod Pred: "<<query_pred<<" val at pred: "<<randstrobes[query_pred].hash<<endl;
    // all_query<<query_pred<<endl;
    uint64_t hash = randstrobes[query_pred].hash;
    if(hash == query_val) {
        while(1){
            if(!query_pred) return 0;
            if(randstrobes[--query_pred].hash != query_val) return query_pred+1;
        }
        // return GetFirst(randstrobes, query_val, query_pred); // first index is guaranteed        
    }
    // Paul's suggestion was to do seq search starting from the prediciton directly
    if(hash < query_val){
        int64_t last_idx = min(int64_t(query_pred + epsilon), int64_t(total_data_points - 1));
        for(int64_t idx = query_pred+1; idx<=last_idx; idx++){
            if(randstrobes[idx].hash == query_val) return idx; // guranteed to be the first index
        }
        return -1;
        // return BS_Seed(randstrobes, query_val, query_pred, min(uint64_t(query_pred + epsilon), fSize - 1));
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
        // return BS_Seed(randstrobes, query_val, max(query_pred - epsilon-1, int64_t(0)), query_pred);
    }
}


    // strobealign version
int64_t pla_index::query(const vector<RefRandstrobe> &randstrobes, const uint64_t query_val) const{
    uint64_t lp_bit_value = query_val >> shift_bits;
    int64_t brkpnt, query_pred;
    int64_t ind_brk_beg_idx = indirection_vec[lp_bit_value];
    uint64_t ind_point_val = static_cast<uint64_t>(brk_uniform_diff_packed.GetValueAt(ind_brk_beg_idx))+max_data_ratio*ind_brk_beg_idx;   
    // logger.debug()<<knot_bs_thres<<endl;
    if(lp_bit_value == 0 && ind_point_val > query_val) return -1;    
    else if(ind_point_val == query_val){
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

            // lp_bit_value--;
            // int64_t ind_brk_end_idx;
            // while(1){
            //     ind_brk_end_idx = indirection_vec[lp_bit_value];
            //     if(ind_brk_end_idx != ind_brk_beg_idx) break;
            //     lp_bit_value--;
            // }
            brkpnt = binary_search(ind_point_val, ind_brk_beg_idx, query_val);
        }
        // cout<<query_val<<" "<<brkpnt<<" "<<direction<<endl;
        const uint64_t brk_beg_kval = static_cast<uint64_t>(brk_uniform_diff_packed.GetValueAt(brkpnt))+max_data_ratio*brkpnt;
        if(brk_beg_kval == query_val) {
            query_pred = y_range_beg.access(brkpnt) - epsilon;
        }
        else{
            // int64_t brk_beg_kval, brk_end_kval, brk_beg_sa_indx, brk_end_sa_indx;
            
            const uint64_t brk_end_kval = static_cast<uint64_t>(brk_uniform_diff_packed.GetValueAt(brkpnt+1))+max_data_ratio*(brkpnt+1);
            const int64_t brk_beg_sa_indx = int64_t(y_range_beg.access(brkpnt)) - epsilon;
            
            const int64_t diff_val = sa_diff_dac_vec[brkpnt];
            const int64_t brk_end_sa_indx = (diff_val & 1)
                                    ? y_range_beg.access(brkpnt + 1) - (diff_val >> 1) - epsilon
                                    : y_range_beg.access(brkpnt + 1) + (diff_val >> 1) - epsilon;
            // logger.debug()<<"Interval: \n"<<brk_beg_kval<<"\n"<<query_val<<"\n"<<brk_end_kval<<endl;
            // logger.debug()<<brkpnt<<" "<<max_data_ratio<<endl;
            // logger.debug()<<brk_beg_kval<<" "
            //     <<brk_end_kval<<" "
            //     <<brk_beg_sa_indx<<" "
            //     <<brk_end_sa_indx<<endl;
            // query_pred = round(brk_beg_sa_indx + 
            //         ((double)(query_val - brk_beg_kval)/(brk_end_kval - brk_beg_kval) )*
            //         (brk_end_sa_indx - brk_beg_sa_indx));
            query_pred = round(brk_beg_sa_indx + 
                (brk_end_sa_indx - brk_beg_sa_indx)*
                (((query_val) - (brk_beg_kval)) /
                (double)(brk_end_kval - brk_beg_kval)) );
            // logger.debug()<<"Here, Pred: "<<query_pred<<endl;
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

    // ofstream query_file("query_new.txt", ios_base::app);
    // query_file<<query_val<<" "<<query_pred<<endl;

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
        // return BS_Seed(randstrobes, query_val, query_pred, min(uint64_t(query_pred + epsilon), fSize - 1));
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

void pla_index::build_rep_stretch_pla_index_strobe(const vector<RefRandstrobe> &randstrobes){
    vector<uint64_t>  brk_sa_indx_vec, brk_sa_end_vec, brk_kval_vec;
    vector<int64_t> brk_diff_vec;
    bool isFirst = true;
    int64_t prev_knot_ei=0, diff_from_curr_start;
    vector<canonical_segment> out(1);

    uint64_t bits;
    uint64_t largest = randstrobes[randstrobes.size()-1].hash;
    if(ceil(log2(largest)) == floor(log2(largest))) bits = ceil(log2(largest)) + 1;
    else bits = ceil(log2(largest));
    if(bits < lp_bits){
        throw std::logic_error("lp cannot be greater than log2(largest)");
    }
    shift_bits = bits - lp_bits;
    total_data_points = randstrobes.size();

    auto in_fun = [&](auto i) { return std::pair<uint64_t, uint64_t>(randstrobes[i].hash, i); };
        
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
        return get_pred_diff(act_idx, pred_idx, randstrobes);
    };

    uint64_t seg_count = make_segmentation_rep_pla_strobe(total_data_points, epsilon, in_fun, out_fun, get_err);

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

    logger.debug()<<"dic size: "<<dic_size<<endl;

    sdsl::dac_vector_dp<> temp(brk_diff_vec);
    sa_diff_dac_vec = temp;
    y_range_beg.encode(brk_sa_indx_vec.data(), brk_sa_indx_vec.size());
    encode_knots(brk_kval_vec);
}


const vector<uint64_t> pla_index::get_processed_ind_vec() const{
    vector<uint64_t> ind_diff_vec;
    int64_t _count = 0;
    for(int64_t i=1; i<indirection_vec.size(); i++){
        uint64_t diff = indirection_vec[i] - indirection_vec[i-1];
        if(diff == 0) _count++;
        ind_diff_vec.push_back(diff);
    }
    return ind_diff_vec;
}

void pla_index::Load(std::istream& is, const vector<RefRandstrobe> &randstrobes){
    total_data_points = randstrobes.size();
    uint64_t largest = randstrobes[total_data_points-1].hash;
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
    if (total_bits == 64) max_data_ratio = ~UINT64_C(0);
    else max_data_ratio = ((1ULL << total_bits) - 1)/(dic_size-1);
    
    sa_diff_dac_vec.load(is);
    
    // cout<<"ef load\n";
    y_range_beg.load(is);
    // cout<<"max err load\n";
    uint16_t max_error_t;
    essentials_arank::load_pod(is, max_error_t);
    // is.read(reinterpret_cast<char*>(&max_error_t), sizeof(max_error_t));
    // err = fread(&max_error_t, sizeof(uint16_t), 1, fp);
    epsilon = (int64_t)max_error_t;    
}

void pla_index::save_unpacked(string dic_fn) const{
    dic_fn = dic_fn + ".unpacked";
    std::ofstream os(dic_fn, std::ios::binary);
    uint8_t t_lp_bits = lp_bits;
    essentials_arank::save_pod(os, t_lp_bits);
    os.write(reinterpret_cast<char const*>(indirection_vec.data()),
                  indirection_vec.size() * sizeof(uint64_t));
    uint32_t nBrkpnts = dic_size;
    essentials_arank::save_pod(os, nBrkpnts);
        
    sa_diff_dac_vec.serialize(os);
    
    y_range_beg.save(os);
    uint16_t t_err_th = epsilon;
    essentials_arank::save_pod(os, t_err_th);    
}

void pla_index::save_bit_info(string dic_fn, CBitPacking& ind_diff_pack) const{
    int pos = dic_fn.find(".bin");
    string pref = dic_fn.substr(0, pos);
    string suff = dic_fn.substr(pos);
    string bit_fn = pref+"_bitInfo"+suff;
    ofstream bit_info(bit_fn.c_str());
    bit_info<<"#Breakpoints: "<<dic_size<<endl;
    bit_info<<"Indirection vector: each: "<<uint64_t(ind_diff_pack.GetPacketSize())<<endl;
    bit_info<<"knot kmers: each: "<<uint64_t(brk_uniform_diff_packed.GetPacketSize())<<endl;
    bit_info<<"Y beg each: "<<ceil(y_range_beg.num_bits()/dic_size)<<endl;    
    bit_info<<"Y end diff: each: "<<
        ceil((size_in_bytes(sa_diff_dac_vec)*8.0)/sa_diff_dac_vec.size())<<endl;
    bit_info.close();
}

void pla_index::Save(std::ostream& os) const{
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


    sa_diff_dac_vec.serialize(os);

    y_range_beg.save(os);
    uint16_t t_err_th = epsilon;
    essentials_arank::save_pod(os, t_err_th);
    // os.write(reinterpret_cast<char const*>(&t_err_th),
    //               sizeof(t_err_th));
    string dic_fn = "stroalign_meta.bin";
    save_unpacked(dic_fn);
    save_bit_info(dic_fn, ind_diff_pack);
}