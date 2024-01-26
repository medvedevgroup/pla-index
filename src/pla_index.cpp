#include "pla_index.hpp"

// constructor at index build
pla_index::pla_index(int64_t eps, int64_t lp_bits, uint64_t largest,
            bool isFirstIndxRet, uint64_t total_points, INDX_TYPE it):
            epsilon(eps),
            lp_bits(lp_bits), isFirstIndexReturned(isFirstIndxRet),
            indx_type(it),
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

// constructor at index query
pla_index::pla_index(int64_t knot_bs_thres, uint64_t largest, 
    uint64_t total_points, string dic_fn, INDX_TYPE it):
    knot_bs_thres(knot_bs_thres), 
    total_data_points(total_points),
    indx_type(it),
    brk_uniform_diff_packed(true), diff_packd(true)
{
    uint64_t bits;
    if(ceil(log2(largest)) == floor(log2(largest))) bits = ceil(log2(largest)) + 1;
    else bits = ceil(log2(largest));
    this->Load(dic_fn, bits);
}

void pla_index::insert_to_ind_vec(int64_t brk_beg_kval, int64_t brk_indx){
    int64_t lp_bit_value = brk_beg_kval >> shift_bits;
    if(indirection_vec.size() > lp_bit_value) return;
    for(int64_t iv=indirection_vec.size(); iv<lp_bit_value+1; iv++){
        indirection_vec.emplace_back(brk_indx);
    }
}

void pla_index::encode_knots(const vector<int64_t> &brk_kval_vec){
    vector<int64_t> brk_uniform_diff_vec;
    uint64_t max_bits = shift_bits + lp_bits;
    
    max_data_ratio = ((1ULL << (max_bits)) - 1)/(dic_size - 1);
    
    for(int64_t i=0; i<dic_size; i++){
        int64_t diff = int64_t(brk_kval_vec[i]) - i*max_data_ratio;
        brk_uniform_diff_vec.emplace_back(diff);
    }
    brk_uniform_diff_packed.BuilidSignedPackedVector(brk_uniform_diff_vec);
}

void pla_index::SetBit(uint64_t &number, uint64_t pos){
    uint64_t mask = 1ULL << (64-pos); // pos = [1, 64]
    number |= mask;
}


int64_t pla_index::binary_search(int64_t lo, int64_t hi, const int64_t kval) const{
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

void pla_index::Load(string dic_fn, int64_t total_bits){
    std::ifstream is(dic_fn, std::ios::binary);
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
    cout<<"#breakpoints: "<<dic_size<<endl;
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
    essentials_arank::load_pod(is, isFirstIndexReturned);
    if(isFirstIndexReturned){
        int64_t bitvec_size = total_data_points/64 + 1;
        index_bitvec.resize(bitvec_size);
        is.read(reinterpret_cast<char*>(index_bitvec.data()), 
            static_cast<std::streamsize>(sizeof(uint64_t) * bitvec_size));  
        // err = fread(&index_bitvec[0], sizeof(uint64_t), bitvec_size, fp);
    }
}

void pla_index::save_unpacked(string dic_fn){
    dic_fn = dic_fn + ".unpacked";
    std::ofstream os(dic_fn, std::ios::binary);
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
    essentials_arank::save_pod(os, isFirstIndexReturned);
    if(isFirstIndexReturned){
        os.write(reinterpret_cast<char const*>(index_bitvec.data()),
                  index_bitvec.size() * sizeof(uint64_t));
    }
}

void pla_index::save_bit_info(string dic_fn, CBitPacking& ind_diff_pack){
    int pos = dic_fn.find(".bin");
    string pref = dic_fn.substr(0, pos);
    string suff = dic_fn.substr(pos);
    string bit_fn = pref+"_bitInfo"+suff;
    ofstream bit_info(bit_fn.c_str());
    bit_info<<"#Breakpoints: "<<dic_size<<endl;
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

uint64_t pla_index::get_size_in_bytes(){
    uint64_t total_size_in_bits = indirection_vec.size()*sizeof(uint64_t) + 
            uint64_t(brk_uniform_diff_packed.GetPacketSize())*brk_uniform_diff_packed.GetNumElements() +
            y_range_beg.num_bits();
    if(indx_type == BASIC_PLA){
        total_size_in_bits += uint64_t(diff_packd.GetPacketSize())*diff_packd.GetNumElements();
    }
    else if(indx_type == REPEAT_PLA){
        total_size_in_bits += size_in_bytes(sa_diff_dac_vec)*8.0;
    }
    if(isFirstIndexReturned) total_size_in_bits += index_bitvec.size()*sizeof(uint64_t);
    return total_size_in_bits/8;
}

void pla_index::Save(string dic_fn){
    std::ofstream os(dic_fn, std::ios::binary);
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
    essentials_arank::save_pod(os, isFirstIndexReturned);
    if(isFirstIndexReturned){
        os.write(reinterpret_cast<char const*>(index_bitvec.data()),
                  index_bitvec.size() * sizeof(uint64_t));
    }

    save_unpacked(dic_fn);
    save_bit_info(dic_fn, ind_diff_pack);
}