#include "pla_index.hpp"

// constructor at index build
pla_index::pla_index(int64_t eps, int64_t lp_bits, uint64_t largest,
            bool isFirstIndxRet, uint64_t total_points, INDX_TYPE it):
            epsilon(eps),
            lp_bits(lp_bits), isFirstIndexReturned(isFirstIndxRet),
            indx_type(it),
            total_data_points(total_points),
            brk_uniform_diff_packed(true), diff_packd(true),
            err_dict_packed(true)
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
    brk_uniform_diff_packed(true), diff_packd(true),
    err_dict_packed(true)
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
    uint64_t total_bytes = essentials::load(f, dic_fn.c_str());
    std::ifstream is(dic_fn, std::ios::binary);
    std::streampos current_position = is.tellg();
    is.seekg(total_bytes, std::ios::beg);

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
    diff_packd.Load_is(is, dic_size-1);    
    
    y_range_beg.load(is);
    // cout<<"max err load\n";
    uint16_t max_error_t;
    essentials_arank::load_pod(is, max_error_t);
    epsilon = (int64_t)max_error_t;

    uint64_t unique_kmers;
    essentials_arank::load_pod(is, unique_kmers);
    err_dict_packed.Load_is(is, unique_kmers);
    
}

void pla_index::save_unpacked(string dic_fn){
    dic_fn = dic_fn + ".unpacked";
    essentials::save(f, dic_fn.c_str());
    std::ofstream os(dic_fn, std::ios_base::app);
    uint8_t t_lp_bits = lp_bits;
    essentials_arank::save_pod(os, t_lp_bits);
    os.write(reinterpret_cast<char const*>(indirection_vec.data()),
                  indirection_vec.size() * sizeof(uint64_t));
    uint32_t nBrkpnts = dic_size;
    essentials_arank::save_pod(os, nBrkpnts);
    diff_packd.Save_os(os);
    y_range_beg.save(os);
    uint16_t t_err_th = epsilon;
    essentials_arank::save_pod(os, t_err_th);
    uint64_t uniq_kmers = err_dict_packed.GetNumElements();
    essentials_arank::save_pod(os, uniq_kmers);
    err_dict_packed.Save_os(os);
}

void pla_index::save_bit_info(string dic_fn, CBitPacking& ind_diff_pack){
    int pos = dic_fn.find(".bin");
    string pref = dic_fn.substr(0, pos);
    string suff = dic_fn.substr(pos);
    string bit_fn = pref+"_bitInfo"+suff;
    ofstream bit_info(bit_fn.c_str());
    bit_info<<"#Breakpoints: "<<dic_size<<endl;
    bit_info<<"#Unique kmers: "<<err_dict_packed.GetNumElements()<<endl;
    bit_info<<"Indirection vector: each: "<<uint64_t(ind_diff_pack.GetPacketSize())<<endl;
    bit_info<<"knot kmers: each: "<<uint64_t(brk_uniform_diff_packed.GetPacketSize())<<endl;    
    bit_info<<"Y beg each: "<<ceil(y_range_beg.num_bits()/dic_size)<<endl;
    bit_info<<"Y end diff: each: "<<uint64_t(diff_packd.GetPacketSize())<<endl;        
    bit_info<<"True err each: "<<uint64_t(err_dict_packed.GetPacketSize())<<endl;
    bit_info<<"MPHF bits per key: "<<static_cast<double>(f.num_bits()) / f.num_keys()<<endl;;

    bit_info.close();
}

uint64_t pla_index::get_total_size_in_bytes(){
    uint64_t total_size_in_bits = indirection_vec.size()*sizeof(uint64_t) +
            uint64_t(brk_uniform_diff_packed.GetPacketSize())*brk_uniform_diff_packed.GetNumElements() +
            y_range_beg.num_bits() + uint64_t(diff_packd.GetPacketSize())*diff_packd.GetNumElements() +
            uint64_t(err_dict_packed.GetPacketSize())*err_dict_packed.GetNumElements() +
            f.num_bits();
    return total_size_in_bits / 8;
}

void pla_index::Save(string dic_fn){
    essentials::save(f, dic_fn.c_str());
    std::ofstream os(dic_fn, std::ios_base::app);
    uint8_t t_lp_bits = lp_bits;
    essentials_arank::save_pod(os, t_lp_bits);
    vector<uint64_t> ind_diff_vec = get_processed_ind_vec();
    CBitPacking ind_diff_pack(false);
    ind_diff_pack.BuilidPackedVector(ind_diff_vec);
    ind_diff_pack.Save_os(os);
    uint32_t nBrkpnts = dic_size;
    
    essentials_arank::save_pod(os, nBrkpnts);
    
    brk_uniform_diff_packed.Save_os(os);
    diff_packd.Save_os(os);
    y_range_beg.save(os);
    uint16_t t_err_th = epsilon;
    essentials_arank::save_pod(os, t_err_th);

    uint64_t uniq_kmers = err_dict_packed.GetNumElements();
    essentials_arank::save_pod(os, uniq_kmers);
    err_dict_packed.Save_os(os);


    // save_unpacked(dic_fn);
    // save_bit_info(dic_fn, ind_diff_pack);
}