#include "BitPacking.h"

CBitPacking::CBitPacking(){}

/**
* Constructor be used for packing elements from a given vector
*
* @param isSigned whether the packet will have a sign bit or not
*/
CBitPacking::CBitPacking(bool isSigned){
    m_bIsDataSigned = isSigned;
}

/**
* @brief to be used from packed vector. This is useful when we 
* want to unpack values from a saved vector.
*
* @param packedData is the already packed vector
* @param packetSize is the number of bits needed to store the maximum of the original value
* @param numElements Number of elements that are packed in the packed data 
*/
CBitPacking::CBitPacking(vector<uint64_t> &packedData, uint64_t packetSize, 
        uint64_t numElemets){
    m_vPackedData = packedData;
    m_iPacketDataSize = packedData.size();
    m_iPacketSize = packetSize;
    m_iNumElements = numElemets;
    m_iLastIndex = -1;
    m_iLastPosition = -1;
}

void CBitPacking::Load(FILE* fp, uint64_t numElements){
    m_iNumElements = numElements;
    uint8_t t_pckt_size;
    int64_t err = fread(&t_pckt_size, sizeof(uint8_t), 1, fp);
    m_iPacketSize = t_pckt_size;
    this->m_iPacketDataSize = GetPacketDataSize();
    m_vPackedData.resize(m_iPacketDataSize);
    // cout<<"packet data size: "<<m_iPacketDataSize<<endl;
    err = fread(&m_vPackedData[0], sizeof(uint64_t), m_iPacketDataSize, fp);
    m_iLastIndex = -1;
    m_iLastPosition = -1;   
    m_sign_mask = ~UINT64_C(0) >> (64 - m_iPacketSize + 1);
}

void CBitPacking::Save(FILE* fp) const{
    uint8_t pckt_size = this->GetPacketSize();
    fwrite(&pckt_size, sizeof(pckt_size), 1, fp);    
    fwrite(&m_vPackedData[0], sizeof(uint64_t), m_vPackedData.size(), fp);
}

void CBitPacking::Save_os(std::ostream &os) const{
    uint8_t pckt_size = this->GetPacketSize();
    os.write(reinterpret_cast<char const*>(&pckt_size),
                  sizeof(pckt_size));
    os.write(reinterpret_cast<char const*>(m_vPackedData.data()),
                  m_vPackedData.size() * sizeof(uint64_t));
}

void CBitPacking::Load_is(std::istream &is, uint64_t numElements){
    m_iNumElements = numElements;
    uint8_t t_pckt_size;
    is.read(reinterpret_cast<char*>(&t_pckt_size), sizeof(t_pckt_size));
    m_iPacketSize = t_pckt_size;
    this->m_iPacketDataSize = GetPacketDataSize();
    m_vPackedData.resize(m_iPacketDataSize);
    // cout<<"packet data size: "<<m_iPacketDataSize<<endl;
    is.read(reinterpret_cast<char*>(m_vPackedData.data()), 
        static_cast<std::streamsize>(sizeof(uint64_t) * m_iPacketDataSize));
    m_iLastIndex = -1;
    m_iLastPosition = -1;   
    m_sign_mask = ~UINT64_C(0) >> (64 - m_iPacketSize + 1);
}

vector<uint64_t> CBitPacking::GetPackedVector() const{
    return m_vPackedData;
}

uint8_t CBitPacking::GetPacketSize() const{
    return m_iPacketSize;
}

void CBitPacking::SetIsDataSigned(bool val){
    m_bIsDataSigned = val;
}

uint64_t CBitPacking::GetNumElements() const{
    return m_iNumElements;
}

uint64_t CBitPacking::GetPacketDataSize() const{
    return ceil((double)(m_iPacketSize * m_iNumElements)/64);
}

void CBitPacking::BuilidPackedVector(vector<uint64_t> &dataVec){
    uint64_t largest = *max_element(dataVec.begin(), dataVec.end());
    if(largest < 2) m_iPacketSize = 1;
    else {
        if(ceil(log2(largest)) != floor(log2(largest)))
            m_iPacketSize = ceil(log2(largest));
        else
            m_iPacketSize = ceil(log2(largest)) + 1;
    }
    // cout<<"Packet size: "<<largest<<" "<<unsigned(m_iPacketSize)<<endl;
    m_iNumElements = dataVec.size();
    m_iPacketDataSize = GetPacketDataSize();
    // cout<<"pcktData size: "<<m_iPacketDataSize<<" #Elements: "<<m_iNumElements<<endl;
    m_vPackedData.resize(m_iPacketDataSize, 0); // initialize everyone to 0
    m_iLastIndex = 0;
    m_iLastPosition = 0;
    for(size_t i=0; i<dataVec.size(); i++){
        SetBits(dataVec[i]);
    }
}


void CBitPacking::BuilidSignedPackedVector(vector<int64_t> &dataVec){
    int64_t _max = *max_element(dataVec.begin(), dataVec.end());
    int64_t _min = *min_element(dataVec.begin(), dataVec.end());
    if(_min < 0) _min *= -1;
    uint64_t max_elem = max(_max, _min);

    m_iNumElements = dataVec.size();
    if(ceil(log2(max_elem)) != floor(log2(max_elem)))
        m_iPacketSize = ceil(log2(max_elem)) + 1; // + 1 for sign bit
    else
        m_iPacketSize = ceil(log2(max_elem)) + 2; // +2 for sign bit and main value is a power of 2
    // m_iPacketSize = ceil(log2(max_elem))+1; // +1 for sign bit
    m_iPacketDataSize = GetPacketDataSize();
    // cout<<"pcktData: "<<m_iPacketDataSize<<endl;
    m_vPackedData.resize(m_iPacketDataSize, 0);
    m_iLastIndex = 0;
    m_iLastPosition = 0;
    
    uint64_t data;
    bool isNeg;
    // uint64_t put_sign_mask = (uint64_t) 1<< (m_iPacketSize-1);
    uint64_t rmv_sign_bitmask = (uint64_t)1<<63;
    rmv_sign_bitmask -= 1;
    int64_t val;
    for(size_t i=0; i<dataVec.size(); i++){
        isNeg = false;
        val = dataVec[i];
        if(val < 0){
            isNeg = true;
            val*=-1;
        } 
        val <<= 1;
        data = val & rmv_sign_bitmask;
        if (isNeg) data = val | 1;
        
        // if(isNeg) data |= put_sign_mask;
        
        SetBits(data);
    }
    m_sign_mask = ~UINT64_C(0) >> (64 - m_iPacketSize + 1);
}

void CBitPacking::SetBits(uint64_t value){    
    int64_t mask;
    if(64-m_iLastPosition >= m_iPacketSize){
        // All bits can be stored at the same element of the vector
        mask = value << (64-(m_iLastPosition)- m_iPacketSize);
        
        m_vPackedData[m_iLastIndex] |= mask;
        // cout<<value <<" "<<bitset<64> (m_vPackedData[m_iLastIndex])<<endl;
        m_iLastPosition += m_iPacketSize;
        if(m_iLastPosition >= 64){
            m_iLastPosition -= 64;
            m_iLastIndex++;
        }
        
    }
    else{
        // use the remaining bits of the current element
        int32_t remBits = 64 - (m_iLastPosition);
        mask = value >> (m_iPacketSize-remBits);
        m_vPackedData[m_iLastIndex] |= mask;
        m_iLastIndex++;
        // mask the bits of the next element with the remaining bits
        mask = value << (64 - m_iPacketSize + remBits);
        m_vPackedData[m_iLastIndex] |= mask;
        m_iLastPosition = m_iPacketSize - remBits;
    }
}

vector<int64_t> CBitPacking::GetUnpackedVector() const{
    vector<int64_t> unpacked(m_iNumElements, 0);
    for(uint64_t i=0; i<m_iNumElements; i++){
        unpacked[i] = GetValueAt(i);
    }
    return unpacked;
}

/**
 * Assumes the value will fit within 63 bit
 * if the packed data is signed, then returns
 * negative depending on the sign bit. In that case
 * uses the remaining bits to calculate the value
 * Otherwise returns the full value
 * 
 * @param pckt the packet number where actual value is
 * compacted
 * @returns actual value encoded in the pckt
*/

int64_t CBitPacking::GetValueAt(uint64_t pckt) const {
    uint64_t packedDataIndex = (pckt * m_iPacketSize) >> 6;  // div 2^6 = 64
    uint64_t value, packedValue, elementPos, mask;

    // Calculate common expressions outside the conditional blocks
    packedValue = m_vPackedData[packedDataIndex];
    elementPos = pckt * m_iPacketSize - (packedDataIndex << 6);

    // Use a more efficient mask
    mask = ~UINT64_C(0) >> elementPos;
    value = packedValue & mask;
    // Simplify the conditional check
    bool fullPacketInside = (pckt + 1) * m_iPacketSize <= ((packedDataIndex + 1) << 6);
    
    if (fullPacketInside) {
        value >>= (64 - elementPos - m_iPacketSize);
    } else {        
        uint64_t remBits = m_iPacketSize - (64 - elementPos);
        packedValue = m_vPackedData[packedDataIndex + 1];

        mask = packedValue >> (64 - remBits);
        value <<= remBits;
        value |= mask;
    }

    if (!m_bIsDataSigned) return value;

    // uint64_t isNegative = value >> (m_iPacketSize - 1);
    bool isNegative = value & 1;
    
    int64_t actualVal = value >> 1;

    return isNegative ? -actualVal : actualVal;
}
