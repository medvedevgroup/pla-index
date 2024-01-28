#pragma once

#include <vector>
#include <algorithm>
#include <string>
#include <cstdlib>
#include <chrono>
#include <cmath>
#include <iostream>
using namespace std;


/**
* Can be used to bitpack a vector of 32 bit values into
* a shorter vector of 64 bit.
* Also, the packed vector can be used at the constructor
* to retrieve the packed values.
*
* Each 32 bit value is packed into log(max_val) bits.
* Can be improved using a boundary vector.
*/
// #pragma pack(push, 1)
class CBitPacking{
private:
    uint64_t m_iPacketSize; // How many bits needed to store the largest element
    uint64_t m_iPacketDataSize;
    uint64_t m_iLastIndex;
    uint64_t m_iLastPosition;
    uint64_t m_iNumElements;
    uint64_t m_sign_mask;
    bool m_bIsDataSigned;
    vector<uint64_t> m_vPackedData;
    void SetBits(uint64_t value);
public:
    CBitPacking(); // default one
    CBitPacking(bool isSigned);
    CBitPacking(vector<uint64_t> &packedData, uint64_t packetSize, uint64_t numElemets);
    void BuilidPackedVector(vector<uint64_t> &dataVec);
    void BuilidSignedPackedVector(vector<int64_t> &dataVec);
    int64_t GetValueAt(uint64_t pckt) const;
    vector<uint64_t> GetPackedVector() const;
    vector<int64_t> GetUnpackedVector() const;
    uint64_t GetNumElements() const;
    uint8_t GetPacketSize() const;
    uint64_t GetPacketDataSize() const;
    void SetIsDataSigned(bool val);
    void Load(FILE* fp, uint64_t numElements);
    void Save(FILE* fp) const;
    void Save_os(std::ostream &os) const;
    void Load_is(std::istream &is, uint64_t numElements);
};
// #pragma pack(pop)
