#ifndef STROBEALIGN_IO_HPP
#define STROBEALIGN_IO_HPP

#include <iostream>
#include <vector>
#include <cstdint>

void write_int_to_ostream(std::ostream& os, int32_t value);
int32_t read_int_from_istream(std::istream& is);

/* Write a vector to an output stream, preceded by its length */
template <typename T>
void write_vector(std::ostream& os, const std::vector<T>& v) {
    uint64_t size = v.size();
    os.write(reinterpret_cast<char*>(&size), sizeof(size));
    os.write(reinterpret_cast<const char*>(v.data()), v.size() * sizeof(T));
}

/* Read a vector written by write_vector */
template <typename T>
void read_vector(std::istream& is, std::vector<T>& v) {
    uint64_t size;
    v.clear();
    is.read(reinterpret_cast<char*>(&size), sizeof(size));
    v.resize(size);
    is.read(reinterpret_cast<char*>(v.data()), size * sizeof(T));
}

#endif
