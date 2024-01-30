#include <iostream>

#include "../include/mm_file/mm_file.hpp"

int main() {
    std::string filename("./tmp.bin");
    static const size_t n = 13;

    {
        // write n uint32_t integers
        mm::file_sink<uint32_t> fout(filename, n);
        std::cout << "mapped " << fout.bytes() << " bytes "
                  << "for " << fout.size() << " integers" << std::endl;

        auto* data = fout.data();
        for (uint32_t i = 0; i != fout.size(); ++i) {
            data[i] = i;
            std::cout << "written " << data[i] << std::endl;
        }

        // test iterator
        for (auto x : fout) {
            std::cout << "written " << x << std::endl;
        }

        fout.close();
    }

    {
        // instruct the kernel that we will read the content
        // of the file sequentially
        int advice = mm::advice::sequential;

        // read the stream as uint16_t integers
        mm::file_source<uint16_t> fin1(filename, advice);
        std::cout << "mapped " << fin1.bytes() << " bytes "
                  << "for " << fin1.size() << " integers" << std::endl;
        auto const* data = fin1.data();
        for (uint32_t i = 0; i != fin1.size(); ++i) {
            std::cout << "read " << data[i] << std::endl;
        }
        fin1.close();

        // test iterator
        mm::file_source<uint16_t> fin2;
        fin2.open(filename, advice);
        for (auto x : fin2) {
            std::cout << "read " << x << std::endl;
        }
        fin2.close();
    }

    std::remove(filename.c_str());

    return 0;
}