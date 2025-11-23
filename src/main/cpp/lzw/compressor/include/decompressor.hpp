#pragma once

#include <unordered_map>
#include <string>
#include "utils.hpp"

//-----------------------------------------------------------------------------------------

namespace my_compress {

class decompressor_t final {
    private:
        std::unordered_map<int, std::string> dictionary_;
        const int init_size_of_dictionary_ = 256;

    public:
        decompressor_t(int init_size_of_dictionary = 256) :
                      init_size_of_dictionary_(init_size_of_dictionary){
            for (int i = 0; i < init_size_of_dictionary_; ++i) {
                dictionary_[i] = std::string{char(i)};
            }
        }
        template<typename it>
        inline std::string decompress(it start, it fin);
        size_t get_dict_size() const {return dictionary_.size();};
        size_t get_init_size_of_dictionary() const {return init_size_of_dictionary_;};
        std::vector<int> get_data(const char* path) const;
        inline void write_res(std::string res, const char* path) const;
};

//-----------------------------------------------------------------------------------------

template<typename it>
std::string decompressor_t::decompress(it start, it fin) {
    std::string output {};

    if (start == fin)
        return output;

    int code = init_size_of_dictionary_ + 1;
    int fst_elem_code = *start;

    std::string fst_elem = dictionary_[fst_elem_code];
    output += fst_elem;
    std::string snd_elem {};
    snd_elem = fst_elem;
    for (auto i = start; i != (fin - 1); ++i) {
        int snd_elem_code = *(i + 1);
        if (!dictionary_.contains(snd_elem_code)) {
			fst_elem = dictionary_[fst_elem_code] + snd_elem;
		}
        else {
            fst_elem = dictionary_[snd_elem_code];
        }
        output += fst_elem;
        snd_elem = "";
        snd_elem += fst_elem[0];
        dictionary_[code] = dictionary_[fst_elem_code] + snd_elem ;
        ++code;
        fst_elem_code = snd_elem_code;
    }
    // std::cout << output << "\n";
    return output;
}

void decompressor_t::write_res(std::string res, const char* path) const {
    std::ofstream out(path, std::ios::binary);

    out << res;
    out.close();
}

std::vector<int> decompressor_t::get_data(const char* path) const {

    std::ifstream data_file(path, std::ios::binary);
    if (!data_file.good()) {
        std::string msg = "Input file for decompression does not exist\n";
        throw std::runtime_error(msg + path);
    }

    size_t f_size = utils::get_file_size(path);
    char* char_buf = new char[f_size];
    data_file.read(char_buf, f_size);
    std::vector<int> data (reinterpret_cast<int*>(char_buf),
                           reinterpret_cast<int*>(char_buf) + f_size / sizeof(int));
    delete [] char_buf;
    data_file.close();
    return data;
}

}
