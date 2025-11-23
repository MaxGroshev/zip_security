#pragma once

#include <unordered_map>
#include <string>
#include "utils.hpp"
#include <iterator>

//-----------------------------------------------------------------------------------------

namespace my_compress {

class compressor_t final {
    private:
        std::unordered_map<std::string, int> dictionary_;
        const int init_size_of_dictionary_ = 256;

    public:
        compressor_t(int init_size_of_dictionary = 256) :
                        init_size_of_dictionary_(init_size_of_dictionary){
            for (int i = 0; i < init_size_of_dictionary_; ++i) {
                dictionary_[std::string{char(i)}] = i;
            }
        }

        compressor_t(const char* path) {
            std::vector<std::string> dict = get_dict(path);
            for (int i = 0; i < dict.size(); ++i) {
                dictionary_[dict[i]] = i;
            }
        }
        template<typename it>
        inline std::vector<int> compress(it start, it fin);
        template<typename it>
        inline std::vector<int>encrypt(it start, it end);

        size_t get_dict_size() const {return dictionary_.size();};
        size_t get_init_size_of_dictionary() const {return init_size_of_dictionary_;};

        inline std::string get_data(const char* path) const;
        inline std::vector<std::string> get_dict(const char* path) const;
        inline void write_res_in_bin(std::vector<int> res, const char* path) const;

};

//-----------------------------------------------------------------------------------------

template<typename it>
std::vector<int> compressor_t::compress(it start, it fin) {
    std::vector<int> output;

    if (start == fin)
        return output;

    std::string fst_elem {};
    std::string snd_elem {};
    std::string concat   {};
    fst_elem += *start;
    int code = dictionary_.size() + 1;
    for (auto i = start; i != fin; ++i) {
        if (i != fin)
            snd_elem += *(i + 1);
        concat = fst_elem + snd_elem;
        if (dictionary_.contains(concat) && *i != 0) {
            fst_elem = concat;
        }
        else {
            output.push_back(dictionary_[fst_elem]);
            dictionary_[concat] = code;
            ++code;
            fst_elem = snd_elem;
        }
        snd_elem = "";
    }

    return output;
}

//-----------------------------------------------------------------------------------------

void compressor_t::write_res_in_bin(std::vector<int> res, const char* path) const {
    std::ofstream out(path, std::ios::binary);

    for (auto& elem : res)
        out.write((char*)&elem, sizeof(int));
    out.close();
}

std::vector<std::string> compressor_t::get_dict(const char* path) const {

    std::ifstream dict_file(path);
    if (!dict_file.good()) {
        std::string msg = "Dictionary file does not exist\n";
        throw std::runtime_error(msg + path);
    }

    std::vector<std::string> data;
    std::string tmp_str = {};
    while (std::getline(dict_file, tmp_str)) {
        data.push_back(std::move(tmp_str));
    }
    dict_file.close();

    return data;
}

}
