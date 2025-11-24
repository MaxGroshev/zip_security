#pragma once

#include "time_control.hpp"
#include "compressor.hpp"
#include "decompressor.hpp"
#include "config.hpp"

//-----------------------------------------------------------------------------------------

namespace my_compress {

class lzw_t final {
    private:
        set_up_t config_ = {};
        compressor_t compressor;
        decompressor_t decompressor;

    public:
        lzw_t(set_up_t config, uint32_t init_encode_dictionary_size,
                               uint32_t init_decode_dictionary_size) :
            config_(config),
            compressor(init_encode_dictionary_size),
            decompressor(init_decode_dictionary_size) {}

        lzw_t(set_up_t config) :
            config_(config),
            compressor(config.get_dict_dir()),
            decompressor(0) {}


//COMPRESSOR_WRAPER------------------------------------------------------------------------

        std::vector<uint32_t> compress(const std::string& input_dat) {
            return compressor.compress(input_dat.begin(), input_dat.end());
        };
        std::vector<uint32_t> compress() {
            std::string input_dat = get_data_for_compress();
            std::vector<uint32_t> res = compressor.compress(input_dat.begin(), input_dat.end());
            std::string write_dir = config_.get_dest_dir();
            compressor.write_res_in_bin(res, write_dir.c_str());
            return res;
        };

//DECOMPRESSOR_WRAPER----------------------------------------------------------------------

        std::string decompress(const std::vector<uint32_t>& input_dat) {
            return decompressor.decompress(input_dat.begin(), input_dat.end());
        };
        std::string decompress() {
            std::vector<uint32_t> input_dat = get_data_for_decompress();
            utils::dump_vect(input_dat);
            std::string res = decompressor.decompress(input_dat.begin(), input_dat.end());
            decompressor.write_res(res, config_.get_dest_dir());
            return res;
        };

//GETTERS----------------------------------------------------------------------------------

        size_t get_compress_dict_size()     const {return compressor.get_dict_size();};
        size_t get_decompress_dict_size()   const {return decompressor.get_dict_size();};
        bool   is_train_mode()              const {return config_.is_train_mode();};
        bool   is_dict_mode()               const {return config_.is_dict_mode();};

        std::string get_data_for_compress() const {
            return utils::read_from_file_into_string(config_.get_src_dir());
        };
        std::vector<uint32_t> get_data_for_decompress() const {
            return decompressor.get_data(config_.get_src_dir());
        };
        void compress_show_res() const {
            size_t comp_size =  utils::get_file_size(config_.get_dest_dir());
            size_t input_size =  utils::get_file_size(config_.get_src_dir());
            std::clog << "File compressed uint32_to: "<< config_.get_dest_dir() << std::endl;
            std::clog << input_size / 1024 << " Kbts" <<
                        " ==> " << comp_size / 1024 << " Kbts" << "\n";
            std::clog << "Compression factor: " << double(comp_size * 100) /
                                                   double(input_size) << " %\n";
        }
};
}
//-----------------------------------------------------------------------------------------

