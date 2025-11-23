#pragma once

#include <unordered_map>
#include "config.hpp"
#include <algorithm>
#include "utils.hpp"

//-----------------------------------------------------------------------------------------

class dict_train_t {

    using um_vector = typename std::vector<std::pair<std::string, int>>;

    private:
        set_up_t config_;
        std::unordered_map<std::string, int> dictionary_;
    public:
        dict_train_t(set_up_t config) : config_(config) {}

        void train() {
            std::string train_dir = config_.get_src_dir();
            std::vector<std::string> samples =
                                        utils::get_file_name_from_dir(train_dir.c_str());

            for (const auto& sample : samples) {
                std::string sample_dat = utils::read_from_file_into_string(sample.c_str()); //const char* ?
                train_on_sample(sample_dat.begin(), sample_dat.end());
            }

            um_vector data = sort_dictionary();
            write_dictionary_to_file(data);
        };

        template<typename it>
        void train_on_sample(it start, it fin) {

            std::string fst_elem {};
            std::string snd_elem {};
            std::string concat   {};
            fst_elem += *start;
            for (auto i = start; i != fin; ++i) {
                if (i != fin)
                    snd_elem += *(i + 1);
                concat = fst_elem + snd_elem;

                auto um_elem = dictionary_.find(concat);
                if (um_elem != dictionary_.end()) {
                    fst_elem = concat;
                    um_elem->second++;
                }
                else {
                    dictionary_[concat] = 1;
                    fst_elem = snd_elem;
                }
                snd_elem = "";
            }
        }

        um_vector sort_dictionary() {
            um_vector data (dictionary_.begin(), dictionary_.end());
            std::sort(data.begin(), data.end(),
                      [](auto fst, auto snd){return (fst.second > snd.second);});

            return data;
        }

        void write_dictionary_to_file(um_vector& data) const {
            std::ofstream out(config_.get_dest_dir(), std::ios::binary);

            for (int i = 0; i < std::min(int(data.size()), config_.get_maxdict()); i++) {
                data[i].first += "\n";
                out.write(data[i].first.data(), data[i].first.size());//probably too low
            }
            out.close();
        }

//-----------------------------------------------------------------------------------------

};
