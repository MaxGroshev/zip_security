#include "lzw.hpp"
#include "utils.hpp"
#include "dict_train.hpp"

using namespace my_compress;

//-----------------------------------------------------------------------------------------

int main(int argc, char** argv) {

    auto start_time = time_control::chrono_cur_time ();
    try {
        set_up_t set_up {argc, argv};

        if (set_up.is_train_mode()) {
            dict_train_t trainer {set_up};
            trainer.train();
        }
        else if (set_up.is_dict_mode()) {
            lzw_t dict_encoder {set_up};
            std::vector<int> res = dict_encoder.compress();
            dict_encoder.compress_show_res();
        }
        else {
            int encoder_init_size = 256;
            int decoder_init_size = 0;
            lzw_t encoder {set_up, encoder_init_size, decoder_init_size};
            std::vector<int> res = encoder.compress();
            encoder.compress_show_res();
        }

    } catch(std::runtime_error& err) {
        std::cout << "ERROR: " << err.what() << '\n';
    } catch(...) {
        std::cout << "Unexpected error\n";
    }
    auto end_time = time_control::chrono_cur_time ();
    time_control::show_run_time(start_time, end_time, "Compression wall time: ");

    return 0;
}
