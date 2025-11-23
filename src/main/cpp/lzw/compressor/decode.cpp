#include "lzw.hpp"
#include "utils.hpp"

using namespace my_compress;

//-----------------------------------------------------------------------------------------

int main(int argc, char** argv) {

    auto start_time = time_control::chrono_cur_time ();
    try {
        set_up_t set_up {argc, argv};
        int encoder_init_size = 0;
        int decoder_init_size = 256;
        lzw_t decoder {set_up, encoder_init_size, decoder_init_size};
        std::string res = decoder.decompress();
    } catch(std::runtime_error& err) {
        std::cout << "ERROR: " << err.what() << '\n';
    } catch(...) {
        std::cout << "Unexpected error\n";
    }
    auto end_time = time_control::chrono_cur_time ();
    time_control::show_run_time(start_time, end_time, "Decompression wall time: ");

    return 0;
}
