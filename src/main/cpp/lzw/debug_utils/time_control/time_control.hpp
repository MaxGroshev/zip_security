#pragma once

#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <cstdlib>
#include <chrono>
#include <iomanip>

//-----------------------------------------------------------------------------------------

namespace time_control {

using namespace std::literals;
using chrono_time_type = std::chrono::time_point<std::chrono::steady_clock>;

inline auto chrono_cur_time() {
    return std::chrono::steady_clock::now ();
};

inline auto chrono_time_in_time_t() {
    const auto now = std::chrono::system_clock::now();

    std::time_t t_now  = std::chrono::system_clock::to_time_t(now);
    return std::put_time(std::localtime(&t_now), "%F %T.\n");
};

inline void show_run_time(const chrono_time_type& start,
                          const chrono_time_type& end, const char* msg = nullptr) {
    std::clog << "----------------------------------------------\n";

    if (!msg)
        std::clog << "Total run time: ";
    else
        std::clog << msg;
    std::clog << ((end - start)) /1ms << " ms\n";
    std::clog << "----------------------------------------------\n";
};

}
