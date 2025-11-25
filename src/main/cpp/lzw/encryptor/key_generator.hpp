#pragma once

#include <cassert>
#include <iostream>
#include <vector>
#include <array>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <unordered_map>
#include <string>
#include <iterator>
#include <random>
#include "array"
#include <sstream>
#include <iomanip> 

//-----------------------------------------------------------------------------------------

namespace my_compress {

template<typename T>
class key_generator_t {
    
private:
    constexpr static ssize_t key_size = 32;

public:
    std::array<T, key_size> generate() const {
        std::random_device rd;
        std::uniform_int_distribution<int> dist(0, 255);
        
        std::array<T, key_size> key{};
        for (size_t i = 0; i < key_size; ++i) {
            key[i] = static_cast<T>(dist(rd));
        }
        return key;
    }

    std::string uint8_vector_to_hex_string(const std::array<T, key_size>& data) const {
        std::stringstream ss;
        ss << std::hex << std::setfill('0');

        for (uint8_t byte : data) {
            ss << std::setw(2) << static_cast<int>(byte);
        }
        return ss.str();
    }

    std::array<T, key_size> hex_string_to_uint8_vector(const std::string& str_key) const  {
        assert(str_key.length == key_size * 2);

        std::array<T, key_size> key{};
        T high_part = 0; 
        T low_part  = 0; 

        for (int i = 0; i < key_size * 2; i += 2) {
            high_part = hex_char_to_int(str_key[i]);
            low_part = hex_char_to_int(str_key[i + 1]);
            assert(high_part != -1 && low_part != -1 && "incorrect hex\n");
            
            T byte_value = static_cast<T>((high_part << 4) | low_part);
            auto arr_index = i / 2;
            
            key[arr_index] = byte_value;
        }
        
        return key;
    }

    private: 
        T hex_char_to_int(char c) const {
            if (c >= '0' && c <= '9') {
                return c - '0';
            } else if (c >= 'a' && c <= 'f') {
                return c - 'a' + 10;
            } else if (c >= 'A' && c <= 'F') {
                return c - 'A' + 10;
            }
            return -1; // Invalid hex character
        }
};

}