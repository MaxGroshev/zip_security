#pragma once

#include <iostream>
#include <vector>
#include <array>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <unordered_map>
#include <string>
#include <iterator>

//-----------------------------------------------------------------------------------------

namespace my_compress {

class chacha20_t {
    public:
        static const int KEY_SIZE = 32;
        static const int NONCE_SIZE = 12;
        static const int BLOCK_SIZE = 64;

        chacha20_t(const std::array<uint32_t, KEY_SIZE>& key, 
                const std::array<uint32_t, NONCE_SIZE>& nonce, 
                uint32_t counter = 1) {
            state_setup(key, nonce, counter);
        }

        void crypt(const std::vector<uint32_t>& input, std::vector<uint32_t>& output) {
            output.resize(input.size());
            std::vector<uint32_t> key_stream_block(BLOCK_SIZE);
            
            size_t numBlocks = input.size() / BLOCK_SIZE;
            size_t remainder = input.size() % BLOCK_SIZE;

            for (size_t i = 0; i < numBlocks; ++i) {
                generate_block(key_stream_block);
                for (size_t j = 0; j < BLOCK_SIZE; ++j) {
                    output[i * BLOCK_SIZE + j] = input[i * BLOCK_SIZE + j] ^ key_stream_block[j];
                }
            }

            if (remainder > 0) {
                generate_block(key_stream_block);
                for (size_t j = 0; j < remainder; ++j) {
                    output[numBlocks * BLOCK_SIZE + j] = input[numBlocks * BLOCK_SIZE + j] ^ key_stream_block[j];
                }
            }
        }

    private:
        uint32_t state[16];

        static inline uint32_t rotl(uint32_t x, int n) {
            return (x << n) | (x >> (32 - n));
        }

        static inline uint32_t pack4(const uint32_t* a) {
            return uint32_t(a[0]) | (uint32_t(a[1]) << 8) | 
                (uint32_t(a[2]) << 16) | (uint32_t(a[3]) << 24);
        }

        static inline void unpack4(uint32_t src, std::vector<uint32_t>& dst, size_t offset) {
            dst[offset + 0] = (src >> 0) & 0xff;
            dst[offset + 1] = (src >> 8) & 0xff;
            dst[offset + 2] = (src >> 16) & 0xff;
            dst[offset + 3] = (src >> 24) & 0xff;
        }

        void state_setup(const std::array<uint32_t, KEY_SIZE>& key, 
                        const std::array<uint32_t, NONCE_SIZE>& nonce, 
                        uint32_t counter) {
            state[0] = 0x61707865;
            state[1] = 0x3320646e;
            state[2] = 0x79622d32;
            state[3] = 0x6b206574;

            for (int i = 0; i < 8; ++i) {
                state[4 + i] = pack4(&key[i * 4]);
            }

            state[12] = counter;

            for (int i = 0; i < 3; ++i) {
                state[13 + i] = pack4(&nonce[i * 4]);
            }
        }

        // The Quarter Round function
        static void quarter_round(uint32_t& a, uint32_t& b, uint32_t& c, uint32_t& d) {
            a += b; d ^= a; d = rotl(d, 16);
            c += d; b ^= c; b = rotl(b, 12);
            a += b; d ^= a; d = rotl(d, 8);
            c += d; b ^= c; b = rotl(b, 7);
        }

        void generate_block(std::vector<uint32_t>& outputBlock) {
            uint32_t workingState[16];
            std::memcpy(workingState, state, sizeof(state));

            for (int i = 0; i < 10; ++i) {
                quarter_round(workingState[0], workingState[4], workingState[8],  workingState[12]);
                quarter_round(workingState[1], workingState[5], workingState[9],  workingState[13]);
                quarter_round(workingState[2], workingState[6], workingState[10], workingState[14]);
                quarter_round(workingState[3], workingState[7], workingState[11], workingState[15]);

                quarter_round(workingState[0], workingState[5], workingState[10], workingState[15]);
                quarter_round(workingState[1], workingState[6], workingState[11], workingState[12]);
                quarter_round(workingState[2], workingState[7], workingState[8],  workingState[13]);
                quarter_round(workingState[3], workingState[4], workingState[9],  workingState[14]);
            }

            for (int i = 0; i < 16; ++i) {
                workingState[i] += state[i];
            }

            for (int i = 0; i < 16; ++i) {
                unpack4(workingState[i], outputBlock, i * 4);
            }

            state[12]++;
        }
};
}
