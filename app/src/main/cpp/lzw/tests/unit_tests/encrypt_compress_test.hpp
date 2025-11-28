#pragma once

#include <gtest/gtest.h>

#include "lzw.hpp"
#include "chacha20.hpp"

using namespace my_compress;

//-----------------------------------------------------------------------------------------

class encrypt_test : public ::testing::Test {
    protected:
        uint32_t counter = 1;
        std::string text = R"("He deals the cards as a meditation
                            And those he plays never suspect
                            He doesn't play for the money he wins
                            He doesn't play for respect

                            [Verse 2]
                            He deals the cards to find the answer
                            The sacred geometry of chance
                            The hidden law of a probable outcome
                            The numbers lead a dance

                            [Chorus]
                            I know that the spades are the swords of a soldier
                            I know that the clubs are weapons of war
                            I know that diamonds mean money for this art
                            But that's not the shape of my heart"})";
    void SetUp() {
    }
};

//-----------------------------------------------------------------------------------------

TEST_F(encrypt_test, encrypt_test1) {


    lzw_t lzw{set_up_t{}, 256, 256};
    auto compressed = lzw.compress(text);

    std::array<uint32_t, 32> key;
    for(uint32_t i=0; i<32; i++) {
        key[i] = i;
    }
    std::array<uint32_t, 12> nonce = {0,0,0,0,0,0,0,0x4a,0,0,0,0};

    chacha20_t chacha(key, nonce, counter);
    std::vector<uint32_t> ciphertext;
    chacha.crypt(compressed, ciphertext);

    chacha20_t chacha_decrypt(key, nonce, counter);
    std::vector<uint32_t> deciphertext;
    chacha_decrypt.crypt(ciphertext, deciphertext);
    
    auto res = lzw.decompress(deciphertext);
    
    // std::cout << res << std::endl;
    // std::cout << compressed.size() << " : " << ciphertext.size() << std::endl;
    ASSERT_TRUE(compressed.size() == ciphertext.size());
    ASSERT_TRUE(res == text);
};