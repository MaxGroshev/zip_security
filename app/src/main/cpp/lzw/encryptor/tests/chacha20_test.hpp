#pragma once

#include <gtest/gtest.h>

#include "chacha20.hpp"

using namespace my_compress;

//-----------------------------------------------------------------------------------------

class encrypt_test : public ::testing::Test {
    protected:
        uint32_t counter = 1;
        std::string text = "Ladies and gentlemen of the class of '99: If I could offer you only one tip for the future, sunscreen would be it.";

    void SetUp() {
    }
};

//-----------------------------------------------------------------------------------------

TEST_F(encrypt_test, encrypt_test1) {

    std::array<uint32_t, 32> key;
    for(int i = 0; i < 32; i++) {
        key[i] = i;
    }
    std::array<uint32_t, 12> nonce = {0,0,0,0,0,0,0,0x4a,0,0,0,0};
    std::vector<uint32_t> plaintext(text.begin(), text.end());

    chacha20_t chacha(key, nonce, counter);
    std::vector<uint32_t> ciphertext;
    chacha.crypt(plaintext, ciphertext);

    chacha20_t chacha_decrypt(key, nonce, counter);
    std::vector<uint32_t> deciphertext;
    chacha_decrypt.crypt(ciphertext, deciphertext);

    std::string recoveredStr(deciphertext.begin(), deciphertext.end());
    std::cout << recoveredStr << std::endl;

    ASSERT_TRUE(plaintext == deciphertext);
}