
#pragma once

#include <gtest/gtest.h>

#include "chacha20.hpp"
#include "key_gen_test.hpp"
#include "key_generator.hpp"

using namespace my_compress;

//-----------------------------------------------------------------------------------------

class keygen_test : public ::testing::Test {
    protected:
        uint32_t counter = 1;
        std::string text = "Ladies and gentlemen of the class of '99: If I could offer you only one tip for the future, sunscreen would be it.";

    void SetUp() {
    }
};

//-----------------------------------------------------------------------------------------

TEST_F(keygen_test, compare_generated_and_translated_keys) {

    key_generator_t<uint32_t> gn{};
    auto key1 = gn.generate();
    auto str_key = gn.uint8_vector_to_hex_string(key1);
    auto key = gn.hex_string_to_uint8_vector(str_key);
    
    // std::cout <<  str_key << std::endl;
    // for(auto byte : key1) {
    //     std::cout << std::hex << std::setw(2) << std::setfill('0') << (int)byte;
    // }
    // std::cout << std::dec << std::endl;
    
    // for (size_t i = 0; i < myString.length(); ++i) {
    //     byteArray[i] = static_cast<uint8_t>(myString[i]);
    // }

    std::array<uint32_t, 12> nonce = {0,0,0,0,0,0,0,0x4a,0,0,0,0};
    std::vector<uint32_t> plaintext(text.begin(), text.end());

    chacha20_t chacha(key, nonce, counter);
    std::vector<uint32_t> ciphertext;
    chacha.crypt(plaintext, ciphertext);

    chacha20_t chacha_decrypt(key, nonce, counter);
    std::vector<uint32_t> deciphertext;
    chacha_decrypt.crypt(ciphertext, deciphertext);


    ASSERT_TRUE(plaintext == deciphertext);
}