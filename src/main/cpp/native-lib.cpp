#include <string>
#include <iostream>
#include <algorithm>
#include <cerrno>
#include <cstring>
#include <fstream>

#include <jni.h>
#include <android/log.h>

#include "lzw/compressor/include/compressor.hpp"
#include "lzw/compressor/include/decompressor.hpp"
#include "lzw/encryptor/chacha20.hpp"
#include "lzw/utils/utils.hpp"

#define LOG_TAG "MyNativeApp"
#define LOGD(...) __android_log_print(ANDROID_LOG_DEBUG, LOG_TAG, __VA_ARGS__)

extern "C" JNIEXPORT jstring JNICALL
Java_com_example_zipper_MainActivity_archiveAndSecure(
        JNIEnv* env,
        jobject /* this */,
        jstring read_from, jstring save_to) {
    LOGD(__PRETTY_FUNCTION__);

    const char *native_read_from = env->GetStringUTFChars(read_from, 0);
    const char *native_save_to = env->GetStringUTFChars(save_to, 0);

    auto text = utils::read_from_file_into_string(native_read_from);
    my_compress::compressor_t compressor{};
    auto res = compressor.compress(text.begin(), text.end());

    std::array<uint32_t, 32> key;
    for(int i=0; i<32; i++) {
        key[i] = i;
    }
    std::array<uint32_t, 12> nonce = {0,0,0,0,0,0,0,0x4a,0,0,0,0};
    
    my_compress::chacha20_t chacha(key, nonce, 1);
    std::vector<uint32_t> ciphertext;
    std::vector<uint32_t> plaintext(text.begin(), text.end());
    chacha.crypt(plaintext, ciphertext);

    try {
        utils::write_int_data_into_bin_file(ciphertext, native_save_to);
    } catch(std::string& err) {
        LOGD("%s", err.c_str());
    }

    return env->NewStringUTF(native_save_to);
}

extern "C" JNIEXPORT jstring JNICALL
Java_com_example_zipper_MainActivity_unarchiveAndOpen(
        JNIEnv* env,
        jobject /* this */,
        jstring read_from, jstring save_to) {

    LOGD(__PRETTY_FUNCTION__);

    const char *native_path_to = env->GetStringUTFChars(read_from, 0);
    std::vector<uint32_t> data{};
    try {
        data = utils::read_from_bin_file_into_int(native_path_to);
    } catch(std::string& err) {
        LOGD("%s", err.c_str());
    }

    std::array<uint32_t, 32> key;
    for(int i=0; i<32; i++) {
        key[i] = i;
    }
    std::array<uint32_t, 12> nonce = {0,0,0,0,0,0,0,0x4a,0,0,0,0};

    my_compress::chacha20_t chacha(key, nonce, 1);
    std::vector<uint32_t> deciphertext;
    std::vector<uint32_t> plaintext(data.begin(), data.end());
    chacha.crypt(plaintext, deciphertext);

    my_compress::decompressor_t decompressor{};
    auto res = decompressor.decompress(deciphertext.begin(), deciphertext.end());

    return env->NewStringUTF(res.c_str());
}