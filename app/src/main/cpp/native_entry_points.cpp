#include "native_entry_points.hpp"

extern "C" JNIEXPORT jstring JNICALL
Java_nativecpp_LzwArchiver_archiveAndSecure(
        JNIEnv* env,
        jobject /* this */,
        jstring read_from, jstring save_to) {
    LOGD(__PRETTY_FUNCTION__);

    const char *native_read_from = env->GetStringUTFChars(read_from, 0);
    const char *native_save_to = env->GetStringUTFChars(save_to, 0);

    auto text = utils::read_from_file_into_string(native_read_from);
    my_compress::compressor_t compressor{};
    auto compressed_data = compressor.compress(text.begin(), text.end());

    my_compress::key_generator_t<uint32_t> gn{};
    auto key = gn.generate();
    std::array<uint32_t, 12> nonce = {0,0,0,0,0,0,0,0x4a,0,0,0,0};
    
    my_compress::chacha20_t chacha(key, nonce, 1);
    std::vector<uint32_t> ciphertext;
    std::vector<uint32_t> uint_vector(compressed_data.begin(), compressed_data.end());
    chacha.crypt(uint_vector, ciphertext);

    try {
        utils::write_int_data_into_bin_file(ciphertext, native_save_to);
    } catch(std::runtime_error &err) {
        LOGE("%s", err.what());
    }

    auto str_key = gn.uint8_vector_to_hex_string(key);
    LOGD("%s", str_key.c_str());
    return env->NewStringUTF(str_key.c_str());
}

extern "C" JNIEXPORT jstring JNICALL
Java_nativecpp_LzwArchiver_unarchiveAndOpen(
        JNIEnv* env,
        jobject /* this */,
        jstring read_from, jstring save_to, jstring my_key) {

    LOGD(__PRETTY_FUNCTION__);

    const char *native_read_from = env->GetStringUTFChars(read_from, 0);
    const char *native_save_to = env->GetStringUTFChars(save_to, 0);

    const char *native_key = env->GetStringUTFChars(my_key, 0);
    std::string str_key(native_key);
    env->ReleaseStringUTFChars(my_key, native_key);
    LOGD("%s", str_key.c_str());

    std::vector<uint32_t> data{};
    try {
        data = utils::read_from_bin_file_into_int(native_read_from);
    } catch(std::runtime_error &err) {
        LOGE("%s", err.what());
    }
    my_compress::key_generator_t<uint32_t> gn {};

    std::array<uint32_t, 32> key = gn.hex_string_to_uint8_vector(str_key);
    std::array<uint32_t, 12> nonce = {0,0,0,0,0,0,0,0x4a,0,0,0,0};

    my_compress::chacha20_t chacha(key, nonce, 1);
    std::vector<uint32_t> deciphertext;
    std::vector<uint32_t> plaintext(data.begin(), data.end());
    chacha.crypt(plaintext, deciphertext);

    my_compress::decompressor_t decompressor{};
    auto res = decompressor.decompress(deciphertext.begin(), deciphertext.end());
    utils::write_string_into_file(native_save_to, res);

    return env->NewStringUTF(res.c_str());
}