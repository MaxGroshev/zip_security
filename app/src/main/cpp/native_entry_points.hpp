#pragma once

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
#include "lzw/encryptor/key_generator.hpp"
#include "lzw/utils/utils.hpp"

#define LOG_TAG "MyNativeApp"
#define LOGD(...) __android_log_print(ANDROID_LOG_DEBUG, LOG_TAG, __VA_ARGS__)
#define LOGE(...) __android_log_print(ANDROID_LOG_ERROR, LOG_TAG, __VA_ARGS__)
#define LOGI(...) __android_log_print(ANDROID_LOG_INFO, LOG_TAG, __VA_ARGS__)
