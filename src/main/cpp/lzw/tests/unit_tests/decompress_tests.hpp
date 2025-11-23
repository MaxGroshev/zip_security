#pragma once

using namespace my_compress;

//-----------------------------------------------------------------------------------------

class lzw_test2 : public ::testing::Test {
    protected:
        std::vector<int> decompress_data1  = {66, 65, 128, 129, 65, 132};
        std::string correct_res1 = {"BABAABAAA"};

        std::vector<int> decompress_data2 = {84, 79, 66, 69, 79, 82, 78, 79, 84,
                                         128, 130, 132, 137, 131, 133, 135};
        std::string correct_res2 = {"TOBEORNOTTOBEORTOBEORNOT"};

    void SetUp() {
    }
};

//-----------------------------------------------------------------------------------------

TEST_F(lzw_test2, decompress_test1) {

    lzw_t lzw{set_up_t{}, 127, 127};
    auto res = lzw.decompress(decompress_data1);

    ASSERT_TRUE(res == correct_res1);
}

TEST_F(lzw_test2, decompress_test2) {

    lzw_t lzw{set_up_t{}, 127, 127};
    auto res = lzw.decompress(decompress_data2);

    ASSERT_TRUE(res == correct_res2);
}
