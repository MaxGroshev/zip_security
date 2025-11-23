#pragma once

using namespace my_compress;

//-----------------------------------------------------------------------------------------

class lzw_test : public ::testing::Test {
    protected:
        std::string compress_data1 = {"BABAABAAA"};
        std::vector<int> correct_res1  = {'B', 'A', 128, 129, 'A', 132};

        std::string compress_data2 = {"TOBEORNOTTOBEORTOBEORNOT"};
        std::vector<int> correct_res2 = {84, 79, 66, 69, 79, 82, 78, 79, 84,
                                         128, 130, 132, 137, 131, 133, 135};

    void SetUp() {
    }
};

//-----------------------------------------------------------------------------------------

TEST_F(lzw_test, compress_test1) {

    lzw_t lzw{set_up_t{}, 127, 127};
    auto  res = lzw.compress(compress_data1);

    ASSERT_TRUE(res == correct_res1);
}

TEST_F(lzw_test, compress_test2) {

    lzw_t lzw{set_up_t{}, 127, 127};
    auto res = lzw.compress(compress_data2);

    ASSERT_TRUE(res == correct_res2);
}
