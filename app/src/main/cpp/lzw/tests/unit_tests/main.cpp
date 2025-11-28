#include "./unit_tests.hpp"

//-----------------------------------------------------------------------------------------

int main (int argc, char* argv[]) {

    testing::InitGoogleTest(&argc, argv);
    int ret_val = RUN_ALL_TESTS();

    return 0;
}
