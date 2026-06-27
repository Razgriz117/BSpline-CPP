
add_executable(test_bspline test_bspline.cpp)
target_link_libraries(test_bspline
    PRIVATE
        GTest::gtest
        GTest::gtest_main
        bspline_lib
)
add_test(NAME BSplineTests COMMAND test_bspline)

add_executable(test_utils test_utils.cpp)
target_link_libraries(test_utils
    PRIVATE
        GTest::gtest
        GTest::gtest_main
        utils_lib
)
add_test(NAME UtilsTests COMMAND test_utils)
