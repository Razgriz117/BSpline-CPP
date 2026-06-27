
add_executable(test_bspline test_bspline.cpp)

target_link_libraries(test_bspline
    PRIVATE
        GTest::gtest
        GTest::gtest_main
        bspline_lib
)

add_test(NAME BSplineTests COMMAND test_bspline)