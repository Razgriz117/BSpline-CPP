
add_executable(test_bspline test_bspline.cpp)
target_link_libraries(test_bspline
    PRIVATE GTest::gtest GTest::gtest_main bspline_lib)
add_test(NAME BSplineTests COMMAND test_bspline)

add_executable(test_utils test_utils.cpp)
target_link_libraries(test_utils
    PRIVATE GTest::gtest GTest::gtest_main utils_lib)
add_test(NAME UtilsTests COMMAND test_utils)

add_executable(test_tise test_tise.cpp)
target_link_libraries(test_tise
    PRIVATE GTest::gtest GTest::gtest_main tise_lib)
add_test(NAME TISETests COMMAND test_tise)

add_executable(test_time_evolution test_time_evolution.cpp)
target_include_directories(test_time_evolution PRIVATE /usr/include/eigen3)
target_link_libraries(test_time_evolution
    PRIVATE GTest::gtest GTest::gtest_main time_evolution_lib)
add_test(NAME TimeEvolutionTests COMMAND test_time_evolution)
