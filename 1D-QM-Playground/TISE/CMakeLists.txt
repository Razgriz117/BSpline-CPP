project(H-BoundStates LANGUAGES CXX)

cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(BUILD_TESTING "Build the testing tree." OFF)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

if (BUILD_TESTING)
    enable_testing()
    find_package(GTest REQUIRED)
    add_subdirectory(tests)
endif()

# B-spline library (used by both the main executable and tests)
add_library(bspline_lib BSpline.cpp)
target_include_directories(bspline_lib
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)
target_link_libraries(bspline_lib
    PUBLIC
        ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}
)

# Utils library (used by the main executable and tests)
add_library(utils_lib utils/utils.cpp)
target_include_directories(utils_lib
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/utils
        ${CMAKE_CURRENT_SOURCE_DIR}
)
target_link_libraries(utils_lib
    PUBLIC
        bspline_lib
)

# Main executable
add_executable(H-BoundStates main.cpp)

target_compile_features(H-BoundStates PRIVATE cxx_std_17)

# Define output directory for timesteps.
set(TIMESTEPS_OUTPUT_DIR "${CMAKE_BINARY_DIR}/timesteps")
add_custom_target(create_timesteps_dir ALL
    COMMAND ${CMAKE_COMMAND} -E make_directory ${TIMESTEPS_OUTPUT_DIR}
)

target_compile_definitions(H-BoundStates PRIVATE
    TIMESTEPS_DIR="${TIMESTEPS_OUTPUT_DIR}"
)

target_include_directories(H-BoundStates PRIVATE
    /usr/include/eigen3
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(H-BoundStates PRIVATE
    bspline_lib
    utils_lib
    ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}
)
