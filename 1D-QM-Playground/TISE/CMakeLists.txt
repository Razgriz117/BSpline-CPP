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

# B-spline library
add_library(bspline_lib BSpline.cpp)
target_include_directories(bspline_lib
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)
target_link_libraries(bspline_lib
    PUBLIC ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}
)

# Utils library
add_library(utils_lib utils/utils.cpp)
target_include_directories(utils_lib
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/utils
        ${CMAKE_CURRENT_SOURCE_DIR}
)
target_link_libraries(utils_lib PUBLIC bspline_lib)

# TISE solver library (Project Part 1)
add_library(tise_lib tise.cpp)
target_include_directories(tise_lib
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)
target_link_libraries(tise_lib
    PUBLIC bspline_lib ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}
)

# Time evolution library (Project Part 2)
add_library(time_evolution_lib time_evolution.cpp)
target_include_directories(time_evolution_lib
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        /usr/include/eigen3
)
target_link_libraries(time_evolution_lib PUBLIC bspline_lib tise_lib)

# Main executable
add_executable(H-BoundStates main.cpp)

target_compile_features(H-BoundStates PRIVATE cxx_std_17)

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
    tise_lib
    time_evolution_lib
    utils_lib
    ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}
)
