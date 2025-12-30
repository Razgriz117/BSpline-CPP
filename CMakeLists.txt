project(H-BoundStates LANGUAGES CXX)

cmake_minimum_required(VERSION 3.10)

# Set C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Options
option(BUILD_TESTING "Build the testing tree." OFF)

# Find BLAS and LAPACK
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

if (BUILD_TESTING)
     enable_testing()
     find_package(GTest REQUIRED)
     add_subdirectory(tests)
endif()

# Create BSpline library
add_library(bspline_lib BSpline.cpp)
target_include_directories(bspline_lib
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        /usr/include/eigen3
)
target_link_libraries(bspline_lib
    PUBLIC
        ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}
)

# Create executable
add_executable(H-BoundStates main.cpp)

target_compile_features(H-BoundStates PRIVATE 
     cxx_std_17
)

target_link_libraries(H-BoundStates PRIVATE bspline_lib)
