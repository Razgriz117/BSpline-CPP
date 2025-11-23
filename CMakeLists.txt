project(H-BoundStates LANGUAGES CXX)

cmake_minimum_required(VERSION 3.10)

# Set C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Find BLAS and LAPACK
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Create executable
add_executable(H-BoundStates main.cpp BSpline.cpp)

target_compile_features(H-BoundStates PRIVATE 
     cxx_std_17
)

target_include_directories(H-BoundStates PRIVATE 
     /usr/include/eigen3
     ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(H-BoundStates PRIVATE 
     ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}
)