cmake_minimum_required(VERSION 3.16)

project(H-BoundStates LANGUAGES CXX)

# Set C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Options
option(BUILD_TESTING "Build the testing tree." OFF)

# Otherwise a plain `cmake -S . -B build` compiles with no optimization at all.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

include(FetchContent)
option(FETCH_DEPENDENCIES
       "Download and build any dependency not found on the system." ON)
set(CMAKE_POLICY_VERSION_MINIMUM 3.5)

# --- BLAS/LAPACK ------------------------------------------------------------
# Required, and deliberately so: main.cpp calls LAPACK's DSBGV directly, and
# the banded-storage layout it sets up is part of what this project teaches.
# That is also why this dependency cannot be fetched and built like the others
# -- LAPACK is Fortran, so building it from source needs a Fortran compiler.
#
# (The full solver in ../1D-QM-Playground/ made the opposite trade and went
# through Eigen instead, so that it builds with no Fortran toolchain at all.)
find_package(BLAS QUIET)
find_package(LAPACK QUIET)
if(NOT LAPACK_FOUND)
    message(FATAL_ERROR
        "LAPACK/BLAS not found, and they cannot be built from source without a "
        "Fortran compiler.\n"
        "  Linux  : sudo apt-get install liblapack-dev libblas-dev\n"
        "  macOS  : already present as Apple's Accelerate framework; if this "
        "still fails, install Xcode command line tools "
        "(xcode-select --install)\n"
        "  Windows: install via vcpkg (vcpkg install lapack, then configure "
        "with -DCMAKE_TOOLCHAIN_FILE=<vcpkg>/scripts/buildsystems/vcpkg.cmake), "
        "or build inside WSL, or use a conda environment providing "
        "liblapack.\n")
endif()

# --- GoogleTest -------------------------------------------------------------
if(BUILD_TESTING)
    enable_testing()
    find_package(GTest QUIET)
    if(NOT GTest_FOUND)
        if(NOT FETCH_DEPENDENCIES)
            message(FATAL_ERROR
                "GTest not found and FETCH_DEPENDENCIES=OFF. Install it "
                "(apt: libgtest-dev, brew: googletest) or allow fetching.")
        endif()
        message(STATUS "GTest not found on system -- fetching source")
        # Without this, GoogleTest links the static CRT on MSVC while this
        # project links the dynamic one, and the test binaries fail to link.
        set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
        set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
        FetchContent_Declare(googletest
            GIT_REPOSITORY https://github.com/google/googletest.git
            GIT_TAG        v1.14.0
            GIT_SHALLOW    TRUE)
        FetchContent_MakeAvailable(googletest)
    endif()
    add_subdirectory(tests)
endif()

# Create BSpline library
add_library(bspline_lib BSpline.cpp)
target_include_directories(bspline_lib
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
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
