project(H-BoundStates LANGUAGES CXX)

cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(BUILD_TESTING "Build the testing tree." OFF)

# Without this a bare `cmake -S TISE -B TISE/build` builds with NO optimization
# flags at all, which for a continuum sweep (100+ energies, each a Coulomb-wave
# evaluation) is the difference between seconds and minutes. Only applies when
# the user has not chosen a build type themselves.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(ENABLE_COVERAGE "Build with --coverage instrumentation (gcov/gcovr)" OFF)
if(ENABLE_COVERAGE)
    add_compile_options(--coverage -O0 -g)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} --coverage")
endif()

find_package(BLAS REQUIRED)
# Eigen3 via its own CMake config rather than a hardcoded /usr/include/eigen3,
# which is Debian/Ubuntu-only -- Homebrew puts it under /opt/homebrew/include,
# MacPorts under /opt/local/include, conda under $CONDA_PREFIX/include. The
# imported target carries the right include directory on every platform.
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(LAPACK REQUIRED)
# nlohmann_json: parses the JSON-encoded piecewise potential passed as input.
find_package(nlohmann_json 3.10.5 REQUIRED)
# muparser: evaluates the per-piece potential expression strings at runtime.
find_package(PkgConfig REQUIRED)
pkg_check_modules(MUPARSER REQUIRED IMPORTED_TARGET muparser)
add_library(muparser::muparser ALIAS PkgConfig::MUPARSER)

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
target_link_libraries(tise_lib PUBLIC 
    bspline_lib ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}
    muparser::muparser
)

# Time evolution library (Project Part 2)
add_library(time_evolution_lib time_evolution.cpp)
target_include_directories(time_evolution_lib
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)
target_link_libraries(time_evolution_lib PUBLIC bspline_lib tise_lib Eigen3::Eigen)

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
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(H-BoundStates PRIVATE
    Eigen3::Eigen
    tise_lib
    time_evolution_lib
    utils_lib
    nlohmann_json::nlohmann_json
    ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}
)

# TISE solver CLI (Controller<->TISE contract, docs/SDD.md §7.2.1).
# BUILD_TISE_SOLVER=OFF is the way to opt out of tise_solver deliberately.
# Left ON (the default), a missing yaml-cpp is a hard error: this is the binary
# controller.py drives and the one students actually run, and the previous
# message(WARNING) let the build "succeed" without it, turning a missing
# dependency into a much later, much less obvious "could not execute" from the
# controller.
option(BUILD_TISE_SOLVER "Build the config-driven tise_solver binary." ON)
if (BUILD_TISE_SOLVER)
    find_package(yaml-cpp REQUIRED)
    add_executable(tise_solver tise_solver_main.cpp)
    target_compile_features(tise_solver PRIVATE cxx_std_17)
    target_link_libraries(tise_solver PRIVATE
        tise_lib
        yaml-cpp
        nlohmann_json::nlohmann_json
    )
endif()