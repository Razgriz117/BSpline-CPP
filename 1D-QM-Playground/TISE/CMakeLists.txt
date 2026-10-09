# ---------------------------------------------------------------------------
# Portability contract (see ../README.md "Build"):
#
#   A student with a C++17 compiler, CMake >= 3.16, and a network connection
#   can build this on macOS, Linux or Windows with no package manager, no
#   admin rights, and no platform-specific instructions.
#
# Everything below exists to keep that true. Each dependency is looked for on
# the system first (fast, and what CI and developer machines use) and, only if
# missing, fetched and built from source. Every one of them is either
# header-only or a small, self-contained CMake project.
#
# Notably absent: LAPACK/BLAS. The solver used to call LAPACK's banded dsbgv
# for the generalized eigenproblem, but LAPACK is Fortran -- there is no
# build-from-source path for it that does not drag in a Fortran toolchain,
# which on Windows is exactly the barrier this file exists to remove. That one
# call now goes through Eigen instead (see tise.cpp's
# solveGeneralizedEigenproblem); at the basis sizes this code runs at the two
# are indistinguishable in speed.
# ---------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.16)

project(H-BoundStates LANGUAGES CXX)

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

# ===========================================================================
# Dependency resolution
# ===========================================================================
# Fetched dependencies live under the build tree (FetchContent's default).
# Relocating them to a shared directory in the source tree looks tempting --
# it would survive `rm -rf build` -- but FETCHCONTENT_BASE_DIR holds each
# dependency's *build* directory as well as its sources, so two build trees
# configured from this project would then share and overwrite one another's
# dependency builds.
include(FetchContent)

option(TISE_FETCH_DEPENDENCIES
       "Download and build any dependency not found on the system." ON)

# Several of the fetched projects still declare cmake_minimum_required() with a
# version CMake 4.x refuses outright. This is the documented escape hatch and
# affects only those subprojects' policy defaults, not ours.
set(CMAKE_POLICY_VERSION_MINIMUM 3.5)

# Fetched dependencies are built static so the resulting binaries run from
# anywhere without a library search path -- students copy build/tise_solver
# around, and a dangling .so/.dylib/.dll reference is a confusing failure. Only
# applies when the user has not made this choice themselves.
# Set as a cache entry, not a plain variable: muparser still calls
# option(BUILD_SHARED_LIBS ...) under CMP0077 OLD, which clears a plain
# variable of that name and would give us a shared build anyway.
if(NOT DEFINED BUILD_SHARED_LIBS)
    set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build fetched dependencies shared.")
endif()

# --- Eigen (header-only) ---------------------------------------------------
# The generalized eigensolver behind tise_lib, and the linear algebra behind
# time_evolution_lib.
find_package(Eigen3 3.3 QUIET NO_MODULE)
if(NOT Eigen3_FOUND)
    if(NOT TISE_FETCH_DEPENDENCIES)
        message(FATAL_ERROR
            "Eigen3 not found and TISE_FETCH_DEPENDENCIES=OFF. Install it "
            "(apt: libeigen3-dev, brew: eigen) or allow fetching.")
    endif()
    message(STATUS "Eigen3 not found on system -- fetching source")
    # SOURCE_SUBDIR names a directory with no CMakeLists.txt, which tells
    # FetchContent_MakeAvailable to download without running Eigen's own build
    # system. Eigen is header-only, so an INTERFACE target over the source tree
    # is all that is needed, and this skips Eigen's heavy (and CMake-4-hostile)
    # configure entirely.
    # The release tarball rather than a git clone: even a shallow clone of
    # Eigen transfers the whole ~90 MB working tree, where the archive is
    # 2.6 MB. The hash pins exactly what gets built.
    FetchContent_Declare(eigen
        URL      https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
        URL_HASH SHA256=8586084f71f9bde545ee7fa6d00288b264a2b7ac3607b974e54d13e7162c1c72
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
        SOURCE_SUBDIR  no-such-dir)
    FetchContent_MakeAvailable(eigen)
    add_library(eigen_fetched INTERFACE)
    target_include_directories(eigen_fetched SYSTEM INTERFACE ${eigen_SOURCE_DIR})
    add_library(Eigen3::Eigen ALIAS eigen_fetched)
endif()

# --- nlohmann_json (header-only) -------------------------------------------
# Parses the JSON-encoded piecewise potential and writes the warnings.json
# sidecar.
find_package(nlohmann_json 3.10.5 QUIET)
if(NOT nlohmann_json_FOUND)
    if(NOT TISE_FETCH_DEPENDENCIES)
        message(FATAL_ERROR
            "nlohmann_json not found and TISE_FETCH_DEPENDENCIES=OFF. Install "
            "it (apt: nlohmann-json3-dev, brew: nlohmann-json) or allow fetching.")
    endif()
    message(STATUS "nlohmann_json not found on system -- fetching source")
    set(JSON_BuildTests OFF CACHE INTERNAL "")
    FetchContent_Declare(nlohmann_json
        GIT_REPOSITORY https://github.com/nlohmann/json.git
        GIT_TAG        v3.11.3
        GIT_SHALLOW    TRUE)
    FetchContent_MakeAvailable(nlohmann_json)
endif()

# --- muparser --------------------------------------------------------------
# Evaluates the per-piece potential expression strings at runtime.
#
# Resolution order matters here. The previous build used pkg-config only, which
# is why this project could not be built on Windows at all: pkg-config is not
# part of a normal MSVC toolchain. CMake config mode is tried first (muparser
# >= 2.3 installs muparserConfig.cmake, including via Homebrew and vcpkg),
# pkg-config second for older Linux installs, source last.
set(TISE_MUPARSER_TARGET "")
find_package(muparser CONFIG QUIET)
if(muparser_FOUND)
    set(TISE_MUPARSER_TARGET muparser::muparser)
    message(STATUS "muparser: found via CMake config")
else()
    find_package(PkgConfig QUIET)
    if(PkgConfig_FOUND)
        pkg_check_modules(MUPARSER QUIET IMPORTED_TARGET muparser)
    endif()
    if(MUPARSER_FOUND)
        set(TISE_MUPARSER_TARGET PkgConfig::MUPARSER)
        message(STATUS "muparser: found via pkg-config")
    else()
        if(NOT TISE_FETCH_DEPENDENCIES)
            message(FATAL_ERROR
                "muparser not found and TISE_FETCH_DEPENDENCIES=OFF. Install it "
                "(apt: libmuparser-dev, brew: muparser) or allow fetching.")
        endif()
        message(STATUS "muparser not found on system -- fetching source")
        set(ENABLE_SAMPLES OFF CACHE INTERNAL "")
        set(ENABLE_OPENMP  OFF CACHE INTERNAL "")
        FetchContent_Declare(muparser
            GIT_REPOSITORY https://github.com/beltoforion/muparser.git
            GIT_TAG        v2.3.4
            GIT_SHALLOW    TRUE)
        # muparser registers its own ParserTest with CTest whenever
        # BUILD_TESTING is on, which would put a dependency's test in the
        # middle of `ctest` output for this project and make a student think
        # it was one of ours. The cache entry has to be overridden rather
        # than shadowed by a normal variable: muparser calls include(CTest),
        # whose own option(BUILD_TESTING ...) runs under CMP0077 OLD and
        # clears any normal variable of that name before muparser tests it.
        set(TISE_BUILD_TESTING_SAVED "${BUILD_TESTING}")
        set(BUILD_TESTING OFF CACHE BOOL "Build the testing tree." FORCE)
        FetchContent_MakeAvailable(muparser)
        set(BUILD_TESTING "${TISE_BUILD_TESTING_SAVED}"
            CACHE BOOL "Build the testing tree." FORCE)
        set(TISE_MUPARSER_TARGET muparser)
    endif()
endif()

# --- yaml-cpp --------------------------------------------------------------
# Parses config.yaml. Only tise_solver needs it.
option(BUILD_TISE_SOLVER "Build the config-driven tise_solver binary." ON)
if(BUILD_TISE_SOLVER)
    set(TISE_YAMLCPP_TARGET "")
    find_package(yaml-cpp QUIET)
    if(yaml-cpp_FOUND)
        # 0.8 exports yaml-cpp::yaml-cpp; older packages export a bare yaml-cpp.
        if(TARGET yaml-cpp::yaml-cpp)
            set(TISE_YAMLCPP_TARGET yaml-cpp::yaml-cpp)
        else()
            set(TISE_YAMLCPP_TARGET yaml-cpp)
        endif()
        message(STATUS "yaml-cpp: found on system")
    else()
        if(NOT TISE_FETCH_DEPENDENCIES)
            message(FATAL_ERROR
                "yaml-cpp not found and TISE_FETCH_DEPENDENCIES=OFF. Install it "
                "(apt: libyaml-cpp-dev, brew: yaml-cpp) or allow fetching. "
                "Pass -DBUILD_TISE_SOLVER=OFF to skip the binary that needs it.")
        endif()
        message(STATUS "yaml-cpp not found on system -- fetching source")
        set(YAML_CPP_BUILD_TESTS    OFF CACHE INTERNAL "")
        set(YAML_CPP_BUILD_TOOLS    OFF CACHE INTERNAL "")
        set(YAML_CPP_BUILD_CONTRIB  OFF CACHE INTERNAL "")
        set(YAML_CPP_INSTALL        OFF CACHE INTERNAL "")
        FetchContent_Declare(yaml-cpp
            GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
            GIT_TAG        0.8.0
            GIT_SHALLOW    TRUE)
        FetchContent_MakeAvailable(yaml-cpp)
        set(TISE_YAMLCPP_TARGET yaml-cpp::yaml-cpp)
    endif()
endif()

# --- GoogleTest ------------------------------------------------------------
if(BUILD_TESTING)
    enable_testing()
    find_package(GTest QUIET)
    if(NOT GTest_FOUND)
        if(NOT TISE_FETCH_DEPENDENCIES)
            message(FATAL_ERROR
                "GTest not found and TISE_FETCH_DEPENDENCIES=OFF. Install it "
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

# ===========================================================================
# Targets
# ===========================================================================

# B-spline library
add_library(bspline_lib BSpline.cpp)
target_include_directories(bspline_lib
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
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
    bspline_lib
    ${TISE_MUPARSER_TARGET}
    Eigen3::Eigen
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
)

# TISE solver CLI (Controller<->TISE contract, docs/SDD.md 7.2.1).
# BUILD_TISE_SOLVER=OFF is the way to opt out of tise_solver deliberately; the
# option itself is declared above, next to the yaml-cpp lookup it gates.
if(BUILD_TISE_SOLVER)
    add_executable(tise_solver tise_solver_main.cpp)
    target_compile_features(tise_solver PRIVATE cxx_std_17)
    target_link_libraries(tise_solver PRIVATE
        tise_lib
        ${TISE_YAMLCPP_TARGET}
        nlohmann_json::nlohmann_json
    )
endif()
