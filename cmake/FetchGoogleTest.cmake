# Finds GoogleTest on the system, or fetches and builds it from source if
# not found and fetching is allowed. Shared between 1D-QM-Playground/TISE/
# CMakeLists.txt and student-base-code/CMakeLists.txt (PR #18 review,
# finding 4.3): the two copies of this block had already diverged in their
# fetch-gate option name within the same PR that introduced the
# duplication -- factored out here so a future GTest version bump or
# CRT-linkage fix can't be applied to one copy and forgotten in the other.
#
# Call as:
#   tise_fetch_or_find_gtest(ALLOW_FETCH ${SOME_PROJECT_OWN_FETCH_OPTION})
#
# where SOME_PROJECT_OWN_FETCH_OPTION is whichever option() the calling
# CMakeLists.txt already declares for this purpose (its name is left
# untouched -- only this shared logic, not either project's public flag,
# is being deduplicated). The caller must `include(FetchContent)` and
# `enable_testing()` itself before calling this, and `add_subdirectory(tests)`
# itself afterward -- this macro only resolves the GTest::gtest/
# GTest::gtest_main targets, it does not know about either project's own
# test directory.
macro(tise_fetch_or_find_gtest)
    cmake_parse_arguments(TISE_GTEST "" "ALLOW_FETCH" "" ${ARGN})

    find_package(GTest QUIET)
    if(NOT GTest_FOUND)
        if(NOT TISE_GTEST_ALLOW_FETCH)
            message(FATAL_ERROR
                "GTest not found and fetching dependencies is disabled. "
                "Install it (apt: libgtest-dev, brew: googletest) or allow "
                "fetching.")
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
endmacro()
