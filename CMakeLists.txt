file(READ Makefile _MAKEFILE_CONTENTS)
string(REGEX MATCH "KFVER=([0-9]+)\n" KFVER_MATCH "${_MAKEFILE_CONTENTS}")
if(NOT KFVER_MATCH)
    message(FATAL_ERROR "Cannot extract version from Makefile")
endif()
set(MAKEFILE_EXTRACTED_VERSION "${CMAKE_MATCH_1}")

cmake_minimum_required(VERSION 3.3)
project(kissfft VERSION "${MAKEFILE_EXTRACTED_VERSION}")

option(KISSFFT_SIMD "Build kissfft with SIMD")
option(KISSFFT_FLOAT "Build kissfft with float type" ON)

set(KISSFFT_FIXED_WIDTH_VALID 16 32)
set(KISSFFT_FIXED_WIDTH "16" CACHE STRING "Width of kissfft fixed integer (16/32)")
set_property(CACHE KISSFFT_FIXED_WIDTH PROPERTY STRINGS ${KISSFFT_FIXED_WIDTH_VALID})
if(NOT KISSFFT_FIXED_WIDTH IN_LIST KISSFFT_FIXED_WIDTH_VALID)
    message(FATAL_ERROR "KISSFFT_FIXED_WIDTH must be 16 or 32")
endif()

set(KISSFFT_FLOAT_TYPE_VALID float double)
set(KISSFFT_FLOAT_TYPE "float" CACHE STRING "Floating type of kissfft (float/double)")
set_property(CACHE KISSFFT_FLOAT_TYPE PROPERTY STRINGS ${KISSFFT_FLOAT_TYPE_VALID})
if(NOT KISSFFT_FLOAT_TYPE IN_LIST KISSFFT_FLOAT_TYPE_VALID)
    message(FATAL_ERROR "KISSFFT_FLOAT_TYPE must be float or double")
endif()

add_library(kissfft
  kiss_fft.c)
add_library(kissfft::kissfft ALIAS kissfft)

target_include_directories(kissfft PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:.>)

set(KISSFFT_SUFFIX)
if(KISSFFT_SIMD)
    target_compile_definitions(kissfft PUBLIC USE_SIMD)
    set(KISSFFT_SUFFIX "${KISSFFT_SUFFIX}_simd")
endif()
if(KISSFFT_FLOAT)
    target_compile_definitions(kissfft PUBLIC kiss_fft_scalar=${KISSFFT_FLOAT_TYPE})
    set(KISSFFT_SUFFIX "${KISSFFT_SUFFIX}_${KISSFFT_FLOAT_TYPE}")
else()
    target_compile_definitions(kissfft PUBLIC KISSFFT_FIXED_POINT=${KISSFFT_FIXED_WIDTH})
    set(KISSFFT_SUFFIX "${KISSFFT_SUFFIX}_int${KISSFFT_FIXED_WIDTH}")
endif()

set_target_properties(kissfft PROPERTIES
    OUTPUT_NAME "kissfft${KISSFFT_SUFFIX}"
    DEFINE_SYMBOL KISS_FFT_BUILD
    EXPORT_NAME "kissfit${KISSFFT_SUFFIX}"
    SOVERSION ${PROJECT_VERSION})
if(BUILD_SHARED_LIBS)
    target_compile_definitions(kissfft PUBLIC KISS_FFT_SHARED)
    set_target_properties(kissfft PROPERTIES
        C_VISIBILITY_PRESET hidden)
endif()

option(KISSFFT_INSTALL "Enable kissfft install" ON)
if (KISSFFT_INSTALL)
    include(GNUInstallDirs)
    install(TARGETS kissfft EXPORT kissfft
        ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
    install(FILES "kiss_fft.h" DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")

    set(KISSFFT_INSTALL_CMAKE "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}" CACHE FILEPATH "Install destination of kissfft cmake modules")
    mark_as_advanced(KISSFFT_INSTALL_CMAKE)

    install(EXPORT kissfft DESTINATION "${KISSFFT_INSTALL_CMAKE}"
        NAMESPACE "kissfft::"
        FILE "${PROJECT_NAME}${KISSFFT_SUFFIX}-targets.cmake")
    include(CMakePackageConfigHelpers)
    configure_package_config_file(kissfft-config.cmake.in kissfft-config.cmake
        INSTALL_DESTINATION "${KISSFFT_INSTALL_CMAKE}")
    write_basic_package_version_file(kissfft-config-version.cmake COMPATIBILITY AnyNewerVersion)
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/kissfft-config.cmake" "${CMAKE_CURRENT_BINARY_DIR}/kissfft-config-version.cmake"
        DESTINATION "${KISSFFT_INSTALL_CMAKE}")
endif()
