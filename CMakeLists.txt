file(READ Makefile _MAKEFILE_CONTENTS)
string(REGEX MATCH "KFVER=([0-9]+)\n" KFVER_MATCH "${_MAKEFILE_CONTENTS}")
if(NOT KFVER_MATCH)
    message(FATAL_ERROR "Cannot extract version from Makefile")
endif()
set(MAKEFILE_EXTRACTED_VERSION "${CMAKE_MATCH_1}")

cmake_minimum_required(VERSION 3.13)
project(kissfft VERSION "${MAKEFILE_EXTRACTED_VERSION}")

option(KISSFFT_SIMD "Build kissfft with SIMD" OFF)
option(KISSFFT_FLOAT "Build kissfft with float type" ON)
option(KISSFFT_OPENMP "Build kissfft with openmp" OFF)

if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    add_compile_options(-ffast-math -fomit-frame-pointer
        -W -Wall -Waggregate-return -Wcast-align -Wcast-qual -Wshadow -Wwrite-strings
        "$<$<COMPILE_LANGUAGE:C>:-Wstrict-prototypes;-Wmissing-prototypes;-Wnested-externs;-Wbad-function-cast>")
endif()

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

target_include_directories(kissfft PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

set(KISSFFT_COMPILE_DEFINITIONS)
if(KISSFFT_SIMD)
    list(APPEND KISSFFT_COMPILE_DEFINITIONS USE_SIMD)
    set(KISSFFT_DATATYPE "simd")
else()
    if(KISSFFT_FLOAT)
        list(APPEND KISSFFT_COMPILE_DEFINITIONS kiss_fft_scalar=${KISSFFT_FLOAT_TYPE})
        set(KISSFFT_DATATYPE "${KISSFFT_FLOAT_TYPE}")
    else()
        list(APPEND KISSFFT_COMPILE_DEFINITIONS KISSFFT_FIXED_POINT=${KISSFFT_FIXED_WIDTH})
        set(KISSFFT_DATATYPE "int${KISSFFT_FIXED_WIDTH}")
    endif()
endif()
if(KISSFFT_OPENMP)
    if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
        target_compile_options(kissfft PRIVATE -fopenmp)
        target_link_options(kissfft PUBLIC -fopenmp)
    else()
        message(FATAL_ERROR "Don't know how to enable openmp for this compiler")
    endif()
endif()
message(STATUS "Building KissFFT with datatype=${KISSFFT_DATATYPE}")

if(BUILD_SHARED_LIBS)
    list(APPEND KISSFFT_COMPILE_DEFINITIONS KISS_FFT_SHARED)
    set_target_properties(kissfft PROPERTIES
        C_VISIBILITY_PRESET hidden)
    set(KISSFFT_EXPORT_SUFFIX "-shared")
else()
    set(KISSFFT_EXPORT_SUFFIX "-static")
endif()
target_compile_definitions(kissfft PUBLIC ${KISSFFT_COMPILE_DEFINITIONS})
set(KISSFFT_OUTPUT_NAME "kissfft_${KISSFFT_DATATYPE}")
set_target_properties(kissfft PROPERTIES
    OUTPUT_NAME "${KISSFFT_OUTPUT_NAME}"
    DEFINE_SYMBOL KISS_FFT_BUILD
    EXPORT_NAME "kissfft-${KISSFFT_DATATYPE}"
    VERSION ${PROJECT_VERSION})
add_library(kissfft::kissfft ALIAS kissfft)
add_library(kissfft::kissff-${KISSFFT_DATATYPE} ALIAS kissfft)
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    target_link_libraries(kissfft PRIVATE m)
endif()

function(add_kissfft_executable NAME)
    add_executable(${NAME} ${ARGN})
    target_link_libraries(${NAME} PRIVATE kissfft::kissfft)
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        target_link_libraries(${NAME} PRIVATE m)
    endif()
    set_target_properties(${NAME} PROPERTIES
        OUTPUT_NAME "${NAME}_${KISSFFT_DATATYPE}")
endfunction()

option(KISSFFT_TOOLS "Build kissfft tools" ON)
if(KISSFFT_TOOLS)
    add_subdirectory(tools)
endif()

option(KISSFFT_TEST "Build and enable kissfft tests" ON)
if(KISSFFT_TEST)
    enable_testing()
    add_subdirectory(test)
endif()

option(KISSFFT_INSTALL "Enable kissfft install" ON)
if (KISSFFT_INSTALL)
    include(GNUInstallDirs)
    install(TARGETS kissfft EXPORT kissfft
        ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
    install(FILES "kiss_fft.h" "kissfft.hh" DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")

    set(KISSFFT_INSTALL_CMAKE "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}" CACHE FILEPATH "Install destination of kissfft cmake modules")
    mark_as_advanced(KISSFFT_INSTALL_CMAKE)

    install(EXPORT kissfft DESTINATION "${KISSFFT_INSTALL_CMAKE}"
        NAMESPACE "kissfft::"
        FILE "${PROJECT_NAME}-${KISSFFT_DATATYPE}${KISSFFT_EXPORT_SUFFIX}-targets.cmake")
    include(CMakePackageConfigHelpers)
    configure_package_config_file(kissfft-config.cmake.in kissfft-config.cmake
        INSTALL_DESTINATION "${KISSFFT_INSTALL_CMAKE}")
    write_basic_package_version_file(kissfft-config-version.cmake COMPATIBILITY AnyNewerVersion)
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/kissfft-config.cmake" "${CMAKE_CURRENT_BINARY_DIR}/kissfft-config-version.cmake"
        DESTINATION "${KISSFFT_INSTALL_CMAKE}")

    set(PKG_KISSFFT_DEFS)
    foreach(_def ${KISSFFT_COMPILE_DEFINITIONS})
        set(PKG_KISSFFT_DEFS "${PKG_KISSFFT_DEFS} -D${_def}")
    endforeach()
    configure_file(kissfft.pc.in "kissfft-${KISSFFT_DATATYPE}.pc" @ONLY)
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/kissfft-${KISSFFT_DATATYPE}.pc"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig")
endif()
