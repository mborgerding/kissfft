cmake_minimum_required(VERSION 2.8)
project(kissfft)

add_library(kissfft
  kiss_fft.c)
add_library(kissfft::kissfft ALIAS kissfft)

target_include_directories(kissfft PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:.>)

set_target_properties(kissfft PROPERTIES
    DEFINE_SYMBOL KISS_FFT_BUILD)
if(BUILD_SHARED_LIBS)
    target_compile_definitions(kissfft PUBLIC KISS_FFT_SHARED)
    set_target_properties(kissfft PROPERTIES
        C_VISIBILITY_PRESET hidden)
endif()

option(KISSFFT_INSTALL "Enable kissfft install" ON)
if (KISSFFT_INSTALL)
    include(GNUInstallDirs)
    install(TARGETS kissfft EXPORT kissfft
        ARCHIVE DESTINATION "${CMAKE_INSTALL_BINDIR}"
        LIBRARY DESTINATION "${CMAKE_INSTALL_BINDIR}"
        RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
    install(FILES "kiss_fft.h" DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
    install(EXPORT kissfft DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
        NAMESPACE "kissfft::"
        FILE "${PROJECT_NAME}-config.cmake")
endif()
