cmake_minimum_required(VERSION 2.8)
project(kissfft)

add_library(kissfft
  kiss_fft.c)

target_include_directories(kissfft PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:.>)
