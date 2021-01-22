add_library(kissfft_help_library OBJECT
    pstats.c
    ../tools/kfc.c
    ../tools/kiss_fftnd.c
    ../tools/kiss_fftr.c
    ../tools/kiss_fftndr.c
)
target_link_libraries(kissfft_help_library PRIVATE kissfft::kissfft)
target_include_directories(kissfft_help_library PRIVATE ../tools)

function(add_kissfft_test_executable NAME)
    add_kissfft_executable(${NAME} ${ARGN} $<TARGET_OBJECTS:kissfft_help_library>)
    target_include_directories(${NAME} PRIVATE ../tools)
    add_test(NAME ${NAME} COMMAND ${NAME})
    set_tests_properties(${NAME} PROPERTIES TIMEOUT 3600)
endfunction()

set(KISSFFT_TEST_NUMFFTS 10000)

add_kissfft_test_executable(bm_kiss benchkiss.c)
# add_test(NAME benchmar COMMAND ${NAME})
# set_tests_properties(${NAME} PROPERTIES TIMEOUT 3600)

include(FindPkgConfig)
if(KISSFFT_FLOAT)
    set(fftw3_pkg fftw3f)
else()
    set(fftw3_pkg fftw3)
endif()
pkg_check_modules(fftw3 REQUIRED IMPORTED_TARGET ${fftw3_pkg})
add_kissfft_test_executable(bm_fftw benchfftw.c)
target_link_libraries(bm_fftw PRIVATE PkgConfig::fftw3)

add_kissfft_test_executable(st twotonetest.c)

add_kissfft_test_executable(tkfc twotonetest.c)
target_compile_definitions(tkfc PRIVATE KFC_TEST)

add_kissfft_test_executable(ffr twotonetest.c)
add_kissfft_test_executable(tr test_real.c)

add_kissfft_test_executable(testcpp testcpp.cc)

find_package(PythonInterp REQUIRED)
add_test(NAME testkiss.py COMMAND "${PYTHON_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/testkiss.py")
set_tests_properties(testkiss.py PROPERTIES
    TIMEOUT 3600
    ENVIRONMENT "DATATYPE=${KISSFFT_DATATYPE}"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
