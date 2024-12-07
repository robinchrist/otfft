cmake_minimum_required(VERSION 3.15)
project(otfftpp CXX)

option(OTFFT_BUILD_SIMPLE_BENCHMARK "Build simple benchmark" OFF)
option(OTFFT_BUILD_TESTS "Build tests" OFF)

add_library(
    otfftpp
)
add_library(otfftpp::otfftpp ALIAS otfftpp)

target_sources(
    otfftpp
    PRIVATE
        src/otfftpp.cpp
    PUBLIC
        FILE_SET 
            HEADERS
        BASE_DIRS
            include
        FILES
            include/otfftpp/msleep.h
            include/otfftpp/otfft_avxdif16.h
            include/otfftpp/otfft_avxdif4.h
            include/otfftpp/otfft_avxdif8.h
            include/otfftpp/otfft_avxdit16.h
            include/otfftpp/otfft_avxdit4.h
            include/otfftpp/otfft_avxdit8.h
            include/otfftpp/otfft_complex.h
            include/otfftpp/otfft_eightstep.h
            include/otfftpp/otfft_gen_delete.h
            include/otfftpp/otfft_gen_fwd0.h
            include/otfftpp/otfft_gen_fwd.h
            include/otfftpp/otfft_gen_fwdu.h
            include/otfftpp/otfft_gen_inv.h
            include/otfftpp/otfft_gen_invn.h
            include/otfftpp/otfft_gen_invu.h
            include/otfftpp/otfft_gen_new.h
            include/otfftpp/otfft_gen_setup.h
            include/otfftpp/otfft.h
            include/otfftpp/otfft_misc.h
            include/otfftpp/otfft_mixedradix.h
            include/otfftpp/otfft_sixstep.h
            include/otfftpp/otfft_sixstepsq.h
)

# Only build-time dep, soon to be replaced
target_include_directories(otfftpp PRIVATE external/enoki/include)

#We can safely set -fno-math-errno, as darwin does not support it anyways!
target_compile_options(otfftpp PRIVATE "-fno-math-errno")
target_compile_definitions(otfftpp PRIVATE "_USE_MATH_DEFINES")


install(
    TARGETS 
        otfftpp
    FILE_SET
        HEADERS #Weird that has to be added explicitly?
)


if(OTFFT_BUILD_SIMPLE_BENCHMARK)
    add_executable(
        benchmark
        benchmark/benchmark.cpp
    )
    target_link_libraries(benchmark otfftpp::otfftpp)
endif()


if(OTFFT_BUILD_TESTS)
    message(STATUS "Generating tests")
    enable_testing()
    add_subdirectory(test)
endif()
