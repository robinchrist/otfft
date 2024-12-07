// Copyright (c) Dewetron 2017
#include "otfftpp/otfft.h"

#undef VALGRIND_MEM_LEAK_DETECTION

#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

namespace
{
    /**
     * Test that spectrum of a delta signal is constant.
     */
    template <std::size_t SIZE>
    void deltaSpectrumTest()
    {
        // delta signal: 1, 0, 0, 0...
        OTFFT_Complex::simd_array<double> test_array(SIZE);
        test_array[0] = 1.0;
        for(size_t index = 1; index < SIZE; ++index) {
            test_array[index] = 0.0;
        }

        OTFFT_Complex::simd_array<OTFFT::complex_t> spectrum(SIZE);
        {
            auto fft = std::make_unique<OTFFT::RFFT>(static_cast<int>(SIZE));
            OTFFT::double_vector fft_in{&test_array[0]};
            OTFFT::complex_vector fft_out{&spectrum[0]};
            fft->fwd0(fft_in, fft_out);
        }

        double absSpectrum[SIZE] = {0};
        std::transform(&spectrum[0], &spectrum[0]+SIZE, absSpectrum, [] (OTFFT::complex_t x) {
            return std::sqrt(OTFFT::norm(x));
        });

        double expected[SIZE];
        std::fill_n(expected, SIZE, 1.0);
        for (std::size_t idx{0}; idx < SIZE / 2; ++idx)
        {
            if (std::fabs(expected[idx]) < 1e-10)
            {
                BOOST_CHECK_SMALL(absSpectrum[idx], 1e-8);
            }
            else
            {
                BOOST_CHECK_CLOSE(expected[idx], absSpectrum[idx], .1);
            }
        }
    }

    /**
     * Test that spectrum of a constant signal is a scaled delta.
     */
    template <std::size_t SIZE>
    void constSpectrumTest()
    {
        OTFFT_Complex::simd_array<double> testArray(SIZE);
        for(size_t index = 0; index < SIZE; ++index) {
            testArray[index] = 1.0;
        }

        OTFFT_Complex::simd_array<OTFFT::complex_t> spectrum(SIZE);
        {
            auto fft = std::make_unique<OTFFT::RFFT>(static_cast<int>(SIZE));
            OTFFT::double_vector fft_in{&testArray[0]};
            OTFFT::complex_vector fft_out{&spectrum[0]};
            fft->fwd0(fft_in, fft_out);
        }

        double absSpectrum[SIZE] = {0};
        std::transform(&spectrum[0], &spectrum[0]+SIZE, absSpectrum, [] (OTFFT::complex_t x) {
            return std::sqrt(OTFFT::norm(x));
        });

        // expecting a delta scaled by SIZE
        double expected[SIZE] = {0};
        expected[0] = SIZE * 1.0;
        for (std::size_t idx{0}; idx < SIZE / 2; ++idx)
        {
            if (std::fabs(expected[idx]) < 1e-10)
            {
                BOOST_CHECK_SMALL(absSpectrum[idx], 1e-8);
            }
            else
            {
                BOOST_CHECK_CLOSE(expected[idx], absSpectrum[idx], .1);
            }
        }
    }

    /**
     * Test that inverse of a delta spectrum is a constant signal.
     */
    template <std::size_t SIZE>
    void deltaInverseTest()
    {
        OTFFT_Complex::simd_array<OTFFT::complex_t> spectrum(SIZE);
        spectrum[0] = SIZE * 1.0;
        for(size_t index = 1; index < SIZE; ++index) {
            spectrum[index] = 0.0;
        }

        OTFFT_Complex::simd_array<double> output(SIZE);
        {
            auto fft = std::make_unique<OTFFT::RFFT>(static_cast<int>(SIZE));
            OTFFT::complex_vector fft_in{&spectrum[0]};
            OTFFT::double_vector fft_out{&output[0]};
            fft->invn(fft_in, fft_out);
        }

        double expected[SIZE];
        std::fill_n(expected, SIZE, 1.0);
        for (std::size_t idx{0}; idx < SIZE / 2; ++idx)
        {
            if (std::fabs(expected[idx]) < 1e-10)
            {
                BOOST_CHECK_SMALL(output[idx], 1e-8);
            }
            else
            {
                BOOST_CHECK_CLOSE(expected[idx], output[idx], .1);
            }
        }
    }

    /**
     * Test that inverse of a constant spectrum is a delta signal.
     */
    template <std::size_t SIZE>
    void constInverseTest()
    {
        OTFFT_Complex::simd_array<OTFFT::complex_t> spectrum(SIZE);
        for(size_t index = 0; index < SIZE; ++index) {
            spectrum[index] = 1.0;
        }

        OTFFT_Complex::simd_array<double> output(SIZE);
        {
            auto fft = std::make_unique<OTFFT::RFFT>(static_cast<int>(SIZE));
            OTFFT::complex_vector fft_in{&spectrum[0]};
            OTFFT::double_vector fft_out{&output[0]};
            fft->invn(fft_in, fft_out);
        }

        double expected[SIZE] = {0};
        expected[0] = 1.0;
        for (std::size_t idx{0}; idx < SIZE / 2; ++idx)
        {
            if (std::fabs(expected[idx]) < 1e-10)
            {
                BOOST_CHECK_SMALL(output[idx], 1e-8);
            }
            else
            {
                BOOST_CHECK_CLOSE(expected[idx], output[idx], .1);
            }
        }
    }

    /**
     * Test that IFFT(FFT(x)) == x.
     */
    template <std::size_t SIZE>
    void identityTest()
    {
        OTFFT_Complex::simd_array<double> testArray(SIZE);
        for(size_t index = 0; index < SIZE; ++index) {
            testArray[index] = 1.0;
        }


        OTFFT_Complex::simd_array<OTFFT::complex_t> spectrum(SIZE);
        {
            auto fft = std::make_unique<OTFFT::RFFT>(static_cast<int>(SIZE));
            OTFFT::double_vector fft_in{&testArray[0]};
            OTFFT::complex_vector fft_out{&spectrum[0]};
            fft->fwd(fft_in, fft_out);
        }

        OTFFT_Complex::simd_array<double> output(SIZE);
        {
            auto fft = std::make_unique<OTFFT::RFFT>(static_cast<int>(SIZE));
            OTFFT::complex_vector fft_in{&spectrum[0]};
            OTFFT::double_vector fft_out{&output[0]};
            fft->inv(fft_in, fft_out);
        }

        for (std::size_t idx{0}; idx < SIZE / 2; ++idx)
        {
            BOOST_CHECK_CLOSE(testArray[idx], output[idx], .1);
        }
    }

    void fftUtilizeReal(std::size_t fft_size)
    {
        OTFFT_Complex::simd_array<double> workspace_real(fft_size);
        for(size_t index = 0; index < fft_size; ++index) {
            workspace_real[index] = 1.0;
        }
        OTFFT_Complex::simd_array<OTFFT::complex_t> workspace_complex(fft_size);
        {
            auto fft = std::make_unique<OTFFT::RFFT>(static_cast<int>(fft_size));
            OTFFT::double_vector fft_in{&workspace_real[0]};
            OTFFT::complex_vector fft_out{&workspace_complex[0]};
            fft->fwd(fft_in, fft_out);
        }

        {
            auto fft = std::make_unique<OTFFT::RFFT>(static_cast<int>(fft_size));
            OTFFT::complex_vector fft_in{&workspace_complex[0]};
            OTFFT::double_vector fft_out{&workspace_real[0]};
            fft->inv(fft_in, fft_out);
        }
    }

    void fftUtilizeComplex(std::size_t fft_size)
    {
        OTFFT_Complex::simd_array<OTFFT::complex_t> workspace(fft_size);
        for(size_t index = 0; index < fft_size; ++index) {
            workspace[index] = OTFFT::complex_t(1.0, 1.0);
        }

        {
            auto fft = std::make_unique<OTFFT::FFT>(static_cast<int>(fft_size));
            OTFFT::complex_vector workspace_ptr{&workspace[0]};
            fft->fwd(workspace_ptr);
        }

        {
            auto fft = std::make_unique<OTFFT::FFT>(static_cast<int>(fft_size));
            OTFFT::complex_vector workspace_ptr{&workspace[0]};
            fft->inv(workspace_ptr);
        }
    }

    void fftUtilizeDCT(std::size_t fft_size)
    {
        OTFFT_Complex::simd_array<double> workspace(fft_size);
        for(size_t index = 0; index < fft_size; ++index) {
            workspace[index] = 1.0;
        }

        {
            auto fft = std::make_unique<OTFFT::DCT>(static_cast<int>(fft_size));
            OTFFT::double_vector workspace_ptr{&workspace[0]};
            fft->fwd(workspace_ptr);
        }

        {
            auto fft = std::make_unique<OTFFT::DCT>(static_cast<int>(fft_size));
            OTFFT::double_vector workspace_ptr{&workspace[0]};
            fft->inv(workspace_ptr);
        }
    }

    void fftUtilizeBluestein(std::size_t fft_size)
    {
        OTFFT_Complex::simd_array<OTFFT::complex_t> workspace(fft_size);
        for(size_t index = 0; index < fft_size; ++index) {
            workspace[index] = OTFFT::complex_t(1.0, 1.0);
        }

        {
            auto fft = std::make_unique<OTFFT::Bluestein>(static_cast<int>(fft_size));
            OTFFT::complex_vector workspace_ptr{&workspace[0]};
            fft->fwd(workspace_ptr);
        }

        {
            auto fft = std::make_unique<OTFFT::Bluestein>(static_cast<int>(fft_size));
            OTFFT::complex_vector workspace_ptr{&workspace[0]};
            fft->inv(workspace_ptr);
        }
    }
}

BOOST_AUTO_TEST_SUITE(otfft_transform_test)

BOOST_AUTO_TEST_CASE(TestDeltaSpectrum)
{
    deltaSpectrumTest<8>();
    deltaSpectrumTest<16>();
    deltaSpectrumTest<256>();
    deltaSpectrumTest<512>();
    deltaSpectrumTest<1024>();
    deltaSpectrumTest<8192>();
    deltaSpectrumTest<16384>();
}

BOOST_AUTO_TEST_CASE(TestConstSpectrum)
{
    constSpectrumTest<8>();
    constSpectrumTest<16>();
    constSpectrumTest<256>();
    constSpectrumTest<512>();
    constSpectrumTest<1024>();
    constSpectrumTest<8192>();
    constSpectrumTest<16384>();
}

BOOST_AUTO_TEST_CASE(TestDeltaInverse)
{
    deltaInverseTest<8>();
    deltaInverseTest<16>();
    deltaInverseTest<256>();
    deltaInverseTest<512>();
    deltaInverseTest<1024>();
    deltaInverseTest<8192>();
    deltaInverseTest<16384>();
}

BOOST_AUTO_TEST_CASE(TestConstInverse)
{
    constInverseTest<8>();
    constInverseTest<16>();
    constInverseTest<256>();
    constInverseTest<512>();
    constInverseTest<1024>();
    constInverseTest<8192>();
    constInverseTest<16384>();
}

BOOST_AUTO_TEST_CASE(TestIdentity)
{
    identityTest<8>();
    identityTest<16>();
    identityTest<256>();
    identityTest<512>();
    identityTest<1024>();
    identityTest<8192>();
    identityTest<16384>();
}

#ifdef VALGRIND_MEM_LEAK_DETECTION
BOOST_AUTO_TEST_CASE(TestValgrind)
{
    const std::size_t N = static_cast<std::size_t>(std::pow(2, 24));

    for (std::size_t n = 8; n < N;)
    {
        fftUtilizeComplex(n);
        fftUtilizeBluestein(n);

        if ((n & 1) == 0)
        {
            fftUtilizeReal(n);
            fftUtilizeDCT(n);
        }

        const std::size_t log_n = static_cast<std::size_t>(std::log2(n));
        n += static_cast<std::size_t>(std::ceil(std::pow(log_n, std::sqrt(log_n))));
    }

    BOOST_CHECK(true);
}
#endif

BOOST_AUTO_TEST_SUITE_END()
