#include "otfftpp/otfft.h"
#include <iostream>
#include <chrono>

int main() {

    const size_t maxExp = 24;

    for(size_t exp = 0; exp <= maxExp; ++exp) {
        size_t SIZE = 1 << exp;
        OTFFT_Complex::simd_array<OTFFT::complex_t> workspace(SIZE);

        auto FFT = std::make_unique<OTFFT::FFT>(SIZE);

        FFT->fwd(&workspace[0]);
        FFT->fwd(&workspace[0]);
        FFT->fwd(&workspace[0]);
        FFT->fwd(&workspace[0]);

        auto start = std::chrono::high_resolution_clock::now();
        FFT->fwd(&workspace[0]);
        auto end = std::chrono::high_resolution_clock::now();
 
        std::cout << "Size 2^" << exp << " = " << SIZE << "-> " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << "µS" << std::endl;
    }

    return 0;
}