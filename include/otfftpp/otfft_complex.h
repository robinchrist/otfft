/******************************************************************************
*  OTFFT Complex & Memory Allocator Version 11.5e
*
*  Copyright (c) 2019 OK Ojisan(Takuya OKAHISA)
*  Released under the MIT license
*  http://opensource.org/licenses/mit-license.php
******************************************************************************/

#ifndef otfft_complex_h
#define otfft_complex_h

#ifdef _MSC_VER
#if _MSC_VER >= 1900
#define VC_CONSTEXPR 1
#else
#error "This compiler is not supported."
#endif
#endif // _MSC_VER

#if __cplusplus >= 201103L || defined(VC_CONSTEXPR)
#define NOEXCEPT noexcept
#else
#define NOEXCEPT
#endif

#if __GNUC__ >= 3
//#define force_inline  __attribute__((const,always_inline))
//#define force_inline2 __attribute__((pure,always_inline))
//#define force_inline3 __attribute__((always_inline))
//#define force_inline  __attribute__((const))
//#define force_inline2 __attribute__((pure))
//#define force_inline3
#define force_inline
#define force_inline2
#define force_inline3
#else
#define force_inline
#define force_inline2
#define force_inline3
#endif

//=============================================================================
// User Defined Complex Number Class
//=============================================================================

#include <complex>

namespace OTFFT_Complex {

struct complex_t
{
    double Re, Im;

    complex_t() NOEXCEPT : Re(0), Im(0) {}
    complex_t(const double& x) NOEXCEPT : Re(x), Im(0) {}
    complex_t(const double& x, const double& y) NOEXCEPT : Re(x), Im(y) {}
    complex_t(const complex_t& z) NOEXCEPT : Re(z.Re), Im(z.Im) {}
    complex_t(const std::complex<double>& z) NOEXCEPT : Re(z.real()), Im(z.imag()) {}
    operator std::complex<double>() const { return std::complex<double>(Re, Im); }

    complex_t& operator=(const complex_t& z) NOEXCEPT
    {
        Re = z.Re;
        Im = z.Im;
        return *this;
    }

    complex_t& operator+=(const complex_t& z) NOEXCEPT
    {
        Re += z.Re;
        Im += z.Im;
        return *this;
    }

    complex_t& operator-=(const complex_t& z) NOEXCEPT
    {
        Re -= z.Re;
        Im -= z.Im;
        return *this;
    }

    complex_t& operator*=(const double& x) NOEXCEPT
    {
        Re *= x;
        Im *= x;
        return *this;
    }

    complex_t& operator/=(const double& x) NOEXCEPT
    {
        Re /= x;
        Im /= x;
        return *this;
    }

    complex_t& operator*=(const complex_t& z) NOEXCEPT
    {
        const double tmp = Re*z.Re - Im*z.Im;
        Im = Re*z.Im + Im*z.Re;
        Re = tmp;
        return *this;
    }
};

typedef double* __restrict const double_vector;
typedef const double* __restrict const const_double_vector;
typedef complex_t* __restrict const complex_vector;
typedef const complex_t* __restrict const const_complex_vector;

static inline double Re(const complex_t& z) NOEXCEPT force_inline;
static inline double Re(const complex_t& z) NOEXCEPT { return z.Re; }
static inline double Im(const complex_t& z) NOEXCEPT force_inline;
static inline double Im(const complex_t& z) NOEXCEPT { return z.Im; }

static inline double norm(const complex_t& z) NOEXCEPT force_inline;
static inline double norm(const complex_t& z) NOEXCEPT
{
    return z.Re*z.Re + z.Im*z.Im;
}
static inline complex_t conj(const complex_t& z) NOEXCEPT force_inline;
static inline complex_t conj(const complex_t& z) NOEXCEPT
{
    return complex_t(z.Re, -z.Im);
}
static inline complex_t jx(const complex_t& z) NOEXCEPT force_inline;
static inline complex_t jx(const complex_t& z) NOEXCEPT
{
    return complex_t(-z.Im, z.Re);
}
static inline complex_t neg(const complex_t& z) NOEXCEPT force_inline;
static inline complex_t neg(const complex_t& z) NOEXCEPT
{
    return complex_t(-z.Re, -z.Im);
}
static inline complex_t mjx(const complex_t& z) NOEXCEPT force_inline;
static inline complex_t mjx(const complex_t& z) NOEXCEPT
{
    return complex_t(z.Im, -z.Re);
}
static inline complex_t operator+(const complex_t& a, const complex_t& b) NOEXCEPT force_inline;
static inline complex_t operator+(const complex_t& a, const complex_t& b) NOEXCEPT
{
    return complex_t(a.Re + b.Re, a.Im + b.Im);
}
static inline complex_t operator-(const complex_t& a, const complex_t& b) NOEXCEPT force_inline;
static inline complex_t operator-(const complex_t& a, const complex_t& b) NOEXCEPT
{
    return complex_t(a.Re - b.Re, a.Im - b.Im);
}
static inline complex_t operator*(const double& a, const complex_t& b) NOEXCEPT force_inline;
static inline complex_t operator*(const double& a, const complex_t& b) NOEXCEPT
{
    return complex_t(a*b.Re, a*b.Im);
}
static inline complex_t operator*(const complex_t& a, const complex_t& b) NOEXCEPT force_inline;
static inline complex_t operator*(const complex_t& a, const complex_t& b) NOEXCEPT
{
    return complex_t(a.Re*b.Re - a.Im*b.Im, a.Re*b.Im + a.Im*b.Re);
}
static inline complex_t operator/(const complex_t& a, const double& b) NOEXCEPT force_inline;
static inline complex_t operator/(const complex_t& a, const double& b) NOEXCEPT
{
    return complex_t(a.Re/b, a.Im/b);
}
static inline complex_t operator/(const complex_t& a, const complex_t& b) NOEXCEPT force_inline;
static inline complex_t operator/(const complex_t& a, const complex_t& b) NOEXCEPT
{
    const double b2 = b.Re*b.Re + b.Im*b.Im;
    return (a * conj(b)) / b2;
}

static inline complex_t expj(const double& theta) NOEXCEPT force_inline;
static inline complex_t expj(const double& theta) NOEXCEPT
{
    //return complex_t(cos(theta), sin(theta));
    return complex_t(std::polar(1.0, theta));
}

} // namespace OTFFT_Complex

//=============================================================================
// Aligned Memory Allocator
//=============================================================================

#include <new>
#include <memory>
#include <cstring>

#ifdef __MINGW32__
#include <malloc.h>
#endif

namespace OTFFT_Complex {

inline void* generic_aligned_alloc(std::size_t size, std::size_t alignment) noexcept
{
    constexpr size_t N = alignof(void*);

    if (alignment < N) {
        alignment = N;
    }

    std::size_t n = size + alignment - N;
    void* p = std::malloc(sizeof(void*) + n);
    if (p) {
        void* p2 = static_cast<char*>(p) + sizeof(p);
        auto * const ap = static_cast<char*>(std::align(alignment, size, p2, n));
        std::memcpy(ap - sizeof(p), &p, sizeof(p));

        p = ap;
    }
    return p;
}

inline void generic_aligned_free(void* ptr) noexcept
{
    if(ptr) {
        void* rp;
        std::memcpy(&rp, static_cast<char*>(ptr) - sizeof(void*), sizeof(void*));
        std::free(rp);
    }
}

#ifdef __AVX512F__
static inline void* simd_malloc(const size_t n) { return generic_aligned_alloc(n, 64); }
#else
    #ifdef __AVX__
    static inline void* simd_malloc(const size_t n) { return generic_aligned_alloc(n, 32); }
    #else
    static inline void* simd_malloc(const size_t n) { return generic_aligned_alloc(n, 16); }
    #endif
#endif

static inline void simd_free(void* p) { generic_aligned_free(p); }

template <class T> 
class simd_array {

public:
    simd_array() NOEXCEPT : alignedArray(nullptr) {}

    simd_array(int n)
    {
        setup(n);
    }

    void setup(int n)
    {
        T* newP = (T*) simd_malloc(n*sizeof(T));
        if (newP == nullptr) throw std::bad_alloc();
        //Destroys the old array and replaces it with the new one
        alignedArray.reset(newP);
    }

    void destroy() { 
        alignedArray.reset();
    }

    T& operator[](int i) NOEXCEPT { 
        return alignedArray.get()[i];
    }
    const T& operator[](int i) const NOEXCEPT { 
        return alignedArray.get()[i];
    }
    T* operator&() const NOEXCEPT = delete;

    T* get() const noexcept {
        return alignedArray.get();
    }

private:

    struct SIMDDeleter {
        void operator()(T* ptr) { 
            if(ptr) simd_free(ptr);
        }
    };

    std::unique_ptr<T, SIMDDeleter> alignedArray;
};

} // namespace OTFFT_Complex

//=============================================================================

#endif // otfft_complex_h
