/******************************************************************************
*  OTFFT Miscellaneous Routines Version 11.5e
*
*  Copyright (c) 2019 OK Ojisan(Takuya OKAHISA)
*  Released under the MIT license
*  http://opensource.org/licenses/mit-license.php
******************************************************************************/

#ifndef otfft_misc_h
#define otfft_misc_h

//=============================================================================
// Customization Options
//=============================================================================

//#define DO_SINGLE_THREAD 1
//#define USE_UNALIGNED_MEMORY 1

//=============================================================================

#include <cmath>
#include <enoki/array.h>
#include <enoki/complex.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

#ifndef M_SQRT2
#define M_SQRT2 1.41421356237309504876378807303183294
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2 0.707106781186547524400844362104849039
#endif

#include "otfft_complex.h"

namespace OTFFT_MISC {

using namespace OTFFT_Complex;

static const double H1X =  0.923879532511286762010323247995557949;
static const double H1Y = -0.382683432365089757574419179753100195;

enum scaling_mode { scale_1 = 0, scale_unitary, scale_length };

static inline complex_t v8x(const complex_t& z) NOEXCEPT force_inline;
static inline complex_t v8x(const complex_t& z) NOEXCEPT
{
    return complex_t(M_SQRT1_2*(z.Re-z.Im), M_SQRT1_2*(z.Re+z.Im));
}
static inline complex_t w8x(const complex_t& z) NOEXCEPT force_inline;
static inline complex_t w8x(const complex_t& z) NOEXCEPT
{
    return complex_t(M_SQRT1_2*(z.Re+z.Im), M_SQRT1_2*(z.Im-z.Re));
}

} // namespace OTFFT_MISC

//=============================================================================
// constexpr sqrt
//=============================================================================

#if __cplusplus >= 201103L || defined(VC_CONSTEXPR)
namespace OTFFT_MISC {

constexpr double sqrt_aux(double a, double x, double y)
{
    return x == y ? x : sqrt_aux(a, (x + a/x)/2, x);
}

constexpr double mysqrt(double x) { return sqrt_aux(x, x/2, x); }

} // namespace OTFFT_MISC
#endif

//=============================================================================
// FFT Weight Initialize Routine
//=============================================================================

namespace OTFFT_MISC {

    static void init_Wt(const int r, const int N, complex_vector W) noexcept
    {
        if (N < r) return;
        const int Nr = N/r;
        const double theta = -2*M_PI/N;

        for (int p = 0; p < Nr; p++) {
            for (int k = 1; k < r; k++) {
                W[p + (k-1)*Nr] = W[N + r*p + k] = expj(theta * k*p);
            }
        }
    }

    template <int r, int N, int k>
    const_complex_vector twid(const_complex_vector W, int p)
    {
        constexpr int Nr = N/r;
        constexpr int d = (k-1)*Nr;
        return W + p + d;
    }

    template <int r, int N, int k>
    const_complex_vector twidT(const_complex_vector W, int p)
    {
        constexpr int d = N + k;
        return W + r*p + d;
    }

} // namespace OTFFT_MISC

//=============================================================================
// SSE2/SSE3
//=============================================================================

namespace OTFFT_MISC {

    using Vec2d = enoki::Array<double, 2>;

    static inline Vec2d cmplx(const double& x, const double& y) NOEXCEPT force_inline;
    static inline Vec2d cmplx(const double& x, const double& y) NOEXCEPT
    {
        return Vec2d {x, y};
    }

    static inline Vec2d getpz(const complex_t& z) NOEXCEPT force_inline;
    static inline Vec2d getpz(const complex_t& z) NOEXCEPT
    {
    #ifdef USE_UNALIGNED_MEMORY
        return enoki::load_unaligned<Vec2d>(&z.Re);
    #else
        return enoki::load<Vec2d>(&z.Re);
    #endif
    }
    static inline Vec2d getpz(const_double_vector x) NOEXCEPT force_inline;
    static inline Vec2d getpz(const_double_vector x) NOEXCEPT
    {
    #ifdef USE_UNALIGNED_MEMORY
        return enoki::load_unaligned<Vec2d>(x);
    #else
        return enoki::load<Vec2d>(x);
    #endif
    }

    static inline void setpz(complex_t& z, const Vec2d x) NOEXCEPT force_inline3;
    static inline void setpz(complex_t& z, const Vec2d x) NOEXCEPT
    {
    #ifdef USE_UNALIGNED_MEMORY
        enoki::store_unaligned<Vec2d>(&z.Re, x);
    #else
        enoki::store<Vec2d>(&z.Re, x);
    #endif
    }
    static inline void setpz(double_vector x, const Vec2d z) NOEXCEPT force_inline3;
    static inline void setpz(double_vector x, const Vec2d z) NOEXCEPT
    {
    #ifdef USE_UNALIGNED_MEMORY
        enoki::store_unaligned<Vec2d>(x, z);
    #else
        enoki::store<Vec2d>(x, z);
    #endif
    }

    static inline Vec2d cnjpz(const Vec2d xy) NOEXCEPT force_inline;
    //Conjugate packed complex?
    static inline Vec2d cnjpz(const Vec2d xy) NOEXCEPT
    {
        return Vec2d { xy[0], -xy[1] };
    }
    static inline Vec2d jxpz(const Vec2d xy) NOEXCEPT force_inline;
    //Multiply packed complex with j / Rotate by 90° (pi/2) counterclockwise?
    static inline Vec2d jxpz(const Vec2d xy) NOEXCEPT
    {
        return Vec2d { -xy[1], xy[0] };
    }

    static inline Vec2d addpz(const Vec2d a, const Vec2d b) NOEXCEPT force_inline;
    static inline Vec2d addpz(const Vec2d a, const Vec2d b) NOEXCEPT
    {
        return a + b;
    }
    static inline Vec2d subpz(const Vec2d a, const Vec2d b) NOEXCEPT force_inline;
    static inline Vec2d subpz(const Vec2d a, const Vec2d b) NOEXCEPT
    {
        return a - b;
    }
    static inline Vec2d mulpd(const Vec2d a, const Vec2d b) NOEXCEPT force_inline;
    static inline Vec2d mulpd(const Vec2d a, const Vec2d b) NOEXCEPT
    {
        return a * b;
    }

    #if __cplusplus >= 201103L || defined(VC_CONSTEXPR)
    template <int N, int mode> static inline Vec2d scalepz(const Vec2d z) NOEXCEPT force_inline;
    template <int N, int mode> static inline Vec2d scalepz(const Vec2d z) NOEXCEPT
    {
        constexpr double scale =
            mode == scale_1       ? 1.0           :
            mode == scale_unitary ? 1.0/mysqrt(N) :
            mode == scale_length  ? 1.0/N         : 0.0;
        const Vec2d sv = { scale, scale };
        return mode == scale_1 ? z : mulpd(sv, z);
    }
    #endif

} // namespace OTFFT_MISC

namespace OTFFT_MISC {

    static inline Vec2d haddpz(const Vec2d ab, const Vec2d xy) NOEXCEPT force_inline;
    //arg ab = (a, b), xy = (x, y) -> returns (a + b, x + z)
    static inline Vec2d haddpz(const Vec2d ab, const Vec2d xy) NOEXCEPT
    {
        //Automatically optimised to vhaddpd by compiler
        return Vec2d { ab[0]+ab[1], xy[0]+xy[1] }; // (a + b, x + y)
    }

    static inline Vec2d mulpz(const Vec2d ab, const Vec2d xy) NOEXCEPT force_inline;
    //Complex multiplication, ab = (a, b) = a + i*b, xy = (x, y) = x + i*y -> returns (a + i*b) * (x + i*y) = (a*x-b*y, a*y+b*x)
    static inline Vec2d mulpz(const Vec2d ab, const Vec2d xy) NOEXCEPT
    {
        enoki::Complex<float> mulRes = enoki::Complex<float>(ab[0], ab[1]) * enoki::Complex<float>(xy[0], xy[1]);
        return Vec2d { enoki::real(mulRes), enoki::imag(mulRes)};
    }

    static inline Vec2d v8xpz(const Vec2d xy) NOEXCEPT force_inline;
    //Computes (1/sqrt(2)) * (x-y, x+y)?
    static inline Vec2d v8xpz(const Vec2d xy) NOEXCEPT
    {
        const Vec2d rr { M_SQRT1_2, M_SQRT1_2 };
        return rr * (xy + jxpz(xy));
    }

} // namespace OTFFT_MISC

//-----------------------------------------------------------------------------

namespace OTFFT_MISC {

    static inline Vec2d w8xpz(const Vec2d xy) NOEXCEPT force_inline;
    //Calculates 1/sqrt(2) * (x + y, x - y)?
    static inline Vec2d w8xpz(const Vec2d xy) NOEXCEPT
    {
        return { M_SQRT1_2*(xy[0] + xy[1]), M_SQRT1_2*(xy[1] - xy[0]) };
    }

    static inline Vec2d h1xpz(const Vec2d xy) NOEXCEPT force_inline;
    static inline Vec2d h1xpz(const Vec2d xy) NOEXCEPT
    {
        static const Vec2d h1 = { H1X, H1Y };
        return mulpz(h1, xy);
    }

    static inline Vec2d h3xpz(const Vec2d xy) NOEXCEPT force_inline;
    static inline Vec2d h3xpz(const Vec2d xy) NOEXCEPT
    {
        static const Vec2d h3 = { -H1Y, -H1X };
        return mulpz(h3, xy);
    }

    static inline Vec2d hfxpz(const Vec2d xy) NOEXCEPT force_inline;
    static inline Vec2d hfxpz(const Vec2d xy) NOEXCEPT
    {
        static const Vec2d hf = { H1X, -H1Y };
        return mulpz(hf, xy);
    }

    static inline Vec2d hdxpz(const Vec2d xy) NOEXCEPT force_inline;
    static inline Vec2d hdxpz(const Vec2d xy) NOEXCEPT
    {
        static const Vec2d hd = { -H1Y, H1X };
        return mulpz(hd, xy);
    }

} // namespace OTFFT_MISC


namespace OTFFT_MISC {

    //struct Vec4d { Vec2d lo, hi; };
    using Vec4d = enoki::Array<double, 4>;

    static inline Vec4d getpz2(const_complex_vector z) NOEXCEPT force_inline;
    static inline Vec4d getpz2(const_complex_vector z) NOEXCEPT
    {
        #ifdef USE_UNALIGNED_MEMORY
            return enoki::load_unaligned<Vec4d>(&z->Re);
        #else
            return enoki::load<Vec4d>(&z->Re);
        #endif
    }

    static inline void setpz2(complex_vector z, const Vec4d& y) NOEXCEPT force_inline3;
    static inline void setpz2(complex_vector z, const Vec4d& y) NOEXCEPT
    {
        #ifdef USE_UNALIGNED_MEMORY
            enoki::store_unaligned<Vec4d>(&z->Re, y);
        #else
            enoki::store<Vec4d>(&z->Re, y);
        #endif
    }

    static inline Vec4d cnjpz2(const Vec4d& xy) NOEXCEPT force_inline;
    static inline Vec4d cnjpz2(const Vec4d& xy) NOEXCEPT
    {
        return enoki::concat( cnjpz(enoki::low(xy)), cnjpz(enoki::high(xy)) );
    }
    static inline Vec4d jxpz2(const Vec4d& xy) NOEXCEPT force_inline;
    static inline Vec4d jxpz2(const Vec4d& xy) NOEXCEPT
    {
        return enoki::concat( jxpz(enoki::low(xy)), jxpz(enoki::high(xy)) );
    }

    static inline Vec4d addpz2(const Vec4d& a, const Vec4d& b) NOEXCEPT force_inline;
    static inline Vec4d addpz2(const Vec4d& a, const Vec4d& b) NOEXCEPT
    {
        return a + b;
    }
    static inline Vec4d subpz2(const Vec4d& a, const Vec4d& b) NOEXCEPT force_inline;
    static inline Vec4d subpz2(const Vec4d& a, const Vec4d& b) NOEXCEPT
    {
        return a - b;
    }
    static inline Vec4d mulpd2(const Vec4d& a, const Vec4d& b) NOEXCEPT force_inline;
    static inline Vec4d mulpd2(const Vec4d& a, const Vec4d& b) NOEXCEPT
    {
        return a * b;
    }

    static inline Vec4d mulpz2(const Vec4d& a, const Vec4d& b) NOEXCEPT force_inline;
    static inline Vec4d mulpz2(const Vec4d& a, const Vec4d& b) NOEXCEPT
    {
        return enoki::concat(mulpz(enoki::low(a), enoki::low(b)), mulpz(enoki::high(a), enoki::high(b)));
    }

    #if __cplusplus >= 201103L || defined(VC_CONSTEXPR)
    template <int N, int mode> static inline Vec4d scalepz2(const Vec4d z) NOEXCEPT force_inline;
    template <int N, int mode> static inline Vec4d scalepz2(const Vec4d z) NOEXCEPT
    {
        constexpr double scale =
            mode == scale_1       ? 1.0           :
            mode == scale_unitary ? 1.0/mysqrt(N) :
            mode == scale_length  ? 1.0/N         : 0.0;
        const Vec2d sv  = { scale, scale };
        const Vec4d sv2 = enoki::concat(sv, sv);
        return mode == scale_1 ? z : mulpd2(sv2, z);
    }
    #endif

    static inline Vec4d v8xpz2(const Vec4d& xy) NOEXCEPT force_inline;
    static inline Vec4d v8xpz2(const Vec4d& xy) NOEXCEPT
    {
        return enoki::concat( v8xpz(enoki::low(xy)), v8xpz(enoki::high(xy)) );
    }

    static inline Vec4d w8xpz2(const Vec4d& xy) NOEXCEPT force_inline;
    static inline Vec4d w8xpz2(const Vec4d& xy) NOEXCEPT
    {
        return enoki::concat( w8xpz(enoki::low(xy)), w8xpz(enoki::high(xy)) );
    }

    static inline Vec4d h1xpz2(const Vec4d& xy) NOEXCEPT force_inline;
    static inline Vec4d h1xpz2(const Vec4d& xy) NOEXCEPT
    {
        return enoki::concat( h1xpz(enoki::low(xy)), h1xpz(enoki::high(xy)) );
    }

    static inline Vec4d h3xpz2(const Vec4d& xy) NOEXCEPT force_inline;
    static inline Vec4d h3xpz2(const Vec4d& xy) NOEXCEPT
    {
        return enoki::concat( h3xpz(enoki::low(xy)), h3xpz(enoki::high(xy)) );
    }

    static inline Vec4d hfxpz2(const Vec4d& xy) NOEXCEPT force_inline;
    static inline Vec4d hfxpz2(const Vec4d& xy) NOEXCEPT
    {
        return enoki::concat( hfxpz(enoki::low(xy)), hfxpz(enoki::high(xy)) );
    }

    static inline Vec4d hdxpz2(const Vec4d& xy) NOEXCEPT force_inline;
    static inline Vec4d hdxpz2(const Vec4d& xy) NOEXCEPT
    {
        return enoki::concat( hdxpz(enoki::low(xy)), hdxpz(enoki::high(xy)) );
    }

    static inline Vec4d duppz2(const Vec2d x) NOEXCEPT force_inline;
    static inline Vec4d duppz2(const Vec2d x) NOEXCEPT
    {
        return enoki::concat(x, x);
    }

    static inline Vec4d duppz3(const complex_t& z) NOEXCEPT force_inline;
    static inline Vec4d duppz3(const complex_t& z) NOEXCEPT
    {
        const Vec2d x = getpz(z);
        return enoki::concat(x, x);
    }

    static inline Vec4d cat(const Vec2d& a, const Vec2d& b) NOEXCEPT force_inline;
    static inline Vec4d cat(const Vec2d& a, const Vec2d& b) NOEXCEPT
    {
        return enoki::concat( a, b );
    }

    static inline Vec4d catlo(const Vec4d& ax, const Vec4d& by) NOEXCEPT force_inline;
    static inline Vec4d catlo(const Vec4d& ax, const Vec4d& by) NOEXCEPT
    {
        return enoki::concat( enoki::low(ax), enoki::low(by) );
    }

    static inline Vec4d cathi(const Vec4d ax, const Vec4d by) NOEXCEPT force_inline;
    static inline Vec4d cathi(const Vec4d ax, const Vec4d by) NOEXCEPT
    {
        return enoki::concat( enoki::high(ax), enoki::high(by) );
    }

} // namespace OTFFT_MISC

//=============================================================================
// 512 bit Vector Emulation
//=============================================================================

namespace OTFFT_MISC {

    //struct Vec8d { Vec4d lo, hi; };
    using Vec8d = enoki::Array<double, 8>;

    static inline Vec8d getez4(const_complex_vector a) NOEXCEPT force_inline;
    static inline Vec8d getez4(const_complex_vector a) NOEXCEPT
    {
        #ifdef USE_UNALIGNED_MEMORY
            return enoki::load_unaligned<Vec8d>(&a->Re);
        #else
            return enoki::load<Vec8d>(&a->Re);
        #endif
    }

    static inline void setez4(complex_vector a, const Vec8d& z) NOEXCEPT force_inline3;
    static inline void setez4(complex_vector a, const Vec8d& z) NOEXCEPT
    {
        #ifdef USE_UNALIGNED_MEMORY
            enoki::store_unaligned<Vec8d>(&a->Re, z);
        #else
            enoki::store<Vec8d>(&a->Re, z);
        #endif
    }

    static inline Vec8d jxez4(const Vec8d& xy) NOEXCEPT force_inline;
    static inline Vec8d jxez4(const Vec8d& xy) NOEXCEPT
    {
        return enoki::concat( jxpz2(enoki::low(xy)), jxpz2(enoki::high(xy)) );
    }

    static inline Vec8d addez4(const Vec8d& a, const Vec8d& b) NOEXCEPT force_inline;
    static inline Vec8d addez4(const Vec8d& a, const Vec8d& b) NOEXCEPT
    {
        return enoki::concat( addpz2(enoki::low(a), enoki::low(b)), addpz2(enoki::high(a), enoki::high(b)) );
    }
    static inline Vec8d subez4(const Vec8d& a, const Vec8d& b) NOEXCEPT force_inline;
    static inline Vec8d subez4(const Vec8d& a, const Vec8d& b) NOEXCEPT
    {
        return enoki::concat( subpz2(enoki::low(a), enoki::low(b)), subpz2(enoki::high(a), enoki::high(b)) );
    }

    static inline Vec8d mulez4(const Vec8d& a, const Vec8d& b) NOEXCEPT force_inline;
    static inline Vec8d mulez4(const Vec8d& a, const Vec8d& b) NOEXCEPT
    {
        return enoki::concat( mulpz2(enoki::low(a), enoki::low(b)), mulpz2(enoki::high(a), enoki::high(b)) );
    }

    static inline Vec8d v8xez4(const Vec8d& xy) NOEXCEPT force_inline;
    static inline Vec8d v8xez4(const Vec8d& xy) NOEXCEPT
    {
        return enoki::concat( v8xpz2(enoki::low(xy)), v8xpz2(enoki::high(xy)) );
    }

    static inline Vec8d w8xez4(const Vec8d& xy) NOEXCEPT force_inline;
    static inline Vec8d w8xez4(const Vec8d& xy) NOEXCEPT
    {
        return enoki::concat( w8xpz2(enoki::low(xy)), w8xpz2(enoki::high(xy)) );
    }

    static inline Vec8d dupez4(const Vec2d x) NOEXCEPT force_inline;
    static inline Vec8d dupez4(const Vec2d x) NOEXCEPT
    {
        const Vec4d y = duppz2(x);
        return enoki::concat( y, y );
    }

    static inline Vec8d dupez5(const complex_t& z) NOEXCEPT force_inline;
    static inline Vec8d dupez5(const complex_t& z) NOEXCEPT
    {
        return dupez4(getpz(z));
    }

} // namespace OTFFT_MISC

//=============================================================================

#endif // otfft_misc_h 
