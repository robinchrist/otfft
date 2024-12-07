/******************************************************************************
*  OTFFT EightStep Version 11.5e
*
*  Copyright (c) 2019 OK Ojisan(Takuya OKAHISA)
*  Released under the MIT license
*  http://opensource.org/licenses/mit-license.php
******************************************************************************/

#ifndef otfft_eightstep_h
#define otfft_eightstep_h

namespace OTFFT_EightStep { ///////////////////////////////////////////////////

using namespace OTFFT_MISC;
using OTFFT_SixStep::weight_t;
using OTFFT_SixStep::const_index_vector;

///////////////////////////////////////////////////////////////////////////////

template <int log_N, int mode> struct fwdfftr
{
    static constexpr int N = 1 << log_N;
    static constexpr int N0 = 0;
    static constexpr int N1 = N/8;
    static constexpr int N2 = N1*2;
    static constexpr int N3 = N1*3;
    static constexpr int N4 = N1*4;
    static constexpr int N5 = N1*5;
    static constexpr int N6 = N1*6;
    static constexpr int N7 = N1*7;

    static inline void transpose_kernel(
            const int p, complex_vector x, complex_vector y) noexcept
    {
        for (int q = 0; q < 8; q += 2) {
            const Vec4d aA = getpz2(x+p+q*N1+N0);
            const Vec4d bB = getpz2(x+p+q*N1+N1);
            const Vec4d ab = catlo(aA, bB);
            const Vec4d AB = cathi(aA, bB);
            setpz2(y+8*p+q+0, ab);
            setpz2(y+8*p+q+8, AB);
        }
    }

    static inline void fft_and_mult_twiddle_factor_kernel(
            const int p, complex_vector x, complex_vector y, weight_t W) noexcept
    {
        const Vec4d x0 = scalepz2<N,mode>(getpz2(x+p+N0));
        const Vec4d x1 = scalepz2<N,mode>(getpz2(x+p+N1));
        const Vec4d x2 = scalepz2<N,mode>(getpz2(x+p+N2));
        const Vec4d x3 = scalepz2<N,mode>(getpz2(x+p+N3));
        const Vec4d x4 = scalepz2<N,mode>(getpz2(x+p+N4));
        const Vec4d x5 = scalepz2<N,mode>(getpz2(x+p+N5));
        const Vec4d x6 = scalepz2<N,mode>(getpz2(x+p+N6));
        const Vec4d x7 = scalepz2<N,mode>(getpz2(x+p+N7));

        const Vec4d  a04 =       addpz2(x0, x4);
        const Vec4d  s04 =       subpz2(x0, x4);
        const Vec4d  a26 =       addpz2(x2, x6);
        const Vec4d js26 = jxpz2(subpz2(x2, x6));
        const Vec4d  a15 =       addpz2(x1, x5);
        const Vec4d  s15 =       subpz2(x1, x5);
        const Vec4d  a37 =       addpz2(x3, x7);
        const Vec4d js37 = jxpz2(subpz2(x3, x7));

        const Vec4d    a04_p1_a26 =        addpz2(a04,  a26);
        const Vec4d    s04_mj_s26 =        subpz2(s04, js26);
        const Vec4d    a04_m1_a26 =        subpz2(a04,  a26);
        const Vec4d    s04_pj_s26 =        addpz2(s04, js26);
        const Vec4d    a15_p1_a37 =        addpz2(a15,  a37);
        const Vec4d w8_s15_mj_s37 = w8xpz2(subpz2(s15, js37));
        const Vec4d  j_a15_m1_a37 =  jxpz2(subpz2(a15,  a37));
        const Vec4d v8_s15_pj_s37 = v8xpz2(addpz2(s15, js37));

        setpz2(y+p+N0,             addpz2(a04_p1_a26,    a15_p1_a37));
        const Vec4d w1p = getpz2(W+p);
        setpz2(y+p+N1, mulpz2(w1p, addpz2(s04_mj_s26, w8_s15_mj_s37)));
        const Vec4d w2p = mulpz2(w1p,w1p);
        setpz2(y+p+N2, mulpz2(w2p, subpz2(a04_m1_a26,  j_a15_m1_a37)));
        const Vec4d w3p = mulpz2(w1p,w2p);
        setpz2(y+p+N3, mulpz2(w3p, subpz2(s04_pj_s26, v8_s15_pj_s37)));
        const Vec4d w4p = mulpz2(w2p,w2p);
        setpz2(y+p+N4, mulpz2(w4p, subpz2(a04_p1_a26,    a15_p1_a37)));
        const Vec4d w5p = mulpz2(w2p,w3p);
        setpz2(y+p+N5, mulpz2(w5p, subpz2(s04_mj_s26, w8_s15_mj_s37)));
        const Vec4d w6p = mulpz2(w3p,w3p);
        setpz2(y+p+N6, mulpz2(w6p, addpz2(a04_m1_a26,  j_a15_m1_a37)));
        const Vec4d w7p = mulpz2(w3p,w4p);
        setpz2(y+p+N7, mulpz2(w7p, addpz2(s04_pj_s26, v8_s15_pj_s37)));

    }

    void operator()(const_index_vector iv,
            complex_vector x, complex_vector y, weight_t W, weight_t Wm, weight_t Ws) const noexcept
    {
        for (int p = 0; p < N1; p += 2) {
            fft_and_mult_twiddle_factor_kernel(p, x, y, W);
        }
        OTFFT_SixStep::fwdffts<log_N-3,scale_1>::part1(iv, y+N0, x+N0, Wm, Ws);
        OTFFT_SixStep::fwdffts<log_N-3,scale_1>::part1(iv, y+N1, x+N1, Wm, Ws);
        OTFFT_SixStep::fwdffts<log_N-3,scale_1>::part1(iv, y+N2, x+N2, Wm, Ws);
        OTFFT_SixStep::fwdffts<log_N-3,scale_1>::part1(iv, y+N3, x+N3, Wm, Ws);
        OTFFT_SixStep::fwdffts<log_N-3,scale_1>::part1(iv, y+N4, x+N4, Wm, Ws);
        OTFFT_SixStep::fwdffts<log_N-3,scale_1>::part1(iv, y+N5, x+N5, Wm, Ws);
        OTFFT_SixStep::fwdffts<log_N-3,scale_1>::part1(iv, y+N6, x+N6, Wm, Ws);
        OTFFT_SixStep::fwdffts<log_N-3,scale_1>::part1(iv, y+N7, x+N7, Wm, Ws);
        for (int p = 0; p < N1; p += 2) {
            transpose_kernel(p, y, x);
        }
    }
};

///////////////////////////////////////////////////////////////////////////////

template <int log_N, int mode> struct invfftr
{
    static constexpr int N = 1 << log_N;
    static constexpr int N0 = 0;
    static constexpr int N1 = N/8;
    static constexpr int N2 = N1*2;
    static constexpr int N3 = N1*3;
    static constexpr int N4 = N1*4;
    static constexpr int N5 = N1*5;
    static constexpr int N6 = N1*6;
    static constexpr int N7 = N1*7;

    static inline void transpose_kernel(
            const int p, complex_vector x, complex_vector y) noexcept
    {
        fwdfftr<log_N,mode>::transpose_kernel(p, x, y);
    }

    static inline void ifft_and_mult_twiddle_factor_kernel(
            const int p, complex_vector x, complex_vector y, weight_t W) noexcept
    {
        const Vec4d x0 = scalepz2<N,mode>(getpz2(x+p+N0));
        const Vec4d x1 = scalepz2<N,mode>(getpz2(x+p+N1));
        const Vec4d x2 = scalepz2<N,mode>(getpz2(x+p+N2));
        const Vec4d x3 = scalepz2<N,mode>(getpz2(x+p+N3));
        const Vec4d x4 = scalepz2<N,mode>(getpz2(x+p+N4));
        const Vec4d x5 = scalepz2<N,mode>(getpz2(x+p+N5));
        const Vec4d x6 = scalepz2<N,mode>(getpz2(x+p+N6));
        const Vec4d x7 = scalepz2<N,mode>(getpz2(x+p+N7));

        const Vec4d  a04 =       addpz2(x0, x4);
        const Vec4d  s04 =       subpz2(x0, x4);
        const Vec4d  a26 =       addpz2(x2, x6);
        const Vec4d js26 = jxpz2(subpz2(x2, x6));
        const Vec4d  a15 =       addpz2(x1, x5);
        const Vec4d  s15 =       subpz2(x1, x5);
        const Vec4d  a37 =       addpz2(x3, x7);
        const Vec4d js37 = jxpz2(subpz2(x3, x7));

        const Vec4d    a04_p1_a26 =        addpz2(a04,  a26);
        const Vec4d    s04_pj_s26 =        addpz2(s04, js26);
        const Vec4d    a04_m1_a26 =        subpz2(a04,  a26);
        const Vec4d    s04_mj_s26 =        subpz2(s04, js26);
        const Vec4d    a15_p1_a37 =        addpz2(a15,  a37);
        const Vec4d v8_s15_pj_s37 = v8xpz2(addpz2(s15, js37));
        const Vec4d  j_a15_m1_a37 =  jxpz2(subpz2(a15,  a37));
        const Vec4d w8_s15_mj_s37 = w8xpz2(subpz2(s15, js37));

        setpz2(y+p+N0,             addpz2(a04_p1_a26,    a15_p1_a37));
        const Vec4d w1p = cnjpz2(getpz2(W+p));
        setpz2(y+p+N1, mulpz2(w1p, addpz2(s04_pj_s26, v8_s15_pj_s37)));
        const Vec4d w2p = mulpz2(w1p,w1p);
        setpz2(y+p+N2, mulpz2(w2p, addpz2(a04_m1_a26,  j_a15_m1_a37)));
        const Vec4d w3p = mulpz2(w1p,w2p);
        setpz2(y+p+N3, mulpz2(w3p, subpz2(s04_mj_s26, w8_s15_mj_s37)));
        const Vec4d w4p = mulpz2(w2p,w2p);
        setpz2(y+p+N4, mulpz2(w4p, subpz2(a04_p1_a26,    a15_p1_a37)));
        const Vec4d w5p = mulpz2(w2p,w3p);
        setpz2(y+p+N5, mulpz2(w5p, subpz2(s04_pj_s26, v8_s15_pj_s37)));
        const Vec4d w6p = mulpz2(w3p,w3p);
        setpz2(y+p+N6, mulpz2(w6p, subpz2(a04_m1_a26,  j_a15_m1_a37)));
        const Vec4d w7p = mulpz2(w3p,w4p);
        setpz2(y+p+N7, mulpz2(w7p, addpz2(s04_mj_s26, w8_s15_mj_s37)));

    }

    void operator()(const_index_vector iv,
            complex_vector x, complex_vector y, weight_t W, weight_t Wm, weight_t Ws) const noexcept
    {
        for (int p = 0; p < N1; p += 2) {
            ifft_and_mult_twiddle_factor_kernel(p, x, y, W);
        }
        OTFFT_SixStep::invffts<log_N-3,scale_1>::part1(iv, y+N0, x+N0, Wm, Ws);
        OTFFT_SixStep::invffts<log_N-3,scale_1>::part1(iv, y+N1, x+N1, Wm, Ws);
        OTFFT_SixStep::invffts<log_N-3,scale_1>::part1(iv, y+N2, x+N2, Wm, Ws);
        OTFFT_SixStep::invffts<log_N-3,scale_1>::part1(iv, y+N3, x+N3, Wm, Ws);
        OTFFT_SixStep::invffts<log_N-3,scale_1>::part1(iv, y+N4, x+N4, Wm, Ws);
        OTFFT_SixStep::invffts<log_N-3,scale_1>::part1(iv, y+N5, x+N5, Wm, Ws);
        OTFFT_SixStep::invffts<log_N-3,scale_1>::part1(iv, y+N6, x+N6, Wm, Ws);
        OTFFT_SixStep::invffts<log_N-3,scale_1>::part1(iv, y+N7, x+N7, Wm, Ws);
        for (int p = 0; p < N1; p += 2) {
            transpose_kernel(p, y, x);
        }
    }
};

} /////////////////////////////////////////////////////////////////////////////

#endif // otfft_eightstep_h
