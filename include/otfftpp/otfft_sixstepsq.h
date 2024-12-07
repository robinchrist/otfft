/******************************************************************************
*  OTFFT SixStep of Square Version 11.5e
*
*  Copyright (c) 2019 OK Ojisan(Takuya OKAHISA)
*  Released under the MIT license
*  http://opensource.org/licenses/mit-license.php
******************************************************************************/

#ifndef otfft_sixstepsq_h
#define otfft_sixstepsq_h

#include "otfft_avxdif8.h"
namespace OTFFT_FourStep = OTFFT_AVXDIF8;

namespace OTFFT_SixStep { /////////////////////////////////////////////////////

template <int log_N, int mode, bool sng = 0> struct fwdffts
{
    static constexpr int log_n = log_N/2;
    static constexpr int N = 1 << log_N;
    static constexpr int n = 1 << log_n;
    static constexpr int m = n/2;
    static constexpr int L = m*(m+1)/2;

    static void transpose_kernel(const int k, const int p, complex_vector x) noexcept
    {
        if (k != p) {
            const int p_kn = p + k*n;
            const int k_pn = k + p*n;
            const Vec4d ab = getpz2(x+p_kn+0);
            const Vec4d AB = getpz2(x+p_kn+n);
            const Vec4d aA = catlo(ab, AB);
            const Vec4d bB = cathi(ab, AB);
            const Vec4d cC = getpz2(x+k_pn+0);
            const Vec4d dD = getpz2(x+k_pn+n);
            setpz2(x+k_pn+n, bB);
            setpz2(x+k_pn+0, aA);
            const Vec4d cd = catlo(cC, dD);
            const Vec4d CD = cathi(cC, dD);
            setpz2(x+p_kn+n, CD);
            setpz2(x+p_kn+0, cd);
            return;
        }
        else {
            const int k_kn = k + k*n;
            const Vec2d x12 = getpz(x[k_kn+1]);
            const Vec2d x21 = getpz(x[k_kn+n]);
            setpz(x[k_kn+n], x12);
            setpz(x[k_kn+1], x21);
            return;
        }
    }

    static void mult_twiddle_factor_kernel(const int i,
            const int p, const int k, complex_vector x, weight_t W) noexcept
    {
        if (k != p) {
            //const Vec2d w11 = getpz(W[k+0 + p*n+0]);
            //const Vec2d w12 = getpz(W[k+1 + p*n+0]);
            //const Vec2d w21 = getpz(W[k+0 + p*n+n]);
            //const Vec2d w22 = getpz(W[k+1 + p*n+n]);
            //const Vec4d w1 = scalepz2<N,mode>(cat(w11, w12));
            //const Vec4d w2 = scalepz2<N,mode>(cat(w21, w22));
            const Vec4d w1 = scalepz2<N,mode>(getpz2(W+4*i+0));
            const Vec4d w2 = scalepz2<N,mode>(getpz2(W+4*i+2));
            const int k_pn = k + p*n;
            const int p_kn = p + k*n;
            const Vec4d ab = mulpz2(w1, getpz2(x+k_pn+0));
            const Vec4d AB = mulpz2(w2, getpz2(x+k_pn+n));
            const Vec4d aA = catlo(ab, AB);
            const Vec4d bB = cathi(ab, AB);
            const Vec4d cC = getpz2(x+p_kn+0);
            const Vec4d dD = getpz2(x+p_kn+n);
            setpz2(x+p_kn+n, bB);
            setpz2(x+p_kn+0, aA);
            const Vec4d cd = mulpz2(w1, catlo(cC, dD));
            const Vec4d CD = mulpz2(w2, cathi(cC, dD));
            setpz2(x+k_pn+n, CD);
            setpz2(x+k_pn+0, cd);
            return;
        }
        else {
            //const Vec2d w11 = getpz(W[p+0 + p*n+0]);
            //const Vec2d w12 = getpz(W[p+1 + p*n+0]);
            //const Vec2d w22 = getpz(W[p+1 + p*n+n]);
            //const Vec4d w1 = scalepz2<N,mode>(cat(w11, w12));
            //const Vec4d w2 = scalepz2<N,mode>(cat(w12, w22));
            const Vec4d w1 = scalepz2<N,mode>(getpz2(W+4*i+0));
            const Vec4d w2 = scalepz2<N,mode>(getpz2(W+4*i+2));
            const int p_pn = p + p*n;
            const Vec4d ab = mulpz2(w1, getpz2(x+p_pn+0));
            const Vec4d AB = mulpz2(w2, getpz2(x+p_pn+n));
            const Vec4d aA = catlo(ab, AB);
            const Vec4d bB = cathi(ab, AB);
            setpz2(x+p_pn+n, bB);
            setpz2(x+p_pn+0, aA);
            return;
        }
    }

    void operator()(const_index_vector iv,
            complex_vector x, complex_vector y, weight_t W, weight_t Ws) const noexcept
    {
        part1(iv, x, y, W, Ws);
    }

    static void part1(const_index_vector iv,
            complex_vector x, complex_vector y, weight_t W, weight_t Ws) noexcept
    {
        for (int i = 0; i < L; i++) {
            const int k = iv[i].row;
            const int p = iv[i].col;
            transpose_kernel(k, p, x);
        }
        for (int p = 0; p < n; p++) {
            const int pn = p*n;
            OTFFT_FourStep::fwdfft<n,1,0,scale_1>()(x + pn, y + pn, Ws);
        }
        for (int i = 0; i < L; i++) {
            const int k = iv[i].col;
            const int p = iv[i].row;
            mult_twiddle_factor_kernel(i, p, k, x, W);
        }
        for (int k = 0; k < n; k++) {
            const int kn = k*n;
            OTFFT_FourStep::fwdfft<n,1,0,scale_1>()(x + kn, y + kn, Ws);
        }
        for (int i = 0; i < L; i++) {
            const int k = iv[i].row;
            const int p = iv[i].col;
            transpose_kernel(k, p, x);
        }
    }

};

///////////////////////////////////////////////////////////////////////////////

template <int log_N, int mode, bool sng = 0> struct invffts
{
    static constexpr int log_n = log_N/2;
    static constexpr int N = 1 << log_N;
    static constexpr int n = 1 << log_n;
    static constexpr int m = n/2;
    static constexpr int L = m*(m+1)/2;
    
    static inline void transpose_kernel(
            const int k, const int p, complex_vector x) noexcept
    {
        fwdffts<log_N,mode,sng>::transpose_kernel(k, p, x);
    }

    static void mult_twiddle_factor_kernel(const int i,
            const int p, const int k, complex_vector x, weight_t W) noexcept
    {
        if (k != p) {
            //const Vec2d w11 = getpz(W[k+0 + p*n+0]);
            //const Vec2d w12 = getpz(W[k+1 + p*n+0]);
            //const Vec2d w21 = getpz(W[k+0 + p*n+n]);
            //const Vec2d w22 = getpz(W[k+1 + p*n+n]);
            //const Vec4d w1 = scalepz2<N,mode>(cnjpz2(cat(w11, w12)));
            //const Vec4d w2 = scalepz2<N,mode>(cnjpz2(cat(w21, w22)));
            const Vec4d w1 = scalepz2<N,mode>(cnjpz2(getpz2(W+4*i+0)));
            const Vec4d w2 = scalepz2<N,mode>(cnjpz2(getpz2(W+4*i+2)));
            const int k_pn = k + p*n;
            const int p_kn = p + k*n;
            const Vec4d ab = mulpz2(w1, getpz2(x+k_pn+0));
            const Vec4d AB = mulpz2(w2, getpz2(x+k_pn+n));
            const Vec4d aA = catlo(ab, AB);
            const Vec4d bB = cathi(ab, AB);
            const Vec4d dD = getpz2(x+p_kn+n);
            const Vec4d cC = getpz2(x+p_kn+0);
            setpz2(x+p_kn+0, aA);
            setpz2(x+p_kn+n, bB);
            const Vec4d cd = mulpz2(w1, catlo(cC, dD));
            const Vec4d CD = mulpz2(w2, cathi(cC, dD));
            setpz2(x+k_pn+n, CD);
            setpz2(x+k_pn+0, cd);
            return;
        }
        else {
            //const Vec2d w11 = getpz(W[p+0 + p*n+0]);
            //const Vec2d w12 = getpz(W[p+1 + p*n+0]);
            //const Vec2d w22 = getpz(W[p+1 + p*n+n]);
            //const Vec4d w1 = scalepz2<N,mode>(cnjpz2(cat(w11, w12)));
            //const Vec4d w2 = scalepz2<N,mode>(cnjpz2(cat(w12, w22)));
            const Vec4d w1 = scalepz2<N,mode>(cnjpz2(getpz2(W+4*i+0)));
            const Vec4d w2 = scalepz2<N,mode>(cnjpz2(getpz2(W+4*i+2)));
            const int p_pn = p + p*n;
            const Vec4d ab = mulpz2(w1, getpz2(x+p_pn+0));
            const Vec4d AB = mulpz2(w2, getpz2(x+p_pn+n));
            const Vec4d aA = catlo(ab, AB);
            const Vec4d bB = cathi(ab, AB);
            setpz2(x+p_pn+n, bB);
            setpz2(x+p_pn+0, aA);
            return;
        }
    }

    void operator()(const_index_vector iv,
            complex_vector x, complex_vector y, weight_t W, weight_t Ws) const noexcept
    {
        part1(iv, x, y, W, Ws);
    }

    static void part1(const_index_vector iv,
            complex_vector x, complex_vector y, weight_t W, weight_t Ws) noexcept
    {
        for (int i = 0; i < L; i++) {
            const int k = iv[i].row;
            const int p = iv[i].col;
            transpose_kernel(k, p, x);
        }
        for (int p = 0; p < n; p++) {
            const int pn = p*n;
            OTFFT_FourStep::invfft<n,1,0,scale_1>()(x + pn, y + pn, Ws);
        }
        for (int i = 0; i < L; i++) {
            const int k = iv[i].col;
            const int p = iv[i].row;
            mult_twiddle_factor_kernel(i, p, k, x, W);
        }
        for (int k = 0; k < n; k++) {
            const int kn = k*n;
            OTFFT_FourStep::invfft<n,1,0,scale_1>()(x + kn, y + kn, Ws);
        }
        for (int i = 0; i < L; i++) {
            const int k = iv[i].row;
            const int p = iv[i].col;
            transpose_kernel(k, p, x);
        }
    }

};

} /////////////////////////////////////////////////////////////////////////////

#endif // otfft_sixstepsq_h
