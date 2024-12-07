/******************************************************************************
*  OTFFT AVXDIF(Radix-4) Version 11.5e
*
*  Copyright (c) 2019 OK Ojisan(Takuya OKAHISA)
*  Released under the MIT license
*  http://opensource.org/licenses/mit-license.php
******************************************************************************/

#ifndef otfft_avxdif4_h
#define otfft_avxdif4_h

#include "otfft_misc.h"

namespace OTFFT_AVXDIF4 { /////////////////////////////////////////////////////

using namespace OTFFT_MISC;

///////////////////////////////////////////////////////////////////////////////
// Forward Butterfly Operation
///////////////////////////////////////////////////////////////////////////////

template <int n, int s> struct fwdcore
{
    static constexpr int m  = n/4;
    static constexpr int N  = n*s;
    static constexpr int N0 = 0;
    static constexpr int N1 = N/4;
    static constexpr int N2 = N1*2;
    static constexpr int N3 = N1*3;

    void operator()(
            complex_vector x, complex_vector y, const_complex_vector W) const noexcept
    {
        for (int p = 0; p < m; p++) {
            const int sp = s*p;
            const int s4p = 4*sp;
            const Vec8d w1p = dupez5(*twidT<4,N,1>(W,sp));
            const Vec8d w2p = dupez5(*twidT<4,N,2>(W,sp));
            const Vec8d w3p = dupez5(*twidT<4,N,3>(W,sp));
            for (int q = 0; q < s; q += 4) {
                complex_vector xq_sp  = x + q + sp;
                complex_vector yq_s4p = y + q + s4p;

                const Vec8d a = getez4(xq_sp+N0);
                const Vec8d c = getez4(xq_sp+N2);
                const Vec8d  apc =       addez4(a, c);
                const Vec8d  amc =       subez4(a, c);
                const Vec8d b = getez4(xq_sp+N1);
                const Vec8d d = getez4(xq_sp+N3);
                const Vec8d  bpd =       addez4(b, d);
                const Vec8d jbmd = jxez4(subez4(b, d));

                setez4(yq_s4p+s*0,             addez4(apc,  bpd));
                setez4(yq_s4p+s*1, mulez4(w1p, subez4(amc, jbmd)));
                setez4(yq_s4p+s*2, mulez4(w2p, subez4(apc,  bpd)));
                setez4(yq_s4p+s*3, mulez4(w3p, addez4(amc, jbmd)));
            }
        }
    }
};

template <int N> struct fwdcore<N,1>
{
    static constexpr int N0 = 0;
    static constexpr int N1 = N/4;
    static constexpr int N2 = N1*2;
    static constexpr int N3 = N1*3;

    void operator()(
            complex_vector x, complex_vector y, const_complex_vector W) const noexcept
    {
        for (int p = 0; p < N1; p += 2) {
            complex_vector x_p  = x + p;
            complex_vector y_4p = y + 4*p;

            const Vec4d a = getpz2(x_p+N0);
            const Vec4d c = getpz2(x_p+N2);
            const Vec4d  apc =       addpz2(a, c);
            const Vec4d  amc =       subpz2(a, c);
            const Vec4d b = getpz2(x_p+N1);
            const Vec4d d = getpz2(x_p+N3);
            const Vec4d  bpd =       addpz2(b, d);
            const Vec4d jbmd = jxpz2(subpz2(b, d));

            const Vec4d w1p = getpz2(twid<4,N,1>(W,p));
            const Vec4d w2p = getpz2(twid<4,N,2>(W,p));
            const Vec4d w3p = getpz2(twid<4,N,3>(W,p));

            const Vec4d aA =             addpz2(apc,  bpd);
            const Vec4d bB = mulpz2(w1p, subpz2(amc, jbmd));
            const Vec4d cC = mulpz2(w2p, subpz2(apc,  bpd));
            const Vec4d dD = mulpz2(w3p, addpz2(amc, jbmd));
            const Vec4d ab = catlo(aA, bB);
            setpz2(y_4p+0, ab);
            const Vec4d cd = catlo(cC, dD);
            setpz2(y_4p+2, cd);
            const Vec4d AB = cathi(aA, bB);
            setpz2(y_4p+4, AB);
            const Vec4d CD = cathi(cC, dD);
            setpz2(y_4p+6, CD);

        }
    }
};

///////////////////////////////////////////////////////////////////////////////

template <int n, int s, bool eo, int mode> struct fwdend;

//-----------------------------------------------------------------------------

template <int s, bool eo, int mode> struct fwdend<4,s,eo,mode>
{
    static constexpr int N = 4*s;

    void operator()(complex_vector x, complex_vector y) const noexcept
    {
        complex_vector z = eo ? y : x;
        for (int q = 0; q < s; q += 2) {
            complex_vector xq = x + q;
            complex_vector zq = z + q;
            const Vec4d a = scalepz2<N,mode>(getpz2(xq+s*0));
            const Vec4d b = scalepz2<N,mode>(getpz2(xq+s*1));
            const Vec4d c = scalepz2<N,mode>(getpz2(xq+s*2));
            const Vec4d d = scalepz2<N,mode>(getpz2(xq+s*3));
            const Vec4d  apc =       addpz2(a, c);
            const Vec4d  amc =       subpz2(a, c);
            const Vec4d  bpd =       addpz2(b, d);
            const Vec4d jbmd = jxpz2(subpz2(b, d));
            setpz2(zq+s*0, addpz2(apc,  bpd));
            setpz2(zq+s*1, subpz2(amc, jbmd));
            setpz2(zq+s*2, subpz2(apc,  bpd));
            setpz2(zq+s*3, addpz2(amc, jbmd));
        }
    }
};

template <bool eo, int mode> struct fwdend<4,1,eo,mode>
{
    inline void operator()(complex_vector x, complex_vector y) const noexcept
    {
        complex_vector z = eo ? y : x;
        const Vec2d a = scalepz<4,mode>(getpz(x[0]));
        const Vec2d b = scalepz<4,mode>(getpz(x[1]));
        const Vec2d c = scalepz<4,mode>(getpz(x[2]));
        const Vec2d d = scalepz<4,mode>(getpz(x[3]));
        const Vec2d  apc =      addpz(a, c);
        const Vec2d  amc =      subpz(a, c);
        const Vec2d  bpd =      addpz(b, d);
        const Vec2d jbmd = jxpz(subpz(b, d));
        setpz(z[0], addpz(apc,  bpd));
        setpz(z[1], subpz(amc, jbmd));
        setpz(z[2], subpz(apc,  bpd));
        setpz(z[3], addpz(amc, jbmd));
    }
};

//-----------------------------------------------------------------------------

template <int s, bool eo, int mode> struct fwdend<2,s,eo,mode>
{
    static constexpr int N = 2*s;

    void operator()(complex_vector x, complex_vector y) const noexcept
    {
        complex_vector z = eo ? y : x;
        for (int q = 0; q < s; q += 2) {
            complex_vector xq = x + q;
            complex_vector zq = z + q;
            const Vec4d a = scalepz2<N,mode>(getpz2(xq+0));
            const Vec4d b = scalepz2<N,mode>(getpz2(xq+s));
            setpz2(zq+0, addpz2(a, b));
            setpz2(zq+s, subpz2(a, b));
        }
    }
};

template <bool eo, int mode> struct fwdend<2,1,eo,mode>
{
    inline void operator()(complex_vector x, complex_vector y) const noexcept
    {
        complex_vector z = eo ? y : x;
        const Vec2d a = scalepz<2,mode>(getpz(x[0]));
        const Vec2d b = scalepz<2,mode>(getpz(x[1]));
        setpz(z[0], addpz(a, b));
        setpz(z[1], subpz(a, b));
    }
};

///////////////////////////////////////////////////////////////////////////////
// Forward FFT
///////////////////////////////////////////////////////////////////////////////

template <int n, int s, bool eo, int mode> struct fwdfft
{
    inline void operator()(
            complex_vector x, complex_vector y, const_complex_vector W) const noexcept
    {
        fwdcore<n,s>()(x, y, W);
        fwdfft<n/4,4*s,!eo,mode>()(y, x, W);
    }
};

template <int s, bool eo, int mode> struct fwdfft<4,s,eo,mode>
{
    inline void operator()(
            complex_vector x, complex_vector y, const_complex_vector) const noexcept
    {
        fwdend<4,s,eo,mode>()(x, y);
    }
};

template <int s, bool eo, int mode> struct fwdfft<2,s,eo,mode>
{
    inline void operator()(
            complex_vector x, complex_vector y, const_complex_vector) const noexcept
    {
        fwdend<2,s,eo,mode>()(x, y);
    }
};

///////////////////////////////////////////////////////////////////////////////
// Inverse Butterfly Operation
///////////////////////////////////////////////////////////////////////////////

template <int n, int s> struct invcore
{
    static constexpr int m  = n/4;
    static constexpr int N  = n*s;
    static constexpr int N0 = 0;
    static constexpr int N1 = N/4;
    static constexpr int N2 = N1*2;
    static constexpr int N3 = N1*3;

    void operator()(
            complex_vector x, complex_vector y, const_complex_vector W) const noexcept
    {
        for (int p = 0; p < m; p++) {
            const int sp = s*p;
            const int s4p = 4*sp;

            const Vec8d w1p = dupez5(conj(*twidT<4,N,1>(W,sp)));
            const Vec8d w2p = dupez5(conj(*twidT<4,N,2>(W,sp)));
            const Vec8d w3p = dupez5(conj(*twidT<4,N,3>(W,sp)));

            for (int q = 0; q < s; q += 4) {
                complex_vector xq_sp  = x + q + sp;
                complex_vector yq_s4p = y + q + s4p;

                const Vec8d a = getez4(xq_sp+N0);
                const Vec8d c = getez4(xq_sp+N2);
                const Vec8d  apc =       addez4(a, c);
                const Vec8d  amc =       subez4(a, c);
                const Vec8d b = getez4(xq_sp+N1);
                const Vec8d d = getez4(xq_sp+N3);
                const Vec8d  bpd =       addez4(b, d);
                const Vec8d jbmd = jxez4(subez4(b, d));

                setez4(yq_s4p+s*0,             addez4(apc,  bpd));
                setez4(yq_s4p+s*1, mulez4(w1p, addez4(amc, jbmd)));
                setez4(yq_s4p+s*2, mulez4(w2p, subez4(apc,  bpd)));
                setez4(yq_s4p+s*3, mulez4(w3p, subez4(amc, jbmd)));
            }
        }
    }
};

template <int N> struct invcore<N,1>
{
    static constexpr int N0 = 0;
    static constexpr int N1 = N/4;
    static constexpr int N2 = N1*2;
    static constexpr int N3 = N1*3;

    void operator()(
            complex_vector x, complex_vector y, const_complex_vector W) const noexcept
    {
        for (int p = 0; p < N1; p += 2) {
            complex_vector x_p  = x + p;
            complex_vector y_4p = y + 4*p;

            const Vec4d a = getpz2(x_p+N0);
            const Vec4d c = getpz2(x_p+N2);
            const Vec4d  apc =       addpz2(a, c);
            const Vec4d  amc =       subpz2(a, c);
            const Vec4d b = getpz2(x_p+N1);
            const Vec4d d = getpz2(x_p+N3);
            const Vec4d  bpd =       addpz2(b, d);
            const Vec4d jbmd = jxpz2(subpz2(b, d));

            const Vec4d w1p = cnjpz2(getpz2(twid<4,N,1>(W,p)));
            const Vec4d w2p = cnjpz2(getpz2(twid<4,N,2>(W,p)));
            const Vec4d w3p = cnjpz2(getpz2(twid<4,N,3>(W,p)));

            const Vec4d aA =             addpz2(apc,  bpd);
            const Vec4d bB = mulpz2(w1p, addpz2(amc, jbmd));
            const Vec4d cC = mulpz2(w2p, subpz2(apc,  bpd));
            const Vec4d dD = mulpz2(w3p, subpz2(amc, jbmd));
            const Vec4d ab = catlo(aA, bB);
            setpz2(y_4p+0, ab);
            const Vec4d cd = catlo(cC, dD);
            setpz2(y_4p+2, cd);
            const Vec4d AB = cathi(aA, bB);
            setpz2(y_4p+4, AB);
            const Vec4d CD = cathi(cC, dD);
            setpz2(y_4p+6, CD);

        }
    }
};

///////////////////////////////////////////////////////////////////////////////

template <int n, int s, bool eo, int mode> struct invend;

//-----------------------------------------------------------------------------

template <int s, bool eo, int mode> struct invend<4,s,eo,mode>
{
    static constexpr int N = 4*s;

    void operator()(complex_vector x, complex_vector y) const noexcept
    {
        complex_vector z = eo ? y : x;
        for (int q = 0; q < s; q += 2) {
            complex_vector xq = x + q;
            complex_vector zq = z + q;
            const Vec4d a = scalepz2<N,mode>(getpz2(xq+s*0));
            const Vec4d b = scalepz2<N,mode>(getpz2(xq+s*1));
            const Vec4d c = scalepz2<N,mode>(getpz2(xq+s*2));
            const Vec4d d = scalepz2<N,mode>(getpz2(xq+s*3));
            const Vec4d  apc =       addpz2(a, c);
            const Vec4d  amc =       subpz2(a, c);
            const Vec4d  bpd =       addpz2(b, d);
            const Vec4d jbmd = jxpz2(subpz2(b, d));
            setpz2(zq+s*0, addpz2(apc,  bpd));
            setpz2(zq+s*1, addpz2(amc, jbmd));
            setpz2(zq+s*2, subpz2(apc,  bpd));
            setpz2(zq+s*3, subpz2(amc, jbmd));
        }
    }
};

template <bool eo, int mode> struct invend<4,1,eo,mode>
{
    inline void operator()(complex_vector x, complex_vector y) const noexcept
    {
        complex_vector z = eo ? y : x;
        const Vec2d a = scalepz<4,mode>(getpz(x[0]));
        const Vec2d b = scalepz<4,mode>(getpz(x[1]));
        const Vec2d c = scalepz<4,mode>(getpz(x[2]));
        const Vec2d d = scalepz<4,mode>(getpz(x[3]));
        const Vec2d  apc =      addpz(a, c);
        const Vec2d  amc =      subpz(a, c);
        const Vec2d  bpd =      addpz(b, d);
        const Vec2d jbmd = jxpz(subpz(b, d));
        setpz(z[0], addpz(apc,  bpd));
        setpz(z[1], addpz(amc, jbmd));
        setpz(z[2], subpz(apc,  bpd));
        setpz(z[3], subpz(amc, jbmd));
    }
};

//-----------------------------------------------------------------------------

template <int s, bool eo, int mode> struct invend<2,s,eo,mode>
{
    static constexpr int N = 2*s;

    void operator()(complex_vector x, complex_vector y) const noexcept
    {
        complex_vector z = eo ? y : x;
        for (int q = 0; q < s; q += 2) {
            complex_vector xq = x + q;
            complex_vector zq = z + q;
            const Vec4d a = scalepz2<N,mode>(getpz2(xq+0));
            const Vec4d b = scalepz2<N,mode>(getpz2(xq+s));
            setpz2(zq+0, addpz2(a, b));
            setpz2(zq+s, subpz2(a, b));
        }
    }
};

template <bool eo, int mode> struct invend<2,1,eo,mode>
{
    inline void operator()(complex_vector x, complex_vector y) const noexcept
    {
        complex_vector z = eo ? y : x;
        const Vec2d a = scalepz<2,mode>(getpz(x[0]));
        const Vec2d b = scalepz<2,mode>(getpz(x[1]));
        setpz(z[0], addpz(a, b));
        setpz(z[1], subpz(a, b));
    }
};

///////////////////////////////////////////////////////////////////////////////
// Inverse FFT
///////////////////////////////////////////////////////////////////////////////

template <int n, int s, bool eo, int mode> struct invfft
{
    inline void operator()(
            complex_vector x, complex_vector y, const_complex_vector W) const noexcept
    {
        invcore<n,s>()(x, y, W);
        invfft<n/4,4*s,!eo,mode>()(y, x, W);
    }
};

template <int s, bool eo, int mode> struct invfft<4,s,eo,mode>
{
    inline void operator()(
            complex_vector x, complex_vector y, const_complex_vector) const noexcept
    {
        invend<4,s,eo,mode>()(x, y);
    }
};

template <int s, bool eo, int mode> struct invfft<2,s,eo,mode>
{
    inline void operator()(
            complex_vector x, complex_vector y, const_complex_vector) const noexcept
    {
        invend<2,s,eo,mode>()(x, y);
    }
};

///////////////////////////////////////////////////////////////////////////////
// FFT object
///////////////////////////////////////////////////////////////////////////////

struct FFT0
{
    int N, log_N;
    simd_array<complex_t> weight;
    complex_t* __restrict W;

    FFT0() noexcept : N(0), log_N(0), W(0) {}
    FFT0(const int n) { setup(n); }

    void setup(int n)
    {
        for (log_N = 0; n > 1; n >>= 1) log_N++;
        setup2(log_N);
    }

    inline void setup2(const int n)
    {
        log_N = n; N = 1 << n;
        if (N <= 4) W = 0;
        else {
            weight.setup(2*N);
            W = weight.get();
            init_Wt(4, N, W);
        }
    }

    ///////////////////////////////////////////////////////////////////////////

    void fwd(complex_vector x, complex_vector y) const noexcept
    {
        constexpr int mode = scale_length;
        switch (log_N) {
            case  0: break;
            case  1: fwdfft<(1<< 1),1,0,mode>()(x, y, W); break;
            case  2: fwdfft<(1<< 2),1,0,mode>()(x, y, W); break;
            case  3: fwdfft<(1<< 3),1,0,mode>()(x, y, W); break;
            case  4: fwdfft<(1<< 4),1,0,mode>()(x, y, W); break;
            case  5: fwdfft<(1<< 5),1,0,mode>()(x, y, W); break;
            case  6: fwdfft<(1<< 6),1,0,mode>()(x, y, W); break;
            case  7: fwdfft<(1<< 7),1,0,mode>()(x, y, W); break;
            case  8: fwdfft<(1<< 8),1,0,mode>()(x, y, W); break;
            case  9: fwdfft<(1<< 9),1,0,mode>()(x, y, W); break;
            case 10: fwdfft<(1<<10),1,0,mode>()(x, y, W); break;
            case 11: fwdfft<(1<<11),1,0,mode>()(x, y, W); break;
            case 12: fwdfft<(1<<12),1,0,mode>()(x, y, W); break;
            case 13: fwdfft<(1<<13),1,0,mode>()(x, y, W); break;
            case 14: fwdfft<(1<<14),1,0,mode>()(x, y, W); break;
            case 15: fwdfft<(1<<15),1,0,mode>()(x, y, W); break;
            case 16: fwdfft<(1<<16),1,0,mode>()(x, y, W); break;
            case 17: fwdfft<(1<<17),1,0,mode>()(x, y, W); break;
            case 18: fwdfft<(1<<18),1,0,mode>()(x, y, W); break;
            case 19: fwdfft<(1<<19),1,0,mode>()(x, y, W); break;
            case 20: fwdfft<(1<<20),1,0,mode>()(x, y, W); break;
            case 21: fwdfft<(1<<21),1,0,mode>()(x, y, W); break;
            case 22: fwdfft<(1<<22),1,0,mode>()(x, y, W); break;
            case 23: fwdfft<(1<<23),1,0,mode>()(x, y, W); break;
            case 24: fwdfft<(1<<24),1,0,mode>()(x, y, W); break;
        }
    }

    void fwd0(complex_vector x, complex_vector y) const noexcept
    {
        constexpr int mode = scale_1;
        switch (log_N) {
            case  0: break;
            case  1: fwdfft<(1<< 1),1,0,mode>()(x, y, W); break;
            case  2: fwdfft<(1<< 2),1,0,mode>()(x, y, W); break;
            case  3: fwdfft<(1<< 3),1,0,mode>()(x, y, W); break;
            case  4: fwdfft<(1<< 4),1,0,mode>()(x, y, W); break;
            case  5: fwdfft<(1<< 5),1,0,mode>()(x, y, W); break;
            case  6: fwdfft<(1<< 6),1,0,mode>()(x, y, W); break;
            case  7: fwdfft<(1<< 7),1,0,mode>()(x, y, W); break;
            case  8: fwdfft<(1<< 8),1,0,mode>()(x, y, W); break;
            case  9: fwdfft<(1<< 9),1,0,mode>()(x, y, W); break;
            case 10: fwdfft<(1<<10),1,0,mode>()(x, y, W); break;
            case 11: fwdfft<(1<<11),1,0,mode>()(x, y, W); break;
            case 12: fwdfft<(1<<12),1,0,mode>()(x, y, W); break;
            case 13: fwdfft<(1<<13),1,0,mode>()(x, y, W); break;
            case 14: fwdfft<(1<<14),1,0,mode>()(x, y, W); break;
            case 15: fwdfft<(1<<15),1,0,mode>()(x, y, W); break;
            case 16: fwdfft<(1<<16),1,0,mode>()(x, y, W); break;
            case 17: fwdfft<(1<<17),1,0,mode>()(x, y, W); break;
            case 18: fwdfft<(1<<18),1,0,mode>()(x, y, W); break;
            case 19: fwdfft<(1<<19),1,0,mode>()(x, y, W); break;
            case 20: fwdfft<(1<<20),1,0,mode>()(x, y, W); break;
            case 21: fwdfft<(1<<21),1,0,mode>()(x, y, W); break;
            case 22: fwdfft<(1<<22),1,0,mode>()(x, y, W); break;
            case 23: fwdfft<(1<<23),1,0,mode>()(x, y, W); break;
            case 24: fwdfft<(1<<24),1,0,mode>()(x, y, W); break;
        }
    }

    void fwdu(complex_vector x, complex_vector y) const noexcept
    {
        constexpr int mode = scale_unitary;
        switch (log_N) {
            case  0: break;
            case  1: fwdfft<(1<< 1),1,0,mode>()(x, y, W); break;
            case  2: fwdfft<(1<< 2),1,0,mode>()(x, y, W); break;
            case  3: fwdfft<(1<< 3),1,0,mode>()(x, y, W); break;
            case  4: fwdfft<(1<< 4),1,0,mode>()(x, y, W); break;
            case  5: fwdfft<(1<< 5),1,0,mode>()(x, y, W); break;
            case  6: fwdfft<(1<< 6),1,0,mode>()(x, y, W); break;
            case  7: fwdfft<(1<< 7),1,0,mode>()(x, y, W); break;
            case  8: fwdfft<(1<< 8),1,0,mode>()(x, y, W); break;
            case  9: fwdfft<(1<< 9),1,0,mode>()(x, y, W); break;
            case 10: fwdfft<(1<<10),1,0,mode>()(x, y, W); break;
            case 11: fwdfft<(1<<11),1,0,mode>()(x, y, W); break;
            case 12: fwdfft<(1<<12),1,0,mode>()(x, y, W); break;
            case 13: fwdfft<(1<<13),1,0,mode>()(x, y, W); break;
            case 14: fwdfft<(1<<14),1,0,mode>()(x, y, W); break;
            case 15: fwdfft<(1<<15),1,0,mode>()(x, y, W); break;
            case 16: fwdfft<(1<<16),1,0,mode>()(x, y, W); break;
            case 17: fwdfft<(1<<17),1,0,mode>()(x, y, W); break;
            case 18: fwdfft<(1<<18),1,0,mode>()(x, y, W); break;
            case 19: fwdfft<(1<<19),1,0,mode>()(x, y, W); break;
            case 20: fwdfft<(1<<20),1,0,mode>()(x, y, W); break;
            case 21: fwdfft<(1<<21),1,0,mode>()(x, y, W); break;
            case 22: fwdfft<(1<<22),1,0,mode>()(x, y, W); break;
            case 23: fwdfft<(1<<23),1,0,mode>()(x, y, W); break;
            case 24: fwdfft<(1<<24),1,0,mode>()(x, y, W); break;
        }
    }

    inline void fwdn(complex_vector x, complex_vector y) const noexcept
    {
        fwd(x, y);
    }

    ///////////////////////////////////////////////////////////////////////////

    void inv(complex_vector x, complex_vector y) const noexcept
    {
        constexpr int mode = scale_1;
        switch (log_N) {
            case  0: break;
            case  1: invfft<(1<< 1),1,0,mode>()(x, y, W); break;
            case  2: invfft<(1<< 2),1,0,mode>()(x, y, W); break;
            case  3: invfft<(1<< 3),1,0,mode>()(x, y, W); break;
            case  4: invfft<(1<< 4),1,0,mode>()(x, y, W); break;
            case  5: invfft<(1<< 5),1,0,mode>()(x, y, W); break;
            case  6: invfft<(1<< 6),1,0,mode>()(x, y, W); break;
            case  7: invfft<(1<< 7),1,0,mode>()(x, y, W); break;
            case  8: invfft<(1<< 8),1,0,mode>()(x, y, W); break;
            case  9: invfft<(1<< 9),1,0,mode>()(x, y, W); break;
            case 10: invfft<(1<<10),1,0,mode>()(x, y, W); break;
            case 11: invfft<(1<<11),1,0,mode>()(x, y, W); break;
            case 12: invfft<(1<<12),1,0,mode>()(x, y, W); break;
            case 13: invfft<(1<<13),1,0,mode>()(x, y, W); break;
            case 14: invfft<(1<<14),1,0,mode>()(x, y, W); break;
            case 15: invfft<(1<<15),1,0,mode>()(x, y, W); break;
            case 16: invfft<(1<<16),1,0,mode>()(x, y, W); break;
            case 17: invfft<(1<<17),1,0,mode>()(x, y, W); break;
            case 18: invfft<(1<<18),1,0,mode>()(x, y, W); break;
            case 19: invfft<(1<<19),1,0,mode>()(x, y, W); break;
            case 20: invfft<(1<<20),1,0,mode>()(x, y, W); break;
            case 21: invfft<(1<<21),1,0,mode>()(x, y, W); break;
            case 22: invfft<(1<<22),1,0,mode>()(x, y, W); break;
            case 23: invfft<(1<<23),1,0,mode>()(x, y, W); break;
            case 24: invfft<(1<<24),1,0,mode>()(x, y, W); break;
        }
    }

    inline void inv0(complex_vector x, complex_vector y) const noexcept
    {
        inv(x, y);
    }

    void invu(complex_vector x, complex_vector y) const noexcept
    {
        constexpr int mode = scale_unitary;
        switch (log_N) {
            case  0: break;
            case  1: invfft<(1<< 1),1,0,mode>()(x, y, W); break;
            case  2: invfft<(1<< 2),1,0,mode>()(x, y, W); break;
            case  3: invfft<(1<< 3),1,0,mode>()(x, y, W); break;
            case  4: invfft<(1<< 4),1,0,mode>()(x, y, W); break;
            case  5: invfft<(1<< 5),1,0,mode>()(x, y, W); break;
            case  6: invfft<(1<< 6),1,0,mode>()(x, y, W); break;
            case  7: invfft<(1<< 7),1,0,mode>()(x, y, W); break;
            case  8: invfft<(1<< 8),1,0,mode>()(x, y, W); break;
            case  9: invfft<(1<< 9),1,0,mode>()(x, y, W); break;
            case 10: invfft<(1<<10),1,0,mode>()(x, y, W); break;
            case 11: invfft<(1<<11),1,0,mode>()(x, y, W); break;
            case 12: invfft<(1<<12),1,0,mode>()(x, y, W); break;
            case 13: invfft<(1<<13),1,0,mode>()(x, y, W); break;
            case 14: invfft<(1<<14),1,0,mode>()(x, y, W); break;
            case 15: invfft<(1<<15),1,0,mode>()(x, y, W); break;
            case 16: invfft<(1<<16),1,0,mode>()(x, y, W); break;
            case 17: invfft<(1<<17),1,0,mode>()(x, y, W); break;
            case 18: invfft<(1<<18),1,0,mode>()(x, y, W); break;
            case 19: invfft<(1<<19),1,0,mode>()(x, y, W); break;
            case 20: invfft<(1<<20),1,0,mode>()(x, y, W); break;
            case 21: invfft<(1<<21),1,0,mode>()(x, y, W); break;
            case 22: invfft<(1<<22),1,0,mode>()(x, y, W); break;
            case 23: invfft<(1<<23),1,0,mode>()(x, y, W); break;
            case 24: invfft<(1<<24),1,0,mode>()(x, y, W); break;
        }
    }

    void invn(complex_vector x, complex_vector y) const noexcept
    {
        constexpr int mode = scale_length;
        switch (log_N) {
            case  0: break;
            case  1: invfft<(1<< 1),1,0,mode>()(x, y, W); break;
            case  2: invfft<(1<< 2),1,0,mode>()(x, y, W); break;
            case  3: invfft<(1<< 3),1,0,mode>()(x, y, W); break;
            case  4: invfft<(1<< 4),1,0,mode>()(x, y, W); break;
            case  5: invfft<(1<< 5),1,0,mode>()(x, y, W); break;
            case  6: invfft<(1<< 6),1,0,mode>()(x, y, W); break;
            case  7: invfft<(1<< 7),1,0,mode>()(x, y, W); break;
            case  8: invfft<(1<< 8),1,0,mode>()(x, y, W); break;
            case  9: invfft<(1<< 9),1,0,mode>()(x, y, W); break;
            case 10: invfft<(1<<10),1,0,mode>()(x, y, W); break;
            case 11: invfft<(1<<11),1,0,mode>()(x, y, W); break;
            case 12: invfft<(1<<12),1,0,mode>()(x, y, W); break;
            case 13: invfft<(1<<13),1,0,mode>()(x, y, W); break;
            case 14: invfft<(1<<14),1,0,mode>()(x, y, W); break;
            case 15: invfft<(1<<15),1,0,mode>()(x, y, W); break;
            case 16: invfft<(1<<16),1,0,mode>()(x, y, W); break;
            case 17: invfft<(1<<17),1,0,mode>()(x, y, W); break;
            case 18: invfft<(1<<18),1,0,mode>()(x, y, W); break;
            case 19: invfft<(1<<19),1,0,mode>()(x, y, W); break;
            case 20: invfft<(1<<20),1,0,mode>()(x, y, W); break;
            case 21: invfft<(1<<21),1,0,mode>()(x, y, W); break;
            case 22: invfft<(1<<22),1,0,mode>()(x, y, W); break;
            case 23: invfft<(1<<23),1,0,mode>()(x, y, W); break;
            case 24: invfft<(1<<24),1,0,mode>()(x, y, W); break;
        }
    }
};

} /////////////////////////////////////////////////////////////////////////////

#endif // otfft_avxdif4_h
