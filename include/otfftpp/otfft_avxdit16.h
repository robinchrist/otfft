/******************************************************************************
*  OTFFT AVXDIT(Radix-16) Version 11.5e
*
*  Copyright (c) 2019 OK Ojisan(Takuya OKAHISA)
*  Released under the MIT license
*  http://opensource.org/licenses/mit-license.php
******************************************************************************/

#ifndef otfft_avxdit16_h
#define otfft_avxdit16_h

#include "otfft_misc.h"
#include "otfft_avxdit8.h"

namespace OTFFT_AVXDIT16 { ////////////////////////////////////////////////////

using namespace OTFFT_MISC;

///////////////////////////////////////////////////////////////////////////////
// Forward Buffterfly Operation
///////////////////////////////////////////////////////////////////////////////

template <int n, int s> struct fwdcore
{
    static constexpr int m  = n/16;
    static constexpr int N  = n*s;
    static constexpr int N0 = 0;
    static constexpr int N1 = N/16;
    static constexpr int N2 = N1*2;
    static constexpr int N3 = N1*3;
    static constexpr int N4 = N1*4;
    static constexpr int N5 = N1*5;
    static constexpr int N6 = N1*6;
    static constexpr int N7 = N1*7;
    static constexpr int N8 = N1*8;
    static constexpr int N9 = N1*9;
    static constexpr int Na = N1*10;
    static constexpr int Nb = N1*11;
    static constexpr int Nc = N1*12;
    static constexpr int Nd = N1*13;
    static constexpr int Ne = N1*14;
    static constexpr int Nf = N1*15;

    void operator()(
            complex_vector x, complex_vector y, const_complex_vector W) const noexcept
    {
        for (int p = 0; p < m; p++) {
            const int sp   = s*p;
            const int s16p = 16*sp;

            const Vec4d w1p = duppz3(*twidT<16,N, 1>(W,sp));
            const Vec4d w2p = duppz3(*twidT<16,N, 2>(W,sp));
            const Vec4d w3p = duppz3(*twidT<16,N, 3>(W,sp));
            const Vec4d w4p = duppz3(*twidT<16,N, 4>(W,sp));
            const Vec4d w5p = duppz3(*twidT<16,N, 5>(W,sp));
            const Vec4d w6p = duppz3(*twidT<16,N, 6>(W,sp));
            const Vec4d w7p = duppz3(*twidT<16,N, 7>(W,sp));
            const Vec4d w8p = duppz3(*twidT<16,N, 8>(W,sp));
            const Vec4d w9p = duppz3(*twidT<16,N, 9>(W,sp));
            const Vec4d wap = duppz3(*twidT<16,N,10>(W,sp));
            const Vec4d wbp = duppz3(*twidT<16,N,11>(W,sp));
            const Vec4d wcp = duppz3(*twidT<16,N,12>(W,sp));
            const Vec4d wdp = duppz3(*twidT<16,N,13>(W,sp));
            const Vec4d wep = duppz3(*twidT<16,N,14>(W,sp));
            const Vec4d wfp = duppz3(*twidT<16,N,15>(W,sp));

            for (int q = 0; q < s; q += 2) {
                complex_vector xq_sp   = x + q + sp;
                complex_vector yq_s16p = y + q + s16p;

                const Vec4d y0 =             getpz2(yq_s16p+s*0x0);
                const Vec4d y1 = mulpz2(w1p, getpz2(yq_s16p+s*0x1));
                const Vec4d y2 = mulpz2(w2p, getpz2(yq_s16p+s*0x2));
                const Vec4d y3 = mulpz2(w3p, getpz2(yq_s16p+s*0x3));
                const Vec4d y4 = mulpz2(w4p, getpz2(yq_s16p+s*0x4));
                const Vec4d y5 = mulpz2(w5p, getpz2(yq_s16p+s*0x5));
                const Vec4d y6 = mulpz2(w6p, getpz2(yq_s16p+s*0x6));
                const Vec4d y7 = mulpz2(w7p, getpz2(yq_s16p+s*0x7));
                const Vec4d y8 = mulpz2(w8p, getpz2(yq_s16p+s*0x8));
                const Vec4d y9 = mulpz2(w9p, getpz2(yq_s16p+s*0x9));
                const Vec4d ya = mulpz2(wap, getpz2(yq_s16p+s*0xa));
                const Vec4d yb = mulpz2(wbp, getpz2(yq_s16p+s*0xb));
                const Vec4d yc = mulpz2(wcp, getpz2(yq_s16p+s*0xc));
                const Vec4d yd = mulpz2(wdp, getpz2(yq_s16p+s*0xd));
                const Vec4d ye = mulpz2(wep, getpz2(yq_s16p+s*0xe));
                const Vec4d yf = mulpz2(wfp, getpz2(yq_s16p+s*0xf));

                const Vec4d a08 = addpz2(y0, y8); const Vec4d s08 = subpz2(y0, y8);
                const Vec4d a4c = addpz2(y4, yc); const Vec4d s4c = subpz2(y4, yc);
                const Vec4d a2a = addpz2(y2, ya); const Vec4d s2a = subpz2(y2, ya);
                const Vec4d a6e = addpz2(y6, ye); const Vec4d s6e = subpz2(y6, ye);
                const Vec4d a19 = addpz2(y1, y9); const Vec4d s19 = subpz2(y1, y9);
                const Vec4d a5d = addpz2(y5, yd); const Vec4d s5d = subpz2(y5, yd);
                const Vec4d a3b = addpz2(y3, yb); const Vec4d s3b = subpz2(y3, yb);
                const Vec4d a7f = addpz2(y7, yf); const Vec4d s7f = subpz2(y7, yf);

                const Vec4d js4c = jxpz2(s4c);
                const Vec4d js6e = jxpz2(s6e);
                const Vec4d js5d = jxpz2(s5d);
                const Vec4d js7f = jxpz2(s7f);

                const Vec4d a08p1a4c = addpz2(a08, a4c); const Vec4d s08mjs4c = subpz2(s08, js4c);
                const Vec4d a08m1a4c = subpz2(a08, a4c); const Vec4d s08pjs4c = addpz2(s08, js4c);
                const Vec4d a2ap1a6e = addpz2(a2a, a6e); const Vec4d s2amjs6e = subpz2(s2a, js6e);
                const Vec4d a2am1a6e = subpz2(a2a, a6e); const Vec4d s2apjs6e = addpz2(s2a, js6e);
                const Vec4d a19p1a5d = addpz2(a19, a5d); const Vec4d s19mjs5d = subpz2(s19, js5d);
                const Vec4d a19m1a5d = subpz2(a19, a5d); const Vec4d s19pjs5d = addpz2(s19, js5d);
                const Vec4d a3bp1a7f = addpz2(a3b, a7f); const Vec4d s3bmjs7f = subpz2(s3b, js7f);
                const Vec4d a3bm1a7f = subpz2(a3b, a7f); const Vec4d s3bpjs7f = addpz2(s3b, js7f);

                const Vec4d w8_s2amjs6e = w8xpz2(s2amjs6e);
                const Vec4d  j_a2am1a6e =  jxpz2(a2am1a6e);
                const Vec4d v8_s2apjs6e = v8xpz2(s2apjs6e);

                const Vec4d a08p1a4c_p1_a2ap1a6e = addpz2(a08p1a4c,    a2ap1a6e);
                const Vec4d s08mjs4c_pw_s2amjs6e = addpz2(s08mjs4c, w8_s2amjs6e);
                const Vec4d a08m1a4c_mj_a2am1a6e = subpz2(a08m1a4c,  j_a2am1a6e);
                const Vec4d s08pjs4c_mv_s2apjs6e = subpz2(s08pjs4c, v8_s2apjs6e);
                const Vec4d a08p1a4c_m1_a2ap1a6e = subpz2(a08p1a4c,    a2ap1a6e);
                const Vec4d s08mjs4c_mw_s2amjs6e = subpz2(s08mjs4c, w8_s2amjs6e);
                const Vec4d a08m1a4c_pj_a2am1a6e = addpz2(a08m1a4c,  j_a2am1a6e);
                const Vec4d s08pjs4c_pv_s2apjs6e = addpz2(s08pjs4c, v8_s2apjs6e);

                const Vec4d w8_s3bmjs7f = w8xpz2(s3bmjs7f);
                const Vec4d  j_a3bm1a7f =  jxpz2(a3bm1a7f);
                const Vec4d v8_s3bpjs7f = v8xpz2(s3bpjs7f);

                const Vec4d a19p1a5d_p1_a3bp1a7f = addpz2(a19p1a5d,    a3bp1a7f);
                const Vec4d s19mjs5d_pw_s3bmjs7f = addpz2(s19mjs5d, w8_s3bmjs7f);
                const Vec4d a19m1a5d_mj_a3bm1a7f = subpz2(a19m1a5d,  j_a3bm1a7f);
                const Vec4d s19pjs5d_mv_s3bpjs7f = subpz2(s19pjs5d, v8_s3bpjs7f);
                const Vec4d a19p1a5d_m1_a3bp1a7f = subpz2(a19p1a5d,    a3bp1a7f);
                const Vec4d s19mjs5d_mw_s3bmjs7f = subpz2(s19mjs5d, w8_s3bmjs7f);
                const Vec4d a19m1a5d_pj_a3bm1a7f = addpz2(a19m1a5d,  j_a3bm1a7f);
                const Vec4d s19pjs5d_pv_s3bpjs7f = addpz2(s19pjs5d, v8_s3bpjs7f);

                setpz2(xq_sp+N0, addpz2(a08p1a4c_p1_a2ap1a6e,    a19p1a5d_p1_a3bp1a7f));
                setpz2(xq_sp+N8, subpz2(a08p1a4c_p1_a2ap1a6e,    a19p1a5d_p1_a3bp1a7f));

                const Vec4d h1_s19mjs5d_pw_s3bmjs7f = h1xpz2(s19mjs5d_pw_s3bmjs7f);
                setpz2(xq_sp+N1, addpz2(s08mjs4c_pw_s2amjs6e, h1_s19mjs5d_pw_s3bmjs7f));
                setpz2(xq_sp+N9, subpz2(s08mjs4c_pw_s2amjs6e, h1_s19mjs5d_pw_s3bmjs7f));

                const Vec4d w8_a19m1a5d_mj_a3bm1a7f = w8xpz2(a19m1a5d_mj_a3bm1a7f);
                setpz2(xq_sp+N2, addpz2(a08m1a4c_mj_a2am1a6e, w8_a19m1a5d_mj_a3bm1a7f));
                setpz2(xq_sp+Na, subpz2(a08m1a4c_mj_a2am1a6e, w8_a19m1a5d_mj_a3bm1a7f));

                const Vec4d h3_s19pjs5d_mv_s3bpjs7f = h3xpz2(s19pjs5d_mv_s3bpjs7f);
                setpz2(xq_sp+N3, addpz2(s08pjs4c_mv_s2apjs6e, h3_s19pjs5d_mv_s3bpjs7f));
                setpz2(xq_sp+Nb, subpz2(s08pjs4c_mv_s2apjs6e, h3_s19pjs5d_mv_s3bpjs7f));

                const Vec4d  j_a19p1a5d_m1_a3bp1a7f =  jxpz2(a19p1a5d_m1_a3bp1a7f);
                setpz2(xq_sp+N4, subpz2(a08p1a4c_m1_a2ap1a6e,  j_a19p1a5d_m1_a3bp1a7f));
                setpz2(xq_sp+Nc, addpz2(a08p1a4c_m1_a2ap1a6e,  j_a19p1a5d_m1_a3bp1a7f));

                const Vec4d hd_s19mjs5d_mw_s3bmjs7f = hdxpz2(s19mjs5d_mw_s3bmjs7f);
                setpz2(xq_sp+N5, subpz2(s08mjs4c_mw_s2amjs6e, hd_s19mjs5d_mw_s3bmjs7f));
                setpz2(xq_sp+Nd, addpz2(s08mjs4c_mw_s2amjs6e, hd_s19mjs5d_mw_s3bmjs7f));

                const Vec4d v8_a19m1a5d_pj_a3bm1a7f = v8xpz2(a19m1a5d_pj_a3bm1a7f);
                setpz2(xq_sp+N6, subpz2(a08m1a4c_pj_a2am1a6e, v8_a19m1a5d_pj_a3bm1a7f));
                setpz2(xq_sp+Ne, addpz2(a08m1a4c_pj_a2am1a6e, v8_a19m1a5d_pj_a3bm1a7f));

                const Vec4d hf_s19pjs5d_pv_s3bpjs7f = hfxpz2(s19pjs5d_pv_s3bpjs7f);
                setpz2(xq_sp+N7, subpz2(s08pjs4c_pv_s2apjs6e, hf_s19pjs5d_pv_s3bpjs7f));
                setpz2(xq_sp+Nf, addpz2(s08pjs4c_pv_s2apjs6e, hf_s19pjs5d_pv_s3bpjs7f));

            }
        }
    }
};

template <int N> struct fwdcore<N,1>
{
    static constexpr int N0 = 0;
    static constexpr int N1 = N/16;
    static constexpr int N2 = N1*2;
    static constexpr int N3 = N1*3;
    static constexpr int N4 = N1*4;
    static constexpr int N5 = N1*5;
    static constexpr int N6 = N1*6;
    static constexpr int N7 = N1*7;
    static constexpr int N8 = N1*8;
    static constexpr int N9 = N1*9;
    static constexpr int Na = N1*10;
    static constexpr int Nb = N1*11;
    static constexpr int Nc = N1*12;
    static constexpr int Nd = N1*13;
    static constexpr int Ne = N1*14;
    static constexpr int Nf = N1*15;

    void operator()(
            complex_vector x, complex_vector y, const_complex_vector W) const noexcept
    {
        for (int p = 0; p < N1; p += 2) {
            complex_vector x_p   = x + p;
            complex_vector y_16p = y + 16*p;

            const Vec4d w1p = getpz2(twid<16,N, 1>(W,p));
            const Vec4d w2p = getpz2(twid<16,N, 2>(W,p));
            const Vec4d w3p = getpz2(twid<16,N, 3>(W,p));
            const Vec4d w4p = getpz2(twid<16,N, 4>(W,p));
            const Vec4d w5p = getpz2(twid<16,N, 5>(W,p));
            const Vec4d w6p = getpz2(twid<16,N, 6>(W,p));
            const Vec4d w7p = getpz2(twid<16,N, 7>(W,p));
            const Vec4d w8p = getpz2(twid<16,N, 8>(W,p));
            const Vec4d w9p = getpz2(twid<16,N, 9>(W,p));
            const Vec4d wap = getpz2(twid<16,N,10>(W,p));
            const Vec4d wbp = getpz2(twid<16,N,11>(W,p));
            const Vec4d wcp = getpz2(twid<16,N,12>(W,p));
            const Vec4d wdp = getpz2(twid<16,N,13>(W,p));
            const Vec4d wep = getpz2(twid<16,N,14>(W,p));
            const Vec4d wfp = getpz2(twid<16,N,15>(W,p));

            const Vec4d ab  = getpz2(y_16p+0x00);
            const Vec4d cd  = getpz2(y_16p+0x02);
            const Vec4d ef  = getpz2(y_16p+0x04);
            const Vec4d gh  = getpz2(y_16p+0x06);
            const Vec4d ij  = getpz2(y_16p+0x08);
            const Vec4d kl  = getpz2(y_16p+0x0a);
            const Vec4d mn  = getpz2(y_16p+0x0c);
            const Vec4d op  = getpz2(y_16p+0x0e);
            const Vec4d AB  = getpz2(y_16p+0x10);
            const Vec4d CD  = getpz2(y_16p+0x12);
            const Vec4d EF  = getpz2(y_16p+0x14);
            const Vec4d GH  = getpz2(y_16p+0x16);
            const Vec4d IJ  = getpz2(y_16p+0x18);
            const Vec4d KL  = getpz2(y_16p+0x1a);
            const Vec4d MN  = getpz2(y_16p+0x1c);
            const Vec4d OP  = getpz2(y_16p+0x1e);

            const Vec4d y0 =             catlo(ab, AB);
            const Vec4d y1 = mulpz2(w1p, cathi(ab, AB));
            const Vec4d y2 = mulpz2(w2p, catlo(cd, CD));
            const Vec4d y3 = mulpz2(w3p, cathi(cd, CD));
            const Vec4d y4 = mulpz2(w4p, catlo(ef, EF));
            const Vec4d y5 = mulpz2(w5p, cathi(ef, EF));
            const Vec4d y6 = mulpz2(w6p, catlo(gh, GH));
            const Vec4d y7 = mulpz2(w7p, cathi(gh, GH));

            const Vec4d y8 = mulpz2(w8p, catlo(ij, IJ));
            const Vec4d y9 = mulpz2(w9p, cathi(ij, IJ));
            const Vec4d ya = mulpz2(wap, catlo(kl, KL));
            const Vec4d yb = mulpz2(wbp, cathi(kl, KL));
            const Vec4d yc = mulpz2(wcp, catlo(mn, MN));
            const Vec4d yd = mulpz2(wdp, cathi(mn, MN));
            const Vec4d ye = mulpz2(wep, catlo(op, OP));
            const Vec4d yf = mulpz2(wfp, cathi(op, OP));

            const Vec4d a08 = addpz2(y0, y8); const Vec4d s08 = subpz2(y0, y8);
            const Vec4d a4c = addpz2(y4, yc); const Vec4d s4c = subpz2(y4, yc);
            const Vec4d a2a = addpz2(y2, ya); const Vec4d s2a = subpz2(y2, ya);
            const Vec4d a6e = addpz2(y6, ye); const Vec4d s6e = subpz2(y6, ye);
            const Vec4d a19 = addpz2(y1, y9); const Vec4d s19 = subpz2(y1, y9);
            const Vec4d a5d = addpz2(y5, yd); const Vec4d s5d = subpz2(y5, yd);
            const Vec4d a3b = addpz2(y3, yb); const Vec4d s3b = subpz2(y3, yb);
            const Vec4d a7f = addpz2(y7, yf); const Vec4d s7f = subpz2(y7, yf);

            const Vec4d js4c = jxpz2(s4c);
            const Vec4d js6e = jxpz2(s6e);
            const Vec4d js5d = jxpz2(s5d);
            const Vec4d js7f = jxpz2(s7f);

            const Vec4d a08p1a4c = addpz2(a08, a4c); const Vec4d s08mjs4c = subpz2(s08, js4c);
            const Vec4d a08m1a4c = subpz2(a08, a4c); const Vec4d s08pjs4c = addpz2(s08, js4c);
            const Vec4d a2ap1a6e = addpz2(a2a, a6e); const Vec4d s2amjs6e = subpz2(s2a, js6e);
            const Vec4d a2am1a6e = subpz2(a2a, a6e); const Vec4d s2apjs6e = addpz2(s2a, js6e);
            const Vec4d a19p1a5d = addpz2(a19, a5d); const Vec4d s19mjs5d = subpz2(s19, js5d);
            const Vec4d a19m1a5d = subpz2(a19, a5d); const Vec4d s19pjs5d = addpz2(s19, js5d);
            const Vec4d a3bp1a7f = addpz2(a3b, a7f); const Vec4d s3bmjs7f = subpz2(s3b, js7f);
            const Vec4d a3bm1a7f = subpz2(a3b, a7f); const Vec4d s3bpjs7f = addpz2(s3b, js7f);

            const Vec4d w8_s2amjs6e = w8xpz2(s2amjs6e);
            const Vec4d  j_a2am1a6e =  jxpz2(a2am1a6e);
            const Vec4d v8_s2apjs6e = v8xpz2(s2apjs6e);

            const Vec4d a08p1a4c_p1_a2ap1a6e = addpz2(a08p1a4c,    a2ap1a6e);
            const Vec4d s08mjs4c_pw_s2amjs6e = addpz2(s08mjs4c, w8_s2amjs6e);
            const Vec4d a08m1a4c_mj_a2am1a6e = subpz2(a08m1a4c,  j_a2am1a6e);
            const Vec4d s08pjs4c_mv_s2apjs6e = subpz2(s08pjs4c, v8_s2apjs6e);
            const Vec4d a08p1a4c_m1_a2ap1a6e = subpz2(a08p1a4c,    a2ap1a6e);
            const Vec4d s08mjs4c_mw_s2amjs6e = subpz2(s08mjs4c, w8_s2amjs6e);
            const Vec4d a08m1a4c_pj_a2am1a6e = addpz2(a08m1a4c,  j_a2am1a6e);
            const Vec4d s08pjs4c_pv_s2apjs6e = addpz2(s08pjs4c, v8_s2apjs6e);

            const Vec4d w8_s3bmjs7f = w8xpz2(s3bmjs7f);
            const Vec4d  j_a3bm1a7f =  jxpz2(a3bm1a7f);
            const Vec4d v8_s3bpjs7f = v8xpz2(s3bpjs7f);

            const Vec4d a19p1a5d_p1_a3bp1a7f = addpz2(a19p1a5d,    a3bp1a7f);
            const Vec4d s19mjs5d_pw_s3bmjs7f = addpz2(s19mjs5d, w8_s3bmjs7f);
            const Vec4d a19m1a5d_mj_a3bm1a7f = subpz2(a19m1a5d,  j_a3bm1a7f);
            const Vec4d s19pjs5d_mv_s3bpjs7f = subpz2(s19pjs5d, v8_s3bpjs7f);
            const Vec4d a19p1a5d_m1_a3bp1a7f = subpz2(a19p1a5d,    a3bp1a7f);
            const Vec4d s19mjs5d_mw_s3bmjs7f = subpz2(s19mjs5d, w8_s3bmjs7f);
            const Vec4d a19m1a5d_pj_a3bm1a7f = addpz2(a19m1a5d,  j_a3bm1a7f);
            const Vec4d s19pjs5d_pv_s3bpjs7f = addpz2(s19pjs5d, v8_s3bpjs7f);

            setpz2(x_p+N0, addpz2(a08p1a4c_p1_a2ap1a6e,    a19p1a5d_p1_a3bp1a7f));
            setpz2(x_p+N8, subpz2(a08p1a4c_p1_a2ap1a6e,    a19p1a5d_p1_a3bp1a7f));

            const Vec4d h1_s19mjs5d_pw_s3bmjs7f = h1xpz2(s19mjs5d_pw_s3bmjs7f);
            setpz2(x_p+N1, addpz2(s08mjs4c_pw_s2amjs6e, h1_s19mjs5d_pw_s3bmjs7f));
            setpz2(x_p+N9, subpz2(s08mjs4c_pw_s2amjs6e, h1_s19mjs5d_pw_s3bmjs7f));

            const Vec4d w8_a19m1a5d_mj_a3bm1a7f = w8xpz2(a19m1a5d_mj_a3bm1a7f);
            setpz2(x_p+N2, addpz2(a08m1a4c_mj_a2am1a6e, w8_a19m1a5d_mj_a3bm1a7f));
            setpz2(x_p+Na, subpz2(a08m1a4c_mj_a2am1a6e, w8_a19m1a5d_mj_a3bm1a7f));

            const Vec4d h3_s19pjs5d_mv_s3bpjs7f = h3xpz2(s19pjs5d_mv_s3bpjs7f);
            setpz2(x_p+N3, addpz2(s08pjs4c_mv_s2apjs6e, h3_s19pjs5d_mv_s3bpjs7f));
            setpz2(x_p+Nb, subpz2(s08pjs4c_mv_s2apjs6e, h3_s19pjs5d_mv_s3bpjs7f));

            const Vec4d  j_a19p1a5d_m1_a3bp1a7f =  jxpz2(a19p1a5d_m1_a3bp1a7f);
            setpz2(x_p+N4, subpz2(a08p1a4c_m1_a2ap1a6e,  j_a19p1a5d_m1_a3bp1a7f));
            setpz2(x_p+Nc, addpz2(a08p1a4c_m1_a2ap1a6e,  j_a19p1a5d_m1_a3bp1a7f));

            const Vec4d hd_s19mjs5d_mw_s3bmjs7f = hdxpz2(s19mjs5d_mw_s3bmjs7f);
            setpz2(x_p+N5, subpz2(s08mjs4c_mw_s2amjs6e, hd_s19mjs5d_mw_s3bmjs7f));
            setpz2(x_p+Nd, addpz2(s08mjs4c_mw_s2amjs6e, hd_s19mjs5d_mw_s3bmjs7f));

            const Vec4d v8_a19m1a5d_pj_a3bm1a7f = v8xpz2(a19m1a5d_pj_a3bm1a7f);
            setpz2(x_p+N6, subpz2(a08m1a4c_pj_a2am1a6e, v8_a19m1a5d_pj_a3bm1a7f));
            setpz2(x_p+Ne, addpz2(a08m1a4c_pj_a2am1a6e, v8_a19m1a5d_pj_a3bm1a7f));

            const Vec4d hf_s19pjs5d_pv_s3bpjs7f = hfxpz2(s19pjs5d_pv_s3bpjs7f);
            setpz2(x_p+N7, subpz2(s08pjs4c_pv_s2apjs6e, hf_s19pjs5d_pv_s3bpjs7f));
            setpz2(x_p+Nf, addpz2(s08pjs4c_pv_s2apjs6e, hf_s19pjs5d_pv_s3bpjs7f));

        }
    }
};

///////////////////////////////////////////////////////////////////////////////

template <int n, int s, bool eo, int mode> struct fwdend;

//-----------------------------------------------------------------------------

template <int s, bool eo, int mode> struct fwdend<16,s,eo,mode>
{
    static constexpr int N = 16*s;

    void operator()(complex_vector x, complex_vector y) const noexcept
    {
        complex_vector z = eo ? y : x;
        for (int q = 0; q < s; q += 2) {
            complex_vector xq = x + q;
            complex_vector zq = z + q;

            const Vec4d z0 = scalepz2<N,mode>(getpz2(zq+s*0x0));
            const Vec4d z1 = scalepz2<N,mode>(getpz2(zq+s*0x1));
            const Vec4d z2 = scalepz2<N,mode>(getpz2(zq+s*0x2));
            const Vec4d z3 = scalepz2<N,mode>(getpz2(zq+s*0x3));
            const Vec4d z4 = scalepz2<N,mode>(getpz2(zq+s*0x4));
            const Vec4d z5 = scalepz2<N,mode>(getpz2(zq+s*0x5));
            const Vec4d z6 = scalepz2<N,mode>(getpz2(zq+s*0x6));
            const Vec4d z7 = scalepz2<N,mode>(getpz2(zq+s*0x7));
            const Vec4d z8 = scalepz2<N,mode>(getpz2(zq+s*0x8));
            const Vec4d z9 = scalepz2<N,mode>(getpz2(zq+s*0x9));
            const Vec4d za = scalepz2<N,mode>(getpz2(zq+s*0xa));
            const Vec4d zb = scalepz2<N,mode>(getpz2(zq+s*0xb));
            const Vec4d zc = scalepz2<N,mode>(getpz2(zq+s*0xc));
            const Vec4d zd = scalepz2<N,mode>(getpz2(zq+s*0xd));
            const Vec4d ze = scalepz2<N,mode>(getpz2(zq+s*0xe));
            const Vec4d zf = scalepz2<N,mode>(getpz2(zq+s*0xf));

            const Vec4d a08 = addpz2(z0, z8); const Vec4d s08 = subpz2(z0, z8);
            const Vec4d a4c = addpz2(z4, zc); const Vec4d s4c = subpz2(z4, zc);
            const Vec4d a2a = addpz2(z2, za); const Vec4d s2a = subpz2(z2, za);
            const Vec4d a6e = addpz2(z6, ze); const Vec4d s6e = subpz2(z6, ze);
            const Vec4d a19 = addpz2(z1, z9); const Vec4d s19 = subpz2(z1, z9);
            const Vec4d a5d = addpz2(z5, zd); const Vec4d s5d = subpz2(z5, zd);
            const Vec4d a3b = addpz2(z3, zb); const Vec4d s3b = subpz2(z3, zb);
            const Vec4d a7f = addpz2(z7, zf); const Vec4d s7f = subpz2(z7, zf);

            const Vec4d js4c = jxpz2(s4c);
            const Vec4d js6e = jxpz2(s6e);
            const Vec4d js5d = jxpz2(s5d);
            const Vec4d js7f = jxpz2(s7f);

            const Vec4d a08p1a4c = addpz2(a08, a4c); const Vec4d s08mjs4c = subpz2(s08, js4c);
            const Vec4d a08m1a4c = subpz2(a08, a4c); const Vec4d s08pjs4c = addpz2(s08, js4c);
            const Vec4d a2ap1a6e = addpz2(a2a, a6e); const Vec4d s2amjs6e = subpz2(s2a, js6e);
            const Vec4d a2am1a6e = subpz2(a2a, a6e); const Vec4d s2apjs6e = addpz2(s2a, js6e);
            const Vec4d a19p1a5d = addpz2(a19, a5d); const Vec4d s19mjs5d = subpz2(s19, js5d);
            const Vec4d a19m1a5d = subpz2(a19, a5d); const Vec4d s19pjs5d = addpz2(s19, js5d);
            const Vec4d a3bp1a7f = addpz2(a3b, a7f); const Vec4d s3bmjs7f = subpz2(s3b, js7f);
            const Vec4d a3bm1a7f = subpz2(a3b, a7f); const Vec4d s3bpjs7f = addpz2(s3b, js7f);

            const Vec4d w8_s2amjs6e = w8xpz2(s2amjs6e);
            const Vec4d  j_a2am1a6e =  jxpz2(a2am1a6e);
            const Vec4d v8_s2apjs6e = v8xpz2(s2apjs6e);

            const Vec4d a08p1a4c_p1_a2ap1a6e = addpz2(a08p1a4c,    a2ap1a6e);
            const Vec4d s08mjs4c_pw_s2amjs6e = addpz2(s08mjs4c, w8_s2amjs6e);
            const Vec4d a08m1a4c_mj_a2am1a6e = subpz2(a08m1a4c,  j_a2am1a6e);
            const Vec4d s08pjs4c_mv_s2apjs6e = subpz2(s08pjs4c, v8_s2apjs6e);
            const Vec4d a08p1a4c_m1_a2ap1a6e = subpz2(a08p1a4c,    a2ap1a6e);
            const Vec4d s08mjs4c_mw_s2amjs6e = subpz2(s08mjs4c, w8_s2amjs6e);
            const Vec4d a08m1a4c_pj_a2am1a6e = addpz2(a08m1a4c,  j_a2am1a6e);
            const Vec4d s08pjs4c_pv_s2apjs6e = addpz2(s08pjs4c, v8_s2apjs6e);

            const Vec4d w8_s3bmjs7f = w8xpz2(s3bmjs7f);
            const Vec4d  j_a3bm1a7f =  jxpz2(a3bm1a7f);
            const Vec4d v8_s3bpjs7f = v8xpz2(s3bpjs7f);

            const Vec4d a19p1a5d_p1_a3bp1a7f = addpz2(a19p1a5d,    a3bp1a7f);
            const Vec4d s19mjs5d_pw_s3bmjs7f = addpz2(s19mjs5d, w8_s3bmjs7f);
            const Vec4d a19m1a5d_mj_a3bm1a7f = subpz2(a19m1a5d,  j_a3bm1a7f);
            const Vec4d s19pjs5d_mv_s3bpjs7f = subpz2(s19pjs5d, v8_s3bpjs7f);
            const Vec4d a19p1a5d_m1_a3bp1a7f = subpz2(a19p1a5d,    a3bp1a7f);
            const Vec4d s19mjs5d_mw_s3bmjs7f = subpz2(s19mjs5d, w8_s3bmjs7f);
            const Vec4d a19m1a5d_pj_a3bm1a7f = addpz2(a19m1a5d,  j_a3bm1a7f);
            const Vec4d s19pjs5d_pv_s3bpjs7f = addpz2(s19pjs5d, v8_s3bpjs7f);

            const Vec4d h1_s19mjs5d_pw_s3bmjs7f = h1xpz2(s19mjs5d_pw_s3bmjs7f);
            const Vec4d w8_a19m1a5d_mj_a3bm1a7f = w8xpz2(a19m1a5d_mj_a3bm1a7f);
            const Vec4d h3_s19pjs5d_mv_s3bpjs7f = h3xpz2(s19pjs5d_mv_s3bpjs7f);
            const Vec4d  j_a19p1a5d_m1_a3bp1a7f =  jxpz2(a19p1a5d_m1_a3bp1a7f);
            const Vec4d hd_s19mjs5d_mw_s3bmjs7f = hdxpz2(s19mjs5d_mw_s3bmjs7f);
            const Vec4d v8_a19m1a5d_pj_a3bm1a7f = v8xpz2(a19m1a5d_pj_a3bm1a7f);
            const Vec4d hf_s19pjs5d_pv_s3bpjs7f = hfxpz2(s19pjs5d_pv_s3bpjs7f);

            setpz2(xq+s*0x0, addpz2(a08p1a4c_p1_a2ap1a6e,    a19p1a5d_p1_a3bp1a7f));
            setpz2(xq+s*0x1, addpz2(s08mjs4c_pw_s2amjs6e, h1_s19mjs5d_pw_s3bmjs7f));
            setpz2(xq+s*0x2, addpz2(a08m1a4c_mj_a2am1a6e, w8_a19m1a5d_mj_a3bm1a7f));
            setpz2(xq+s*0x3, addpz2(s08pjs4c_mv_s2apjs6e, h3_s19pjs5d_mv_s3bpjs7f));
            setpz2(xq+s*0x4, subpz2(a08p1a4c_m1_a2ap1a6e,  j_a19p1a5d_m1_a3bp1a7f));
            setpz2(xq+s*0x5, subpz2(s08mjs4c_mw_s2amjs6e, hd_s19mjs5d_mw_s3bmjs7f));
            setpz2(xq+s*0x6, subpz2(a08m1a4c_pj_a2am1a6e, v8_a19m1a5d_pj_a3bm1a7f));
            setpz2(xq+s*0x7, subpz2(s08pjs4c_pv_s2apjs6e, hf_s19pjs5d_pv_s3bpjs7f));

            setpz2(xq+s*0x8, subpz2(a08p1a4c_p1_a2ap1a6e,    a19p1a5d_p1_a3bp1a7f));
            setpz2(xq+s*0x9, subpz2(s08mjs4c_pw_s2amjs6e, h1_s19mjs5d_pw_s3bmjs7f));
            setpz2(xq+s*0xa, subpz2(a08m1a4c_mj_a2am1a6e, w8_a19m1a5d_mj_a3bm1a7f));
            setpz2(xq+s*0xb, subpz2(s08pjs4c_mv_s2apjs6e, h3_s19pjs5d_mv_s3bpjs7f));
            setpz2(xq+s*0xc, addpz2(a08p1a4c_m1_a2ap1a6e,  j_a19p1a5d_m1_a3bp1a7f));
            setpz2(xq+s*0xd, addpz2(s08mjs4c_mw_s2amjs6e, hd_s19mjs5d_mw_s3bmjs7f));
            setpz2(xq+s*0xe, addpz2(a08m1a4c_pj_a2am1a6e, v8_a19m1a5d_pj_a3bm1a7f));
            setpz2(xq+s*0xf, addpz2(s08pjs4c_pv_s2apjs6e, hf_s19pjs5d_pv_s3bpjs7f));
        }
    }
};

template <bool eo, int mode> struct fwdend<16,1,eo,mode>
{
    inline void operator()(complex_vector x, complex_vector y) const noexcept
    {
        complex_vector z = eo ? y : x;
        const Vec2d z0 = scalepz<16,mode>(getpz(z[0x0]));
        const Vec2d z1 = scalepz<16,mode>(getpz(z[0x1]));
        const Vec2d z2 = scalepz<16,mode>(getpz(z[0x2]));
        const Vec2d z3 = scalepz<16,mode>(getpz(z[0x3]));
        const Vec2d z4 = scalepz<16,mode>(getpz(z[0x4]));
        const Vec2d z5 = scalepz<16,mode>(getpz(z[0x5]));
        const Vec2d z6 = scalepz<16,mode>(getpz(z[0x6]));
        const Vec2d z7 = scalepz<16,mode>(getpz(z[0x7]));
        const Vec2d z8 = scalepz<16,mode>(getpz(z[0x8]));
        const Vec2d z9 = scalepz<16,mode>(getpz(z[0x9]));
        const Vec2d za = scalepz<16,mode>(getpz(z[0xa]));
        const Vec2d zb = scalepz<16,mode>(getpz(z[0xb]));
        const Vec2d zc = scalepz<16,mode>(getpz(z[0xc]));
        const Vec2d zd = scalepz<16,mode>(getpz(z[0xd]));
        const Vec2d ze = scalepz<16,mode>(getpz(z[0xe]));
        const Vec2d zf = scalepz<16,mode>(getpz(z[0xf]));

        const Vec2d a08 = addpz(z0, z8); const Vec2d s08 = subpz(z0, z8);
        const Vec2d a4c = addpz(z4, zc); const Vec2d s4c = subpz(z4, zc);
        const Vec2d a2a = addpz(z2, za); const Vec2d s2a = subpz(z2, za);
        const Vec2d a6e = addpz(z6, ze); const Vec2d s6e = subpz(z6, ze);
        const Vec2d a19 = addpz(z1, z9); const Vec2d s19 = subpz(z1, z9);
        const Vec2d a5d = addpz(z5, zd); const Vec2d s5d = subpz(z5, zd);
        const Vec2d a3b = addpz(z3, zb); const Vec2d s3b = subpz(z3, zb);
        const Vec2d a7f = addpz(z7, zf); const Vec2d s7f = subpz(z7, zf);

        const Vec2d js4c = jxpz(s4c);
        const Vec2d js6e = jxpz(s6e);
        const Vec2d js5d = jxpz(s5d);
        const Vec2d js7f = jxpz(s7f);

        const Vec2d a08p1a4c = addpz(a08, a4c); const Vec2d s08mjs4c = subpz(s08, js4c);
        const Vec2d a08m1a4c = subpz(a08, a4c); const Vec2d s08pjs4c = addpz(s08, js4c);
        const Vec2d a2ap1a6e = addpz(a2a, a6e); const Vec2d s2amjs6e = subpz(s2a, js6e);
        const Vec2d a2am1a6e = subpz(a2a, a6e); const Vec2d s2apjs6e = addpz(s2a, js6e);
        const Vec2d a19p1a5d = addpz(a19, a5d); const Vec2d s19mjs5d = subpz(s19, js5d);
        const Vec2d a19m1a5d = subpz(a19, a5d); const Vec2d s19pjs5d = addpz(s19, js5d);
        const Vec2d a3bp1a7f = addpz(a3b, a7f); const Vec2d s3bmjs7f = subpz(s3b, js7f);
        const Vec2d a3bm1a7f = subpz(a3b, a7f); const Vec2d s3bpjs7f = addpz(s3b, js7f);

        const Vec2d w8_s2amjs6e = w8xpz(s2amjs6e);
        const Vec2d  j_a2am1a6e =  jxpz(a2am1a6e);
        const Vec2d v8_s2apjs6e = v8xpz(s2apjs6e);

        const Vec2d a08p1a4c_p1_a2ap1a6e = addpz(a08p1a4c,    a2ap1a6e);
        const Vec2d s08mjs4c_pw_s2amjs6e = addpz(s08mjs4c, w8_s2amjs6e);
        const Vec2d a08m1a4c_mj_a2am1a6e = subpz(a08m1a4c,  j_a2am1a6e);
        const Vec2d s08pjs4c_mv_s2apjs6e = subpz(s08pjs4c, v8_s2apjs6e);
        const Vec2d a08p1a4c_m1_a2ap1a6e = subpz(a08p1a4c,    a2ap1a6e);
        const Vec2d s08mjs4c_mw_s2amjs6e = subpz(s08mjs4c, w8_s2amjs6e);
        const Vec2d a08m1a4c_pj_a2am1a6e = addpz(a08m1a4c,  j_a2am1a6e);
        const Vec2d s08pjs4c_pv_s2apjs6e = addpz(s08pjs4c, v8_s2apjs6e);

        const Vec2d w8_s3bmjs7f = w8xpz(s3bmjs7f);
        const Vec2d  j_a3bm1a7f =  jxpz(a3bm1a7f);
        const Vec2d v8_s3bpjs7f = v8xpz(s3bpjs7f);

        const Vec2d a19p1a5d_p1_a3bp1a7f = addpz(a19p1a5d,    a3bp1a7f);
        const Vec2d s19mjs5d_pw_s3bmjs7f = addpz(s19mjs5d, w8_s3bmjs7f);
        const Vec2d a19m1a5d_mj_a3bm1a7f = subpz(a19m1a5d,  j_a3bm1a7f);
        const Vec2d s19pjs5d_mv_s3bpjs7f = subpz(s19pjs5d, v8_s3bpjs7f);
        const Vec2d a19p1a5d_m1_a3bp1a7f = subpz(a19p1a5d,    a3bp1a7f);
        const Vec2d s19mjs5d_mw_s3bmjs7f = subpz(s19mjs5d, w8_s3bmjs7f);
        const Vec2d a19m1a5d_pj_a3bm1a7f = addpz(a19m1a5d,  j_a3bm1a7f);
        const Vec2d s19pjs5d_pv_s3bpjs7f = addpz(s19pjs5d, v8_s3bpjs7f);

        const Vec2d h1_s19mjs5d_pw_s3bmjs7f = h1xpz(s19mjs5d_pw_s3bmjs7f);
        const Vec2d w8_a19m1a5d_mj_a3bm1a7f = w8xpz(a19m1a5d_mj_a3bm1a7f);
        const Vec2d h3_s19pjs5d_mv_s3bpjs7f = h3xpz(s19pjs5d_mv_s3bpjs7f);
        const Vec2d  j_a19p1a5d_m1_a3bp1a7f =  jxpz(a19p1a5d_m1_a3bp1a7f);
        const Vec2d hd_s19mjs5d_mw_s3bmjs7f = hdxpz(s19mjs5d_mw_s3bmjs7f);
        const Vec2d v8_a19m1a5d_pj_a3bm1a7f = v8xpz(a19m1a5d_pj_a3bm1a7f);
        const Vec2d hf_s19pjs5d_pv_s3bpjs7f = hfxpz(s19pjs5d_pv_s3bpjs7f);

        setpz(x[0x0], addpz(a08p1a4c_p1_a2ap1a6e,    a19p1a5d_p1_a3bp1a7f));
        setpz(x[0x1], addpz(s08mjs4c_pw_s2amjs6e, h1_s19mjs5d_pw_s3bmjs7f));
        setpz(x[0x2], addpz(a08m1a4c_mj_a2am1a6e, w8_a19m1a5d_mj_a3bm1a7f));
        setpz(x[0x3], addpz(s08pjs4c_mv_s2apjs6e, h3_s19pjs5d_mv_s3bpjs7f));
        setpz(x[0x4], subpz(a08p1a4c_m1_a2ap1a6e,  j_a19p1a5d_m1_a3bp1a7f));
        setpz(x[0x5], subpz(s08mjs4c_mw_s2amjs6e, hd_s19mjs5d_mw_s3bmjs7f));
        setpz(x[0x6], subpz(a08m1a4c_pj_a2am1a6e, v8_a19m1a5d_pj_a3bm1a7f));
        setpz(x[0x7], subpz(s08pjs4c_pv_s2apjs6e, hf_s19pjs5d_pv_s3bpjs7f));

        setpz(x[0x8], subpz(a08p1a4c_p1_a2ap1a6e,    a19p1a5d_p1_a3bp1a7f));
        setpz(x[0x9], subpz(s08mjs4c_pw_s2amjs6e, h1_s19mjs5d_pw_s3bmjs7f));
        setpz(x[0xa], subpz(a08m1a4c_mj_a2am1a6e, w8_a19m1a5d_mj_a3bm1a7f));
        setpz(x[0xb], subpz(s08pjs4c_mv_s2apjs6e, h3_s19pjs5d_mv_s3bpjs7f));
        setpz(x[0xc], addpz(a08p1a4c_m1_a2ap1a6e,  j_a19p1a5d_m1_a3bp1a7f));
        setpz(x[0xd], addpz(s08mjs4c_mw_s2amjs6e, hd_s19mjs5d_mw_s3bmjs7f));
        setpz(x[0xe], addpz(a08m1a4c_pj_a2am1a6e, v8_a19m1a5d_pj_a3bm1a7f));
        setpz(x[0xf], addpz(s08pjs4c_pv_s2apjs6e, hf_s19pjs5d_pv_s3bpjs7f));
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
        fwdfft<n/16,16*s,!eo,mode>()(y, x, W);
        fwdcore<n,s>()(x, y, W);
    }
};

template <int s, bool eo, int mode> struct fwdfft<16,s,eo,mode>
{
    inline void operator()(
        complex_vector x, complex_vector y, const_complex_vector) const noexcept
    {
        fwdend<16,s,eo,mode>()(x, y);
    }
};

template <int s, bool eo, int mode> struct fwdfft<8,s,eo,mode>
{
    inline void operator()(
        complex_vector x, complex_vector y, const_complex_vector) const noexcept
    {
        OTFFT_AVXDIT8::fwdend<8,s,eo,mode>()(x, y);
    }
};

template <int s, bool eo, int mode> struct fwdfft<4,s,eo,mode>
{
    inline void operator()(
        complex_vector x, complex_vector y, const_complex_vector) const noexcept
    {
        OTFFT_AVXDIT4::fwdend<4,s,eo,mode>()(x, y);
    }
};

template <int s, bool eo, int mode> struct fwdfft<2,s,eo,mode>
{
    inline void operator()(
        complex_vector x, complex_vector y, const_complex_vector) const noexcept
    {
        OTFFT_AVXDIT4::fwdend<2,s,eo,mode>()(x, y);
    }
};

///////////////////////////////////////////////////////////////////////////////
// Inverse Butterfly Operation
///////////////////////////////////////////////////////////////////////////////

template <int n, int s> struct invcore
{
    static constexpr int m  = n/16;
    static constexpr int N  = n*s;
    static constexpr int N0 = 0;
    static constexpr int N1 = N/16;
    static constexpr int N2 = N1*2;
    static constexpr int N3 = N1*3;
    static constexpr int N4 = N1*4;
    static constexpr int N5 = N1*5;
    static constexpr int N6 = N1*6;
    static constexpr int N7 = N1*7;
    static constexpr int N8 = N1*8;
    static constexpr int N9 = N1*9;
    static constexpr int Na = N1*10;
    static constexpr int Nb = N1*11;
    static constexpr int Nc = N1*12;
    static constexpr int Nd = N1*13;
    static constexpr int Ne = N1*14;
    static constexpr int Nf = N1*15;

    void operator()(
            complex_vector x, complex_vector y, const_complex_vector W) const noexcept
    {
        for (int p = 0; p < m; p++) {
            const int sp   = s*p;
            const int s16p = 16*sp;

            const Vec4d w1p = cnjpz2(duppz3(*twidT<16,N, 1>(W,sp)));
            const Vec4d w2p = cnjpz2(duppz3(*twidT<16,N, 2>(W,sp)));
            const Vec4d w3p = cnjpz2(duppz3(*twidT<16,N, 3>(W,sp)));
            const Vec4d w4p = cnjpz2(duppz3(*twidT<16,N, 4>(W,sp)));
            const Vec4d w5p = cnjpz2(duppz3(*twidT<16,N, 5>(W,sp)));
            const Vec4d w6p = cnjpz2(duppz3(*twidT<16,N, 6>(W,sp)));
            const Vec4d w7p = cnjpz2(duppz3(*twidT<16,N, 7>(W,sp)));
            const Vec4d w8p = cnjpz2(duppz3(*twidT<16,N, 8>(W,sp)));
            const Vec4d w9p = cnjpz2(duppz3(*twidT<16,N, 9>(W,sp)));
            const Vec4d wap = cnjpz2(duppz3(*twidT<16,N,10>(W,sp)));
            const Vec4d wbp = cnjpz2(duppz3(*twidT<16,N,11>(W,sp)));
            const Vec4d wcp = cnjpz2(duppz3(*twidT<16,N,12>(W,sp)));
            const Vec4d wdp = cnjpz2(duppz3(*twidT<16,N,13>(W,sp)));
            const Vec4d wep = cnjpz2(duppz3(*twidT<16,N,14>(W,sp)));
            const Vec4d wfp = cnjpz2(duppz3(*twidT<16,N,15>(W,sp)));

            for (int q = 0; q < s; q += 2) {
                complex_vector xq_sp   = x + q + sp;
                complex_vector yq_s16p = y + q + s16p;

                const Vec4d y0 =             getpz2(yq_s16p+s*0x0);
                const Vec4d y1 = mulpz2(w1p, getpz2(yq_s16p+s*0x1));
                const Vec4d y2 = mulpz2(w2p, getpz2(yq_s16p+s*0x2));
                const Vec4d y3 = mulpz2(w3p, getpz2(yq_s16p+s*0x3));
                const Vec4d y4 = mulpz2(w4p, getpz2(yq_s16p+s*0x4));
                const Vec4d y5 = mulpz2(w5p, getpz2(yq_s16p+s*0x5));
                const Vec4d y6 = mulpz2(w6p, getpz2(yq_s16p+s*0x6));
                const Vec4d y7 = mulpz2(w7p, getpz2(yq_s16p+s*0x7));
                const Vec4d y8 = mulpz2(w8p, getpz2(yq_s16p+s*0x8));
                const Vec4d y9 = mulpz2(w9p, getpz2(yq_s16p+s*0x9));
                const Vec4d ya = mulpz2(wap, getpz2(yq_s16p+s*0xa));
                const Vec4d yb = mulpz2(wbp, getpz2(yq_s16p+s*0xb));
                const Vec4d yc = mulpz2(wcp, getpz2(yq_s16p+s*0xc));
                const Vec4d yd = mulpz2(wdp, getpz2(yq_s16p+s*0xd));
                const Vec4d ye = mulpz2(wep, getpz2(yq_s16p+s*0xe));
                const Vec4d yf = mulpz2(wfp, getpz2(yq_s16p+s*0xf));

                const Vec4d a08 = addpz2(y0, y8); const Vec4d s08 = subpz2(y0, y8);
                const Vec4d a4c = addpz2(y4, yc); const Vec4d s4c = subpz2(y4, yc);
                const Vec4d a2a = addpz2(y2, ya); const Vec4d s2a = subpz2(y2, ya);
                const Vec4d a6e = addpz2(y6, ye); const Vec4d s6e = subpz2(y6, ye);
                const Vec4d a19 = addpz2(y1, y9); const Vec4d s19 = subpz2(y1, y9);
                const Vec4d a5d = addpz2(y5, yd); const Vec4d s5d = subpz2(y5, yd);
                const Vec4d a3b = addpz2(y3, yb); const Vec4d s3b = subpz2(y3, yb);
                const Vec4d a7f = addpz2(y7, yf); const Vec4d s7f = subpz2(y7, yf);

                const Vec4d js4c = jxpz2(s4c);
                const Vec4d js6e = jxpz2(s6e);
                const Vec4d js5d = jxpz2(s5d);
                const Vec4d js7f = jxpz2(s7f);

                const Vec4d a08p1a4c = addpz2(a08, a4c); const Vec4d s08mjs4c = subpz2(s08, js4c);
                const Vec4d a08m1a4c = subpz2(a08, a4c); const Vec4d s08pjs4c = addpz2(s08, js4c);
                const Vec4d a2ap1a6e = addpz2(a2a, a6e); const Vec4d s2amjs6e = subpz2(s2a, js6e);
                const Vec4d a2am1a6e = subpz2(a2a, a6e); const Vec4d s2apjs6e = addpz2(s2a, js6e);
                const Vec4d a19p1a5d = addpz2(a19, a5d); const Vec4d s19mjs5d = subpz2(s19, js5d);
                const Vec4d a19m1a5d = subpz2(a19, a5d); const Vec4d s19pjs5d = addpz2(s19, js5d);
                const Vec4d a3bp1a7f = addpz2(a3b, a7f); const Vec4d s3bmjs7f = subpz2(s3b, js7f);
                const Vec4d a3bm1a7f = subpz2(a3b, a7f); const Vec4d s3bpjs7f = addpz2(s3b, js7f);

                const Vec4d w8_s2amjs6e = w8xpz2(s2amjs6e);
                const Vec4d  j_a2am1a6e =  jxpz2(a2am1a6e);
                const Vec4d v8_s2apjs6e = v8xpz2(s2apjs6e);

                const Vec4d a08p1a4c_p1_a2ap1a6e = addpz2(a08p1a4c,    a2ap1a6e);
                const Vec4d s08mjs4c_pw_s2amjs6e = addpz2(s08mjs4c, w8_s2amjs6e);
                const Vec4d a08m1a4c_mj_a2am1a6e = subpz2(a08m1a4c,  j_a2am1a6e);
                const Vec4d s08pjs4c_mv_s2apjs6e = subpz2(s08pjs4c, v8_s2apjs6e);
                const Vec4d a08p1a4c_m1_a2ap1a6e = subpz2(a08p1a4c,    a2ap1a6e);
                const Vec4d s08mjs4c_mw_s2amjs6e = subpz2(s08mjs4c, w8_s2amjs6e);
                const Vec4d a08m1a4c_pj_a2am1a6e = addpz2(a08m1a4c,  j_a2am1a6e);
                const Vec4d s08pjs4c_pv_s2apjs6e = addpz2(s08pjs4c, v8_s2apjs6e);

                const Vec4d w8_s3bmjs7f = w8xpz2(s3bmjs7f);
                const Vec4d  j_a3bm1a7f =  jxpz2(a3bm1a7f);
                const Vec4d v8_s3bpjs7f = v8xpz2(s3bpjs7f);

                const Vec4d a19p1a5d_p1_a3bp1a7f = addpz2(a19p1a5d,    a3bp1a7f);
                const Vec4d s19mjs5d_pw_s3bmjs7f = addpz2(s19mjs5d, w8_s3bmjs7f);
                const Vec4d a19m1a5d_mj_a3bm1a7f = subpz2(a19m1a5d,  j_a3bm1a7f);
                const Vec4d s19pjs5d_mv_s3bpjs7f = subpz2(s19pjs5d, v8_s3bpjs7f);
                const Vec4d a19p1a5d_m1_a3bp1a7f = subpz2(a19p1a5d,    a3bp1a7f);
                const Vec4d s19mjs5d_mw_s3bmjs7f = subpz2(s19mjs5d, w8_s3bmjs7f);
                const Vec4d a19m1a5d_pj_a3bm1a7f = addpz2(a19m1a5d,  j_a3bm1a7f);
                const Vec4d s19pjs5d_pv_s3bpjs7f = addpz2(s19pjs5d, v8_s3bpjs7f);

                setpz2(xq_sp+N0, addpz2(a08p1a4c_p1_a2ap1a6e,    a19p1a5d_p1_a3bp1a7f));
                setpz2(xq_sp+N8, subpz2(a08p1a4c_p1_a2ap1a6e,    a19p1a5d_p1_a3bp1a7f));

                const Vec4d hf_s19pjs5d_pv_s3bpjs7f = hfxpz2(s19pjs5d_pv_s3bpjs7f);
                setpz2(xq_sp+N1, addpz2(s08pjs4c_pv_s2apjs6e, hf_s19pjs5d_pv_s3bpjs7f));
                setpz2(xq_sp+N9, subpz2(s08pjs4c_pv_s2apjs6e, hf_s19pjs5d_pv_s3bpjs7f));

                const Vec4d v8_a19m1a5d_pj_a3bm1a7f = v8xpz2(a19m1a5d_pj_a3bm1a7f);
                setpz2(xq_sp+N2, addpz2(a08m1a4c_pj_a2am1a6e, v8_a19m1a5d_pj_a3bm1a7f));
                setpz2(xq_sp+Na, subpz2(a08m1a4c_pj_a2am1a6e, v8_a19m1a5d_pj_a3bm1a7f));

                const Vec4d hd_s19mjs5d_mw_s3bmjs7f = hdxpz2(s19mjs5d_mw_s3bmjs7f);
                setpz2(xq_sp+N3, addpz2(s08mjs4c_mw_s2amjs6e, hd_s19mjs5d_mw_s3bmjs7f));
                setpz2(xq_sp+Nb, subpz2(s08mjs4c_mw_s2amjs6e, hd_s19mjs5d_mw_s3bmjs7f));

                const Vec4d  j_a19p1a5d_m1_a3bp1a7f =  jxpz2(a19p1a5d_m1_a3bp1a7f);
                setpz2(xq_sp+N4, addpz2(a08p1a4c_m1_a2ap1a6e,  j_a19p1a5d_m1_a3bp1a7f));
                setpz2(xq_sp+Nc, subpz2(a08p1a4c_m1_a2ap1a6e,  j_a19p1a5d_m1_a3bp1a7f));

                const Vec4d h3_s19pjs5d_mv_s3bpjs7f = h3xpz2(s19pjs5d_mv_s3bpjs7f);
                setpz2(xq_sp+N5, subpz2(s08pjs4c_mv_s2apjs6e, h3_s19pjs5d_mv_s3bpjs7f));
                setpz2(xq_sp+Nd, addpz2(s08pjs4c_mv_s2apjs6e, h3_s19pjs5d_mv_s3bpjs7f));

                const Vec4d w8_a19m1a5d_mj_a3bm1a7f = w8xpz2(a19m1a5d_mj_a3bm1a7f);
                setpz2(xq_sp+N6, subpz2(a08m1a4c_mj_a2am1a6e, w8_a19m1a5d_mj_a3bm1a7f));
                setpz2(xq_sp+Ne, addpz2(a08m1a4c_mj_a2am1a6e, w8_a19m1a5d_mj_a3bm1a7f));

                const Vec4d h1_s19mjs5d_pw_s3bmjs7f = h1xpz2(s19mjs5d_pw_s3bmjs7f);
                setpz2(xq_sp+N7, subpz2(s08mjs4c_pw_s2amjs6e, h1_s19mjs5d_pw_s3bmjs7f));
                setpz2(xq_sp+Nf, addpz2(s08mjs4c_pw_s2amjs6e, h1_s19mjs5d_pw_s3bmjs7f));

            }
        }
    }
};

template <int N> struct invcore<N,1>
{
    static constexpr int N0 = 0;
    static constexpr int N1 = N/16;
    static constexpr int N2 = N1*2;
    static constexpr int N3 = N1*3;
    static constexpr int N4 = N1*4;
    static constexpr int N5 = N1*5;
    static constexpr int N6 = N1*6;
    static constexpr int N7 = N1*7;
    static constexpr int N8 = N1*8;
    static constexpr int N9 = N1*9;
    static constexpr int Na = N1*10;
    static constexpr int Nb = N1*11;
    static constexpr int Nc = N1*12;
    static constexpr int Nd = N1*13;
    static constexpr int Ne = N1*14;
    static constexpr int Nf = N1*15;

    void operator()(
            complex_vector x, complex_vector y, const_complex_vector W) const noexcept
    {
        for (int p = 0; p < N1; p += 2) {
            complex_vector x_p   = x + p;
            complex_vector y_16p = y + 16*p;

            const Vec4d w1p = cnjpz2(getpz2(twid<16,N, 1>(W,p)));
            const Vec4d w2p = cnjpz2(getpz2(twid<16,N, 2>(W,p)));
            const Vec4d w3p = cnjpz2(getpz2(twid<16,N, 3>(W,p)));
            const Vec4d w4p = cnjpz2(getpz2(twid<16,N, 4>(W,p)));
            const Vec4d w5p = cnjpz2(getpz2(twid<16,N, 5>(W,p)));
            const Vec4d w6p = cnjpz2(getpz2(twid<16,N, 6>(W,p)));
            const Vec4d w7p = cnjpz2(getpz2(twid<16,N, 7>(W,p)));
            const Vec4d w8p = cnjpz2(getpz2(twid<16,N, 8>(W,p)));
            const Vec4d w9p = cnjpz2(getpz2(twid<16,N, 9>(W,p)));
            const Vec4d wap = cnjpz2(getpz2(twid<16,N,10>(W,p)));
            const Vec4d wbp = cnjpz2(getpz2(twid<16,N,11>(W,p)));
            const Vec4d wcp = cnjpz2(getpz2(twid<16,N,12>(W,p)));
            const Vec4d wdp = cnjpz2(getpz2(twid<16,N,13>(W,p)));
            const Vec4d wep = cnjpz2(getpz2(twid<16,N,14>(W,p)));
            const Vec4d wfp = cnjpz2(getpz2(twid<16,N,15>(W,p)));

            const Vec4d ab  = getpz2(y_16p+0x00);
            const Vec4d cd  = getpz2(y_16p+0x02);
            const Vec4d ef  = getpz2(y_16p+0x04);
            const Vec4d gh  = getpz2(y_16p+0x06);
            const Vec4d ij  = getpz2(y_16p+0x08);
            const Vec4d kl  = getpz2(y_16p+0x0a);
            const Vec4d mn  = getpz2(y_16p+0x0c);
            const Vec4d op  = getpz2(y_16p+0x0e);
            const Vec4d AB  = getpz2(y_16p+0x10);
            const Vec4d CD  = getpz2(y_16p+0x12);
            const Vec4d EF  = getpz2(y_16p+0x14);
            const Vec4d GH  = getpz2(y_16p+0x16);
            const Vec4d IJ  = getpz2(y_16p+0x18);
            const Vec4d KL  = getpz2(y_16p+0x1a);
            const Vec4d MN  = getpz2(y_16p+0x1c);
            const Vec4d OP  = getpz2(y_16p+0x1e);

            const Vec4d y0 =             catlo(ab, AB);
            const Vec4d y1 = mulpz2(w1p, cathi(ab, AB));
            const Vec4d y2 = mulpz2(w2p, catlo(cd, CD));
            const Vec4d y3 = mulpz2(w3p, cathi(cd, CD));
            const Vec4d y4 = mulpz2(w4p, catlo(ef, EF));
            const Vec4d y5 = mulpz2(w5p, cathi(ef, EF));
            const Vec4d y6 = mulpz2(w6p, catlo(gh, GH));
            const Vec4d y7 = mulpz2(w7p, cathi(gh, GH));

            const Vec4d y8 = mulpz2(w8p, catlo(ij, IJ));
            const Vec4d y9 = mulpz2(w9p, cathi(ij, IJ));
            const Vec4d ya = mulpz2(wap, catlo(kl, KL));
            const Vec4d yb = mulpz2(wbp, cathi(kl, KL));
            const Vec4d yc = mulpz2(wcp, catlo(mn, MN));
            const Vec4d yd = mulpz2(wdp, cathi(mn, MN));
            const Vec4d ye = mulpz2(wep, catlo(op, OP));
            const Vec4d yf = mulpz2(wfp, cathi(op, OP));

            const Vec4d a08 = addpz2(y0, y8); const Vec4d s08 = subpz2(y0, y8);
            const Vec4d a4c = addpz2(y4, yc); const Vec4d s4c = subpz2(y4, yc);
            const Vec4d a2a = addpz2(y2, ya); const Vec4d s2a = subpz2(y2, ya);
            const Vec4d a6e = addpz2(y6, ye); const Vec4d s6e = subpz2(y6, ye);
            const Vec4d a19 = addpz2(y1, y9); const Vec4d s19 = subpz2(y1, y9);
            const Vec4d a5d = addpz2(y5, yd); const Vec4d s5d = subpz2(y5, yd);
            const Vec4d a3b = addpz2(y3, yb); const Vec4d s3b = subpz2(y3, yb);
            const Vec4d a7f = addpz2(y7, yf); const Vec4d s7f = subpz2(y7, yf);

            const Vec4d js4c = jxpz2(s4c);
            const Vec4d js6e = jxpz2(s6e);
            const Vec4d js5d = jxpz2(s5d);
            const Vec4d js7f = jxpz2(s7f);

            const Vec4d a08p1a4c = addpz2(a08, a4c); const Vec4d s08mjs4c = subpz2(s08, js4c);
            const Vec4d a08m1a4c = subpz2(a08, a4c); const Vec4d s08pjs4c = addpz2(s08, js4c);
            const Vec4d a2ap1a6e = addpz2(a2a, a6e); const Vec4d s2amjs6e = subpz2(s2a, js6e);
            const Vec4d a2am1a6e = subpz2(a2a, a6e); const Vec4d s2apjs6e = addpz2(s2a, js6e);
            const Vec4d a19p1a5d = addpz2(a19, a5d); const Vec4d s19mjs5d = subpz2(s19, js5d);
            const Vec4d a19m1a5d = subpz2(a19, a5d); const Vec4d s19pjs5d = addpz2(s19, js5d);
            const Vec4d a3bp1a7f = addpz2(a3b, a7f); const Vec4d s3bmjs7f = subpz2(s3b, js7f);
            const Vec4d a3bm1a7f = subpz2(a3b, a7f); const Vec4d s3bpjs7f = addpz2(s3b, js7f);

            const Vec4d w8_s2amjs6e = w8xpz2(s2amjs6e);
            const Vec4d  j_a2am1a6e =  jxpz2(a2am1a6e);
            const Vec4d v8_s2apjs6e = v8xpz2(s2apjs6e);

            const Vec4d a08p1a4c_p1_a2ap1a6e = addpz2(a08p1a4c,    a2ap1a6e);
            const Vec4d s08mjs4c_pw_s2amjs6e = addpz2(s08mjs4c, w8_s2amjs6e);
            const Vec4d a08m1a4c_mj_a2am1a6e = subpz2(a08m1a4c,  j_a2am1a6e);
            const Vec4d s08pjs4c_mv_s2apjs6e = subpz2(s08pjs4c, v8_s2apjs6e);
            const Vec4d a08p1a4c_m1_a2ap1a6e = subpz2(a08p1a4c,    a2ap1a6e);
            const Vec4d s08mjs4c_mw_s2amjs6e = subpz2(s08mjs4c, w8_s2amjs6e);
            const Vec4d a08m1a4c_pj_a2am1a6e = addpz2(a08m1a4c,  j_a2am1a6e);
            const Vec4d s08pjs4c_pv_s2apjs6e = addpz2(s08pjs4c, v8_s2apjs6e);

            const Vec4d w8_s3bmjs7f = w8xpz2(s3bmjs7f);
            const Vec4d  j_a3bm1a7f =  jxpz2(a3bm1a7f);
            const Vec4d v8_s3bpjs7f = v8xpz2(s3bpjs7f);

            const Vec4d a19p1a5d_p1_a3bp1a7f = addpz2(a19p1a5d,    a3bp1a7f);
            const Vec4d s19mjs5d_pw_s3bmjs7f = addpz2(s19mjs5d, w8_s3bmjs7f);
            const Vec4d a19m1a5d_mj_a3bm1a7f = subpz2(a19m1a5d,  j_a3bm1a7f);
            const Vec4d s19pjs5d_mv_s3bpjs7f = subpz2(s19pjs5d, v8_s3bpjs7f);
            const Vec4d a19p1a5d_m1_a3bp1a7f = subpz2(a19p1a5d,    a3bp1a7f);
            const Vec4d s19mjs5d_mw_s3bmjs7f = subpz2(s19mjs5d, w8_s3bmjs7f);
            const Vec4d a19m1a5d_pj_a3bm1a7f = addpz2(a19m1a5d,  j_a3bm1a7f);
            const Vec4d s19pjs5d_pv_s3bpjs7f = addpz2(s19pjs5d, v8_s3bpjs7f);

            setpz2(x_p+N0, addpz2(a08p1a4c_p1_a2ap1a6e,    a19p1a5d_p1_a3bp1a7f));
            setpz2(x_p+N8, subpz2(a08p1a4c_p1_a2ap1a6e,    a19p1a5d_p1_a3bp1a7f));

            const Vec4d hf_s19pjs5d_pv_s3bpjs7f = hfxpz2(s19pjs5d_pv_s3bpjs7f);
            setpz2(x_p+N1, addpz2(s08pjs4c_pv_s2apjs6e, hf_s19pjs5d_pv_s3bpjs7f));
            setpz2(x_p+N9, subpz2(s08pjs4c_pv_s2apjs6e, hf_s19pjs5d_pv_s3bpjs7f));

            const Vec4d v8_a19m1a5d_pj_a3bm1a7f = v8xpz2(a19m1a5d_pj_a3bm1a7f);
            setpz2(x_p+N2, addpz2(a08m1a4c_pj_a2am1a6e, v8_a19m1a5d_pj_a3bm1a7f));
            setpz2(x_p+Na, subpz2(a08m1a4c_pj_a2am1a6e, v8_a19m1a5d_pj_a3bm1a7f));

            const Vec4d hd_s19mjs5d_mw_s3bmjs7f = hdxpz2(s19mjs5d_mw_s3bmjs7f);
            setpz2(x_p+N3, addpz2(s08mjs4c_mw_s2amjs6e, hd_s19mjs5d_mw_s3bmjs7f));
            setpz2(x_p+Nb, subpz2(s08mjs4c_mw_s2amjs6e, hd_s19mjs5d_mw_s3bmjs7f));

            const Vec4d  j_a19p1a5d_m1_a3bp1a7f =  jxpz2(a19p1a5d_m1_a3bp1a7f);
            setpz2(x_p+N4, addpz2(a08p1a4c_m1_a2ap1a6e,  j_a19p1a5d_m1_a3bp1a7f));
            setpz2(x_p+Nc, subpz2(a08p1a4c_m1_a2ap1a6e,  j_a19p1a5d_m1_a3bp1a7f));

            const Vec4d h3_s19pjs5d_mv_s3bpjs7f = h3xpz2(s19pjs5d_mv_s3bpjs7f);
            setpz2(x_p+N5, subpz2(s08pjs4c_mv_s2apjs6e, h3_s19pjs5d_mv_s3bpjs7f));
            setpz2(x_p+Nd, addpz2(s08pjs4c_mv_s2apjs6e, h3_s19pjs5d_mv_s3bpjs7f));

            const Vec4d w8_a19m1a5d_mj_a3bm1a7f = w8xpz2(a19m1a5d_mj_a3bm1a7f);
            setpz2(x_p+N6, subpz2(a08m1a4c_mj_a2am1a6e, w8_a19m1a5d_mj_a3bm1a7f));
            setpz2(x_p+Ne, addpz2(a08m1a4c_mj_a2am1a6e, w8_a19m1a5d_mj_a3bm1a7f));

            const Vec4d h1_s19mjs5d_pw_s3bmjs7f = h1xpz2(s19mjs5d_pw_s3bmjs7f);
            setpz2(x_p+N7, subpz2(s08mjs4c_pw_s2amjs6e, h1_s19mjs5d_pw_s3bmjs7f));
            setpz2(x_p+Nf, addpz2(s08mjs4c_pw_s2amjs6e, h1_s19mjs5d_pw_s3bmjs7f));

        }
    }
};

///////////////////////////////////////////////////////////////////////////////

template <int n, int s, bool eo, int mode> struct invend;

//-----------------------------------------------------------------------------

template <int s, bool eo, int mode> struct invend<16,s,eo,mode>
{
    static constexpr int N = 16*s;

    void operator()(complex_vector x, complex_vector y) const noexcept
    {
        complex_vector z = eo ? y : x;
        for (int q = 0; q < s; q += 2) {
            complex_vector xq = x + q;
            complex_vector zq = z + q;

            const Vec4d z0 = scalepz2<N,mode>(getpz2(zq+s*0x0));
            const Vec4d z1 = scalepz2<N,mode>(getpz2(zq+s*0x1));
            const Vec4d z2 = scalepz2<N,mode>(getpz2(zq+s*0x2));
            const Vec4d z3 = scalepz2<N,mode>(getpz2(zq+s*0x3));
            const Vec4d z4 = scalepz2<N,mode>(getpz2(zq+s*0x4));
            const Vec4d z5 = scalepz2<N,mode>(getpz2(zq+s*0x5));
            const Vec4d z6 = scalepz2<N,mode>(getpz2(zq+s*0x6));
            const Vec4d z7 = scalepz2<N,mode>(getpz2(zq+s*0x7));
            const Vec4d z8 = scalepz2<N,mode>(getpz2(zq+s*0x8));
            const Vec4d z9 = scalepz2<N,mode>(getpz2(zq+s*0x9));
            const Vec4d za = scalepz2<N,mode>(getpz2(zq+s*0xa));
            const Vec4d zb = scalepz2<N,mode>(getpz2(zq+s*0xb));
            const Vec4d zc = scalepz2<N,mode>(getpz2(zq+s*0xc));
            const Vec4d zd = scalepz2<N,mode>(getpz2(zq+s*0xd));
            const Vec4d ze = scalepz2<N,mode>(getpz2(zq+s*0xe));
            const Vec4d zf = scalepz2<N,mode>(getpz2(zq+s*0xf));

            const Vec4d a08 = addpz2(z0, z8); const Vec4d s08 = subpz2(z0, z8);
            const Vec4d a4c = addpz2(z4, zc); const Vec4d s4c = subpz2(z4, zc);
            const Vec4d a2a = addpz2(z2, za); const Vec4d s2a = subpz2(z2, za);
            const Vec4d a6e = addpz2(z6, ze); const Vec4d s6e = subpz2(z6, ze);
            const Vec4d a19 = addpz2(z1, z9); const Vec4d s19 = subpz2(z1, z9);
            const Vec4d a5d = addpz2(z5, zd); const Vec4d s5d = subpz2(z5, zd);
            const Vec4d a3b = addpz2(z3, zb); const Vec4d s3b = subpz2(z3, zb);
            const Vec4d a7f = addpz2(z7, zf); const Vec4d s7f = subpz2(z7, zf);

            const Vec4d js4c = jxpz2(s4c);
            const Vec4d js6e = jxpz2(s6e);
            const Vec4d js5d = jxpz2(s5d);
            const Vec4d js7f = jxpz2(s7f);

            const Vec4d a08p1a4c = addpz2(a08, a4c); const Vec4d s08mjs4c = subpz2(s08, js4c);
            const Vec4d a08m1a4c = subpz2(a08, a4c); const Vec4d s08pjs4c = addpz2(s08, js4c);
            const Vec4d a2ap1a6e = addpz2(a2a, a6e); const Vec4d s2amjs6e = subpz2(s2a, js6e);
            const Vec4d a2am1a6e = subpz2(a2a, a6e); const Vec4d s2apjs6e = addpz2(s2a, js6e);
            const Vec4d a19p1a5d = addpz2(a19, a5d); const Vec4d s19mjs5d = subpz2(s19, js5d);
            const Vec4d a19m1a5d = subpz2(a19, a5d); const Vec4d s19pjs5d = addpz2(s19, js5d);
            const Vec4d a3bp1a7f = addpz2(a3b, a7f); const Vec4d s3bmjs7f = subpz2(s3b, js7f);
            const Vec4d a3bm1a7f = subpz2(a3b, a7f); const Vec4d s3bpjs7f = addpz2(s3b, js7f);

            const Vec4d w8_s2amjs6e = w8xpz2(s2amjs6e);
            const Vec4d  j_a2am1a6e =  jxpz2(a2am1a6e);
            const Vec4d v8_s2apjs6e = v8xpz2(s2apjs6e);

            const Vec4d a08p1a4c_p1_a2ap1a6e = addpz2(a08p1a4c,    a2ap1a6e);
            const Vec4d s08mjs4c_pw_s2amjs6e = addpz2(s08mjs4c, w8_s2amjs6e);
            const Vec4d a08m1a4c_mj_a2am1a6e = subpz2(a08m1a4c,  j_a2am1a6e);
            const Vec4d s08pjs4c_mv_s2apjs6e = subpz2(s08pjs4c, v8_s2apjs6e);
            const Vec4d a08p1a4c_m1_a2ap1a6e = subpz2(a08p1a4c,    a2ap1a6e);
            const Vec4d s08mjs4c_mw_s2amjs6e = subpz2(s08mjs4c, w8_s2amjs6e);
            const Vec4d a08m1a4c_pj_a2am1a6e = addpz2(a08m1a4c,  j_a2am1a6e);
            const Vec4d s08pjs4c_pv_s2apjs6e = addpz2(s08pjs4c, v8_s2apjs6e);

            const Vec4d w8_s3bmjs7f = w8xpz2(s3bmjs7f);
            const Vec4d  j_a3bm1a7f =  jxpz2(a3bm1a7f);
            const Vec4d v8_s3bpjs7f = v8xpz2(s3bpjs7f);

            const Vec4d a19p1a5d_p1_a3bp1a7f = addpz2(a19p1a5d,    a3bp1a7f);
            const Vec4d s19mjs5d_pw_s3bmjs7f = addpz2(s19mjs5d, w8_s3bmjs7f);
            const Vec4d a19m1a5d_mj_a3bm1a7f = subpz2(a19m1a5d,  j_a3bm1a7f);
            const Vec4d s19pjs5d_mv_s3bpjs7f = subpz2(s19pjs5d, v8_s3bpjs7f);
            const Vec4d a19p1a5d_m1_a3bp1a7f = subpz2(a19p1a5d,    a3bp1a7f);
            const Vec4d s19mjs5d_mw_s3bmjs7f = subpz2(s19mjs5d, w8_s3bmjs7f);
            const Vec4d a19m1a5d_pj_a3bm1a7f = addpz2(a19m1a5d,  j_a3bm1a7f);
            const Vec4d s19pjs5d_pv_s3bpjs7f = addpz2(s19pjs5d, v8_s3bpjs7f);

            const Vec4d h1_s19mjs5d_pw_s3bmjs7f = h1xpz2(s19mjs5d_pw_s3bmjs7f);
            const Vec4d w8_a19m1a5d_mj_a3bm1a7f = w8xpz2(a19m1a5d_mj_a3bm1a7f);
            const Vec4d h3_s19pjs5d_mv_s3bpjs7f = h3xpz2(s19pjs5d_mv_s3bpjs7f);
            const Vec4d  j_a19p1a5d_m1_a3bp1a7f =  jxpz2(a19p1a5d_m1_a3bp1a7f);
            const Vec4d hd_s19mjs5d_mw_s3bmjs7f = hdxpz2(s19mjs5d_mw_s3bmjs7f);
            const Vec4d v8_a19m1a5d_pj_a3bm1a7f = v8xpz2(a19m1a5d_pj_a3bm1a7f);
            const Vec4d hf_s19pjs5d_pv_s3bpjs7f = hfxpz2(s19pjs5d_pv_s3bpjs7f);

            setpz2(xq+s*0x0, addpz2(a08p1a4c_p1_a2ap1a6e,    a19p1a5d_p1_a3bp1a7f));
            setpz2(xq+s*0x1, addpz2(s08pjs4c_pv_s2apjs6e, hf_s19pjs5d_pv_s3bpjs7f));
            setpz2(xq+s*0x2, addpz2(a08m1a4c_pj_a2am1a6e, v8_a19m1a5d_pj_a3bm1a7f));
            setpz2(xq+s*0x3, addpz2(s08mjs4c_mw_s2amjs6e, hd_s19mjs5d_mw_s3bmjs7f));
            setpz2(xq+s*0x4, addpz2(a08p1a4c_m1_a2ap1a6e,  j_a19p1a5d_m1_a3bp1a7f));
            setpz2(xq+s*0x5, subpz2(s08pjs4c_mv_s2apjs6e, h3_s19pjs5d_mv_s3bpjs7f));
            setpz2(xq+s*0x6, subpz2(a08m1a4c_mj_a2am1a6e, w8_a19m1a5d_mj_a3bm1a7f));
            setpz2(xq+s*0x7, subpz2(s08mjs4c_pw_s2amjs6e, h1_s19mjs5d_pw_s3bmjs7f));

            setpz2(xq+s*0x8, subpz2(a08p1a4c_p1_a2ap1a6e,    a19p1a5d_p1_a3bp1a7f));
            setpz2(xq+s*0x9, subpz2(s08pjs4c_pv_s2apjs6e, hf_s19pjs5d_pv_s3bpjs7f));
            setpz2(xq+s*0xa, subpz2(a08m1a4c_pj_a2am1a6e, v8_a19m1a5d_pj_a3bm1a7f));
            setpz2(xq+s*0xb, subpz2(s08mjs4c_mw_s2amjs6e, hd_s19mjs5d_mw_s3bmjs7f));
            setpz2(xq+s*0xc, subpz2(a08p1a4c_m1_a2ap1a6e,  j_a19p1a5d_m1_a3bp1a7f));
            setpz2(xq+s*0xd, addpz2(s08pjs4c_mv_s2apjs6e, h3_s19pjs5d_mv_s3bpjs7f));
            setpz2(xq+s*0xe, addpz2(a08m1a4c_mj_a2am1a6e, w8_a19m1a5d_mj_a3bm1a7f));
            setpz2(xq+s*0xf, addpz2(s08mjs4c_pw_s2amjs6e, h1_s19mjs5d_pw_s3bmjs7f));
        }
    }
};

template <bool eo, int mode> struct invend<16,1,eo,mode>
{
    inline void operator()(complex_vector x, complex_vector y) const noexcept
    {
        complex_vector z = eo ? y : x;
        const Vec2d z0 = scalepz<16,mode>(getpz(z[0x0]));
        const Vec2d z1 = scalepz<16,mode>(getpz(z[0x1]));
        const Vec2d z2 = scalepz<16,mode>(getpz(z[0x2]));
        const Vec2d z3 = scalepz<16,mode>(getpz(z[0x3]));
        const Vec2d z4 = scalepz<16,mode>(getpz(z[0x4]));
        const Vec2d z5 = scalepz<16,mode>(getpz(z[0x5]));
        const Vec2d z6 = scalepz<16,mode>(getpz(z[0x6]));
        const Vec2d z7 = scalepz<16,mode>(getpz(z[0x7]));
        const Vec2d z8 = scalepz<16,mode>(getpz(z[0x8]));
        const Vec2d z9 = scalepz<16,mode>(getpz(z[0x9]));
        const Vec2d za = scalepz<16,mode>(getpz(z[0xa]));
        const Vec2d zb = scalepz<16,mode>(getpz(z[0xb]));
        const Vec2d zc = scalepz<16,mode>(getpz(z[0xc]));
        const Vec2d zd = scalepz<16,mode>(getpz(z[0xd]));
        const Vec2d ze = scalepz<16,mode>(getpz(z[0xe]));
        const Vec2d zf = scalepz<16,mode>(getpz(z[0xf]));

        const Vec2d a08 = addpz(z0, z8); const Vec2d s08 = subpz(z0, z8);
        const Vec2d a4c = addpz(z4, zc); const Vec2d s4c = subpz(z4, zc);
        const Vec2d a2a = addpz(z2, za); const Vec2d s2a = subpz(z2, za);
        const Vec2d a6e = addpz(z6, ze); const Vec2d s6e = subpz(z6, ze);
        const Vec2d a19 = addpz(z1, z9); const Vec2d s19 = subpz(z1, z9);
        const Vec2d a5d = addpz(z5, zd); const Vec2d s5d = subpz(z5, zd);
        const Vec2d a3b = addpz(z3, zb); const Vec2d s3b = subpz(z3, zb);
        const Vec2d a7f = addpz(z7, zf); const Vec2d s7f = subpz(z7, zf);

        const Vec2d js4c = jxpz(s4c);
        const Vec2d js6e = jxpz(s6e);
        const Vec2d js5d = jxpz(s5d);
        const Vec2d js7f = jxpz(s7f);

        const Vec2d a08p1a4c = addpz(a08, a4c); const Vec2d s08mjs4c = subpz(s08, js4c);
        const Vec2d a08m1a4c = subpz(a08, a4c); const Vec2d s08pjs4c = addpz(s08, js4c);
        const Vec2d a2ap1a6e = addpz(a2a, a6e); const Vec2d s2amjs6e = subpz(s2a, js6e);
        const Vec2d a2am1a6e = subpz(a2a, a6e); const Vec2d s2apjs6e = addpz(s2a, js6e);
        const Vec2d a19p1a5d = addpz(a19, a5d); const Vec2d s19mjs5d = subpz(s19, js5d);
        const Vec2d a19m1a5d = subpz(a19, a5d); const Vec2d s19pjs5d = addpz(s19, js5d);
        const Vec2d a3bp1a7f = addpz(a3b, a7f); const Vec2d s3bmjs7f = subpz(s3b, js7f);
        const Vec2d a3bm1a7f = subpz(a3b, a7f); const Vec2d s3bpjs7f = addpz(s3b, js7f);

        const Vec2d w8_s2amjs6e = w8xpz(s2amjs6e);
        const Vec2d  j_a2am1a6e =  jxpz(a2am1a6e);
        const Vec2d v8_s2apjs6e = v8xpz(s2apjs6e);

        const Vec2d a08p1a4c_p1_a2ap1a6e = addpz(a08p1a4c,    a2ap1a6e);
        const Vec2d s08mjs4c_pw_s2amjs6e = addpz(s08mjs4c, w8_s2amjs6e);
        const Vec2d a08m1a4c_mj_a2am1a6e = subpz(a08m1a4c,  j_a2am1a6e);
        const Vec2d s08pjs4c_mv_s2apjs6e = subpz(s08pjs4c, v8_s2apjs6e);
        const Vec2d a08p1a4c_m1_a2ap1a6e = subpz(a08p1a4c,    a2ap1a6e);
        const Vec2d s08mjs4c_mw_s2amjs6e = subpz(s08mjs4c, w8_s2amjs6e);
        const Vec2d a08m1a4c_pj_a2am1a6e = addpz(a08m1a4c,  j_a2am1a6e);
        const Vec2d s08pjs4c_pv_s2apjs6e = addpz(s08pjs4c, v8_s2apjs6e);

        const Vec2d w8_s3bmjs7f = w8xpz(s3bmjs7f);
        const Vec2d  j_a3bm1a7f =  jxpz(a3bm1a7f);
        const Vec2d v8_s3bpjs7f = v8xpz(s3bpjs7f);

        const Vec2d a19p1a5d_p1_a3bp1a7f = addpz(a19p1a5d,    a3bp1a7f);
        const Vec2d s19mjs5d_pw_s3bmjs7f = addpz(s19mjs5d, w8_s3bmjs7f);
        const Vec2d a19m1a5d_mj_a3bm1a7f = subpz(a19m1a5d,  j_a3bm1a7f);
        const Vec2d s19pjs5d_mv_s3bpjs7f = subpz(s19pjs5d, v8_s3bpjs7f);
        const Vec2d a19p1a5d_m1_a3bp1a7f = subpz(a19p1a5d,    a3bp1a7f);
        const Vec2d s19mjs5d_mw_s3bmjs7f = subpz(s19mjs5d, w8_s3bmjs7f);
        const Vec2d a19m1a5d_pj_a3bm1a7f = addpz(a19m1a5d,  j_a3bm1a7f);
        const Vec2d s19pjs5d_pv_s3bpjs7f = addpz(s19pjs5d, v8_s3bpjs7f);

        const Vec2d h1_s19mjs5d_pw_s3bmjs7f = h1xpz(s19mjs5d_pw_s3bmjs7f);
        const Vec2d w8_a19m1a5d_mj_a3bm1a7f = w8xpz(a19m1a5d_mj_a3bm1a7f);
        const Vec2d h3_s19pjs5d_mv_s3bpjs7f = h3xpz(s19pjs5d_mv_s3bpjs7f);
        const Vec2d  j_a19p1a5d_m1_a3bp1a7f =  jxpz(a19p1a5d_m1_a3bp1a7f);
        const Vec2d hd_s19mjs5d_mw_s3bmjs7f = hdxpz(s19mjs5d_mw_s3bmjs7f);
        const Vec2d v8_a19m1a5d_pj_a3bm1a7f = v8xpz(a19m1a5d_pj_a3bm1a7f);
        const Vec2d hf_s19pjs5d_pv_s3bpjs7f = hfxpz(s19pjs5d_pv_s3bpjs7f);

        setpz(x[0x0], addpz(a08p1a4c_p1_a2ap1a6e,    a19p1a5d_p1_a3bp1a7f));
        setpz(x[0x1], addpz(s08pjs4c_pv_s2apjs6e, hf_s19pjs5d_pv_s3bpjs7f));
        setpz(x[0x2], addpz(a08m1a4c_pj_a2am1a6e, v8_a19m1a5d_pj_a3bm1a7f));
        setpz(x[0x3], addpz(s08mjs4c_mw_s2amjs6e, hd_s19mjs5d_mw_s3bmjs7f));
        setpz(x[0x4], addpz(a08p1a4c_m1_a2ap1a6e,  j_a19p1a5d_m1_a3bp1a7f));
        setpz(x[0x5], subpz(s08pjs4c_mv_s2apjs6e, h3_s19pjs5d_mv_s3bpjs7f));
        setpz(x[0x6], subpz(a08m1a4c_mj_a2am1a6e, w8_a19m1a5d_mj_a3bm1a7f));
        setpz(x[0x7], subpz(s08mjs4c_pw_s2amjs6e, h1_s19mjs5d_pw_s3bmjs7f));

        setpz(x[0x8], subpz(a08p1a4c_p1_a2ap1a6e,    a19p1a5d_p1_a3bp1a7f));
        setpz(x[0x9], subpz(s08pjs4c_pv_s2apjs6e, hf_s19pjs5d_pv_s3bpjs7f));
        setpz(x[0xa], subpz(a08m1a4c_pj_a2am1a6e, v8_a19m1a5d_pj_a3bm1a7f));
        setpz(x[0xb], subpz(s08mjs4c_mw_s2amjs6e, hd_s19mjs5d_mw_s3bmjs7f));
        setpz(x[0xc], subpz(a08p1a4c_m1_a2ap1a6e,  j_a19p1a5d_m1_a3bp1a7f));
        setpz(x[0xd], addpz(s08pjs4c_mv_s2apjs6e, h3_s19pjs5d_mv_s3bpjs7f));
        setpz(x[0xe], addpz(a08m1a4c_mj_a2am1a6e, w8_a19m1a5d_mj_a3bm1a7f));
        setpz(x[0xf], addpz(s08mjs4c_pw_s2amjs6e, h1_s19mjs5d_pw_s3bmjs7f));
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
        invfft<n/16,16*s,!eo,mode>()(y, x, W);
        invcore<n,s>()(x, y, W);
    }
};

template <int s, bool eo, int mode> struct invfft<16,s,eo,mode>
{
    inline void operator()(
        complex_vector x, complex_vector y, const_complex_vector) const noexcept
    {
        invend<16,s,eo,mode>()(x, y);
    }
};

template <int s, bool eo, int mode> struct invfft<8,s,eo,mode>
{
    inline void operator()(
        complex_vector x, complex_vector y, const_complex_vector) const noexcept
    {
        OTFFT_AVXDIT8::invend<8,s,eo,mode>()(x, y);
    }
};

template <int s, bool eo, int mode> struct invfft<4,s,eo,mode>
{
    inline void operator()(
        complex_vector x, complex_vector y, const_complex_vector) const noexcept
    {
        OTFFT_AVXDIT4::invend<4,s,eo,mode>()(x, y);
    }
};

template <int s, bool eo, int mode> struct invfft<2,s,eo,mode>
{
    inline void operator()(
        complex_vector x, complex_vector y, const_complex_vector) const noexcept
    {
        OTFFT_AVXDIT4::invend<2,s,eo,mode>()(x, y);
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
        if (N <= 16) W = 0;
        else {
            weight.setup(2*N);
            W = weight.get();
            init_Wt(16, N, W);
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

#endif // otfft_avxdit16_h
