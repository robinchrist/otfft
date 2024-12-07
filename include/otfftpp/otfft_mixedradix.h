/******************************************************************************
*  OTFFT Mixed Radix Version 11.5e
*
*  Copyright (c) 2016 OK Ojisan(Takuya OKAHISA)
*  Released under the MIT license
*  http://opensource.org/licenses/mit-license.php
******************************************************************************/

#ifndef otfft_mixedradix_h
#define otfft_mixedradix_h

#include <iostream>
#include "otfft_misc.h"

namespace OTFFT_MixedRadix { //////////////////////////////////////////////////

using namespace OTFFT_MISC;

/*struct cpx {
    Vec2d z;
    cpx(const Vec2d& z) noexcept : z(z) {}
    operator Vec2d() const noexcept { return z; }
};

static inline cpx operator+(const cpx a, const cpx b) noexcept
{
    return addpz(a, b);
}
static inline cpx operator*(const cpx a, const cpx b) noexcept
{
    return mulpz(a, b);
}*/
using cpx = enoki::Complex<float>;

///////////////////////////////////////////////////////////////////////////////
// Forward Butterfly Operation
///////////////////////////////////////////////////////////////////////////////

void fwdend2(const int s, const bool eo,
        complex_vector x, complex_vector y) noexcept
{
    complex_vector z = eo ? y : x;
    if (s >= 2) {
        for (int q = 0; q < s; q += 2) {
            complex_vector xq = x + q;
            complex_vector zq = z + q;
            const Vec4d a = getpz2(xq+0);
            const Vec4d b = getpz2(xq+s);
            setpz2(zq+0, addpz2(a, b));
            setpz2(zq+s, subpz2(a, b));
        }
    }
    else {
        const Vec2d a = getpz(x[0]);
        const Vec2d b = getpz(x[1]);
        setpz(z[0], addpz(a, b));
        setpz(z[1], subpz(a, b));
    }
}

void fwdcore2(const int n, const int s,
        complex_vector x, complex_vector y, const_complex_vector W) noexcept
{
    const int m  = n/2;
    const int N  = n*s;
    const int N0 = 0;
    const int N1 = N/2;
    if (s >= 2) {
        for (int p = 0; p < m; p++) {
            const int sp = s*p;
            const int s2p = 2*sp;
            const Vec4d wp = duppz3(W[sp]);
            for (int q = 0; q < s; q += 2) {
                complex_vector xq_sp  = x + q + sp;
                complex_vector yq_s2p = y + q + s2p;
                const Vec4d a = getpz2(xq_sp+N0);
                const Vec4d b = getpz2(xq_sp+N1);
                setpz2(yq_s2p+s*0,            addpz2(a, b));
                setpz2(yq_s2p+s*1, mulpz2(wp, subpz2(a, b)));
            }
        }
    }
    else {
        for (int p = 0; p < m; p++) {
            complex_vector x_p  = x + p;
            complex_vector y_2p = y + 2*p;
            const Vec2d wp = getpz(W[p]);
            const Vec2d a = getpz(x_p[N0]);
            const Vec2d b = getpz(x_p[N1]);
            setpz(y_2p[0],           addpz(a, b));
            setpz(y_2p[1], mulpz(wp, subpz(a, b)));
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

void fwdend4(const int s, const bool eo,
        complex_vector x, complex_vector y) noexcept
{
    complex_vector z = eo ? y : x;
    if (s >= 2) {
        for (int q = 0; q < s; q += 2) {
            complex_vector xq = x + q;
            complex_vector zq = z + q;
            const Vec4d a = getpz2(xq+s*0);
            const Vec4d b = getpz2(xq+s*1);
            const Vec4d c = getpz2(xq+s*2);
            const Vec4d d = getpz2(xq+s*3);
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
    else {
        const Vec2d a = getpz(x[0]);
        const Vec2d b = getpz(x[1]);
        const Vec2d c = getpz(x[2]);
        const Vec2d d = getpz(x[3]);
        const Vec2d  apc =      addpz(a, c);
        const Vec2d  amc =      subpz(a, c);
        const Vec2d  bpd =      addpz(b, d);
        const Vec2d jbmd = jxpz(subpz(b, d));
        setpz(z[0], addpz(apc,  bpd));
        setpz(z[1], subpz(amc, jbmd));
        setpz(z[2], subpz(apc,  bpd));
        setpz(z[3], addpz(amc, jbmd));
    }
}

void fwdcore4(const int n, const int s,
        complex_vector x, complex_vector y, const_complex_vector W) noexcept
{
    const int m  = n/4;
    const int N  = n*s;
    const int N0 = 0;
    const int N1 = N/4;
    const int N2 = N1*2;
    const int N3 = N1*3;
    if (s >= 2) {
        for (int p = 0; p < m; p++) {
            const int sp = s*p;
            const int s4p = 4*sp;
            const Vec4d w1p = duppz3(W[1*sp]);
            const Vec4d w2p = duppz3(W[2*sp]);
            const Vec4d w3p = duppz3(W[3*sp]);
            for (int q = 0; q < s; q += 2) {
                complex_vector xq_sp  = x + q + sp;
                complex_vector yq_s4p = y + q + s4p;
                const Vec4d a = getpz2(xq_sp+N0);
                const Vec4d b = getpz2(xq_sp+N1);
                const Vec4d c = getpz2(xq_sp+N2);
                const Vec4d d = getpz2(xq_sp+N3);
                const Vec4d  apc =       addpz2(a, c);
                const Vec4d  amc =       subpz2(a, c);
                const Vec4d  bpd =       addpz2(b, d);
                const Vec4d jbmd = jxpz2(subpz2(b, d));
                setpz2(yq_s4p+s*0,             addpz2(apc,  bpd));
                setpz2(yq_s4p+s*1, mulpz2(w1p, subpz2(amc, jbmd)));
                setpz2(yq_s4p+s*2, mulpz2(w2p, subpz2(apc,  bpd)));
                setpz2(yq_s4p+s*3, mulpz2(w3p, addpz2(amc, jbmd)));
            }
        }
    }
    else {
        for (int p = 0; p < m; p++) {
            complex_vector x_p  = x + p;
            complex_vector y_4p = y + 4*p;
            const Vec2d w1p = getpz(W[p]);
            const Vec2d w2p = mulpz(w1p,w1p);
            const Vec2d w3p = mulpz(w1p,w2p);
            const Vec2d a = getpz(x_p[N0]);
            const Vec2d b = getpz(x_p[N1]);
            const Vec2d c = getpz(x_p[N2]);
            const Vec2d d = getpz(x_p[N3]);
            const Vec2d  apc =      addpz(a, c);
            const Vec2d  amc =      subpz(a, c);
            const Vec2d  bpd =      addpz(b, d);
            const Vec2d jbmd = jxpz(subpz(b, d));
            setpz(y_4p[0],            addpz(apc,  bpd));
            setpz(y_4p[1], mulpz(w1p, subpz(amc, jbmd)));
            setpz(y_4p[2], mulpz(w2p, subpz(apc,  bpd)));
            setpz(y_4p[3], mulpz(w3p, addpz(amc, jbmd)));
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

void fwdend8(const int s, const bool eo,
        complex_vector x, complex_vector y) noexcept
{
    complex_vector z = eo ? y : x;
    if (s >= 2) {
        for (int q = 0; q < s; q += 2) {
            complex_vector xq = x + q;
            complex_vector zq = z + q;
            const Vec4d x0 = getpz2(xq+s*0);
            const Vec4d x1 = getpz2(xq+s*1);
            const Vec4d x2 = getpz2(xq+s*2);
            const Vec4d x3 = getpz2(xq+s*3);
            const Vec4d x4 = getpz2(xq+s*4);
            const Vec4d x5 = getpz2(xq+s*5);
            const Vec4d x6 = getpz2(xq+s*6);
            const Vec4d x7 = getpz2(xq+s*7);
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
            setpz2(zq+s*0, addpz2(a04_p1_a26,    a15_p1_a37));
            setpz2(zq+s*1, addpz2(s04_mj_s26, w8_s15_mj_s37));
            setpz2(zq+s*2, subpz2(a04_m1_a26,  j_a15_m1_a37));
            setpz2(zq+s*3, subpz2(s04_pj_s26, v8_s15_pj_s37));
            setpz2(zq+s*4, subpz2(a04_p1_a26,    a15_p1_a37));
            setpz2(zq+s*5, subpz2(s04_mj_s26, w8_s15_mj_s37));
            setpz2(zq+s*6, addpz2(a04_m1_a26,  j_a15_m1_a37));
            setpz2(zq+s*7, addpz2(s04_pj_s26, v8_s15_pj_s37));
        }
    }
    else {
        const Vec2d x0 = getpz(x[0]);
        const Vec2d x1 = getpz(x[1]);
        const Vec2d x2 = getpz(x[2]);
        const Vec2d x3 = getpz(x[3]);
        const Vec2d x4 = getpz(x[4]);
        const Vec2d x5 = getpz(x[5]);
        const Vec2d x6 = getpz(x[6]);
        const Vec2d x7 = getpz(x[7]);
        const Vec2d  a04 =      addpz(x0, x4);
        const Vec2d  s04 =      subpz(x0, x4);
        const Vec2d  a26 =      addpz(x2, x6);
        const Vec2d js26 = jxpz(subpz(x2, x6));
        const Vec2d  a15 =      addpz(x1, x5);
        const Vec2d  s15 =      subpz(x1, x5);
        const Vec2d  a37 =      addpz(x3, x7);
        const Vec2d js37 = jxpz(subpz(x3, x7));
        const Vec2d    a04_p1_a26 =       addpz(a04,  a26);
        const Vec2d    s04_mj_s26 =       subpz(s04, js26);
        const Vec2d    a04_m1_a26 =       subpz(a04,  a26);
        const Vec2d    s04_pj_s26 =       addpz(s04, js26);
        const Vec2d    a15_p1_a37 =       addpz(a15,  a37);
        const Vec2d w8_s15_mj_s37 = w8xpz(subpz(s15, js37));
        const Vec2d  j_a15_m1_a37 =  jxpz(subpz(a15,  a37));
        const Vec2d v8_s15_pj_s37 = v8xpz(addpz(s15, js37));
        setpz(z[0], addpz(a04_p1_a26,    a15_p1_a37));
        setpz(z[1], addpz(s04_mj_s26, w8_s15_mj_s37));
        setpz(z[2], subpz(a04_m1_a26,  j_a15_m1_a37));
        setpz(z[3], subpz(s04_pj_s26, v8_s15_pj_s37));
        setpz(z[4], subpz(a04_p1_a26,    a15_p1_a37));
        setpz(z[5], subpz(s04_mj_s26, w8_s15_mj_s37));
        setpz(z[6], addpz(a04_m1_a26,  j_a15_m1_a37));
        setpz(z[7], addpz(s04_pj_s26, v8_s15_pj_s37));
    }
}

void fwdcore8(const int n, const int s,
        complex_vector x, complex_vector y, const_complex_vector W) noexcept
{
    const int m  = n/8;
    const int N  = n*s;
    const int N0 = 0;
    const int N1 = N/8;
    const int N2 = N1*2;
    const int N3 = N1*3;
    const int N4 = N1*4;
    const int N5 = N1*5;
    const int N6 = N1*6;
    const int N7 = N1*7;
    if (s >= 2) {
        for (int p = 0; p < m; p++) {
            const int sp = s*p;
            const int s8p = 8*sp;
            const Vec4d w1p = duppz3(W[1*sp]);
            const Vec4d w2p = duppz3(W[2*sp]);
            const Vec4d w3p = duppz3(W[3*sp]);
            const Vec4d w4p = mulpz2(w2p,w2p);
            const Vec4d w5p = mulpz2(w2p,w3p);
            const Vec4d w6p = mulpz2(w3p,w3p);
            const Vec4d w7p = mulpz2(w3p,w4p);
            for (int q = 0; q < s; q += 2) {
                complex_vector xq_sp  = x + q + sp;
                complex_vector yq_s8p = y + q + s8p;
                const Vec4d x0 = getpz2(xq_sp+N0);
                const Vec4d x1 = getpz2(xq_sp+N1);
                const Vec4d x2 = getpz2(xq_sp+N2);
                const Vec4d x3 = getpz2(xq_sp+N3);
                const Vec4d x4 = getpz2(xq_sp+N4);
                const Vec4d x5 = getpz2(xq_sp+N5);
                const Vec4d x6 = getpz2(xq_sp+N6);
                const Vec4d x7 = getpz2(xq_sp+N7);
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
                setpz2(yq_s8p+s*0,             addpz2(a04_p1_a26,    a15_p1_a37));
                setpz2(yq_s8p+s*1, mulpz2(w1p, addpz2(s04_mj_s26, w8_s15_mj_s37)));
                setpz2(yq_s8p+s*2, mulpz2(w2p, subpz2(a04_m1_a26,  j_a15_m1_a37)));
                setpz2(yq_s8p+s*3, mulpz2(w3p, subpz2(s04_pj_s26, v8_s15_pj_s37)));
                setpz2(yq_s8p+s*4, mulpz2(w4p, subpz2(a04_p1_a26,    a15_p1_a37)));
                setpz2(yq_s8p+s*5, mulpz2(w5p, subpz2(s04_mj_s26, w8_s15_mj_s37)));
                setpz2(yq_s8p+s*6, mulpz2(w6p, addpz2(a04_m1_a26,  j_a15_m1_a37)));
                setpz2(yq_s8p+s*7, mulpz2(w7p, addpz2(s04_pj_s26, v8_s15_pj_s37)));
            }
        }
    }
    else {
        for (int p = 0; p < m; p++) {
            complex_vector x_p  = x + p;
            complex_vector y_8p = y + 8*p;
            const Vec2d w1p = getpz(W[p]);
            const Vec2d w2p = mulpz(w1p,w1p);
            const Vec2d w3p = mulpz(w1p,w2p);
            const Vec2d w4p = mulpz(w2p,w2p);
            const Vec2d w5p = mulpz(w2p,w3p);
            const Vec2d w6p = mulpz(w3p,w3p);
            const Vec2d w7p = mulpz(w3p,w4p);
            const Vec2d x0 = getpz(x_p[N0]);
            const Vec2d x1 = getpz(x_p[N1]);
            const Vec2d x2 = getpz(x_p[N2]);
            const Vec2d x3 = getpz(x_p[N3]);
            const Vec2d x4 = getpz(x_p[N4]);
            const Vec2d x5 = getpz(x_p[N5]);
            const Vec2d x6 = getpz(x_p[N6]);
            const Vec2d x7 = getpz(x_p[N7]);
            const Vec2d  a04 =      addpz(x0, x4);
            const Vec2d  s04 =      subpz(x0, x4);
            const Vec2d  a26 =      addpz(x2, x6);
            const Vec2d js26 = jxpz(subpz(x2, x6));
            const Vec2d  a15 =      addpz(x1, x5);
            const Vec2d  s15 =      subpz(x1, x5);
            const Vec2d  a37 =      addpz(x3, x7);
            const Vec2d js37 = jxpz(subpz(x3, x7));
            const Vec2d    a04_p1_a26 =       addpz(a04,  a26);
            const Vec2d    s04_mj_s26 =       subpz(s04, js26);
            const Vec2d    a04_m1_a26 =       subpz(a04,  a26);
            const Vec2d    s04_pj_s26 =       addpz(s04, js26);
            const Vec2d    a15_p1_a37 =       addpz(a15,  a37);
            const Vec2d w8_s15_mj_s37 = w8xpz(subpz(s15, js37));
            const Vec2d  j_a15_m1_a37 =  jxpz(subpz(a15,  a37));
            const Vec2d v8_s15_pj_s37 = v8xpz(addpz(s15, js37));
            setpz(y_8p[0],            addpz(a04_p1_a26,    a15_p1_a37));
            setpz(y_8p[1], mulpz(w1p, addpz(s04_mj_s26, w8_s15_mj_s37)));
            setpz(y_8p[2], mulpz(w2p, subpz(a04_m1_a26,  j_a15_m1_a37)));
            setpz(y_8p[3], mulpz(w3p, subpz(s04_pj_s26, v8_s15_pj_s37)));
            setpz(y_8p[4], mulpz(w4p, subpz(a04_p1_a26,    a15_p1_a37)));
            setpz(y_8p[5], mulpz(w5p, subpz(s04_mj_s26, w8_s15_mj_s37)));
            setpz(y_8p[6], mulpz(w6p, addpz(a04_m1_a26,  j_a15_m1_a37)));
            setpz(y_8p[7], mulpz(w7p, addpz(s04_pj_s26, v8_s15_pj_s37)));
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

void fwdend5(const int s, const bool eo,
        complex_vector x, complex_vector y, const_complex_vector W) noexcept
{
    const cpx w1 = getpz(W[s]);
    const cpx w2 = mulpz(w1,w1);
    const cpx w3 = mulpz(w1,w2);
    const cpx w4 = mulpz(w2,w2);
    complex_vector z = eo ? y : x;
    for (int q = 0; q < s; q++) {
        const cpx a = getpz(x[q+s*0]);
        const cpx b = getpz(x[q+s*1]);
        const cpx c = getpz(x[q+s*2]);
        const cpx d = getpz(x[q+s*3]);
        const cpx e = getpz(x[q+s*4]);
        setpz(z[q+s*0], a + b + c + d + e);
        setpz(z[q+s*1], a + w1*b + w2*c + w3*d + w4*e);
        setpz(z[q+s*2], a + w2*b + w4*c + w1*d + w3*e);
        setpz(z[q+s*3], a + w3*b + w1*c + w4*d + w2*e);
        setpz(z[q+s*4], a + w4*b + w3*c + w2*d + w1*e);
    }
}

void fwdcore5(const int n, const int s,
        complex_vector x, complex_vector y, const_complex_vector W) noexcept
{
    const int N  = n*s;
    const int m  = n/5;
    const int N0 = 0;
    const int N1 = N/5;
    const int N2 = N1*2;
    const int N3 = N1*3;
    const int N4 = N1*4;
    const cpx w1 = getpz(W[N1]);
    const cpx w2 = mulpz(w1,w1);
    const cpx w3 = mulpz(w1,w2);
    const cpx w4 = mulpz(w2,w2);
    for (int p = 0; p < m; p++) {
        const int sp = s*p;
        const cpx w1p = getpz(W[sp]);
        const cpx w2p = mulpz(w1p,w1p);
        const cpx w3p = mulpz(w1p,w2p);
        const cpx w4p = mulpz(w2p,w2p);
        for (int q = 0; q < s; q++) {
            const int q_sp = q + sp;
            const cpx a = getpz(x[q_sp+N0]);
            const cpx b = getpz(x[q_sp+N1]);
            const cpx c = getpz(x[q_sp+N2]);
            const cpx d = getpz(x[q_sp+N3]);
            const cpx e = getpz(x[q_sp+N4]);
            const int q_s5p = q + sp*5;
            setpz(y[q_s5p+s*0],  a + b + c + d + e);
            setpz(y[q_s5p+s*1], (a + w1*b + w2*c + w3*d + w4*e)*w1p);
            setpz(y[q_s5p+s*2], (a + w2*b + w4*c + w1*d + w3*e)*w2p);
            setpz(y[q_s5p+s*3], (a + w3*b + w1*c + w4*d + w2*e)*w3p);
            setpz(y[q_s5p+s*4], (a + w4*b + w3*c + w2*d + w1*e)*w4p);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

void fwdend3(const int s, const bool eo,
        complex_vector x, complex_vector y, const_complex_vector W) noexcept
{
    const cpx w1 = getpz(W[s]);
    const cpx w2 = mulpz(w1,w1);
    complex_vector z = eo ? y : x;
    for (int q = 0; q < s; q++) {
        const cpx a = getpz(x[q+s*0]);
        const cpx b = getpz(x[q+s*1]);
        const cpx c = getpz(x[q+s*2]);
        setpz(z[q+s*0], a + b + c);
        setpz(z[q+s*1], a + w1*b + w2*c);
        setpz(z[q+s*2], a + w2*b + w1*c);
    }
}

void fwdcore3(const int n, const int s,
        complex_vector x, complex_vector y, const_complex_vector W) noexcept
{
    const int N  = n*s;
    const int m  = n/3;
    const int N0 = 0;
    const int N1 = N/3;
    const int N2 = N1*2;
    const cpx w1 = getpz(W[N1]);
    const cpx w2 = mulpz(w1,w1);
    for (int p = 0; p < m; p++) {
        const int sp = s*p;
        const cpx w1p = getpz(W[sp]);
        const cpx w2p = mulpz(w1p,w1p);
        for (int q = 0; q < s; q++) {
            const int q_sp = q + sp;
            const cpx a = getpz(x[q_sp+N0]);
            const cpx b = getpz(x[q_sp+N1]);
            const cpx c = getpz(x[q_sp+N2]);
            const int q_s3p = q + sp*3;
            setpz(y[q_s3p+s*0],  a + b + c);
            setpz(y[q_s3p+s*1], (a + w1*b + w2*c)*w1p);
            setpz(y[q_s3p+s*2], (a + w2*b + w1*c)*w2p);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
// Any Size FFT except Radix-2,3,5
///////////////////////////////////////////////////////////////////////////////

void fwdfftany(const int r, const int n, const int s, const bool eo,
        complex_vector x, complex_vector y, const_complex_vector W) noexcept
{
    const Vec2d zero { 0, 0 };
    const int N = n*s;
    int k = r;
    while (n%k != 0) {
        if (k*k > n) { k = n; break; }
        k += 2;
    }
    if (k == n) {
        for (int q = 0; q < s; q++) {
            for (int i = 0; i < k; i++) {
                cpx z = zero;
                for (int j = 0; j < k; j++) {
                    const cpx a   = getpz(x[q+s*j]);
                    const cpx wij = getpz(W[s*((i*j)%k)]);
                    z = z + a*wij;
                }
                setpz(y[q+s*i], z);
            }
        }
        if (!eo) for (int p = 0; p < N; p++) setpz(x[p], getpz(y[p]));
    }
    else {
        const int m  = n/k;
        const int ms = m*s;
        for (int p = 0; p < m; p++) {
            const int sp = s*p;
            for (int q = 0; q < s; q++) {
                const int q_sp  = q + sp;
                const int q_spk = q + sp*k;
                for (int i = 0; i < k; i++) {
                    cpx z = zero;
                    for (int j = 0; j < k; j++) {
                        const cpx a   = getpz(x[q_sp+ms*j]);
                        const cpx wij = getpz(W[ms*((i*j)%k)]);
                        z = z + a*wij;
                    }
                    const cpx wip = getpz(W[i*sp]);
                    setpz(y[q_spk+s*i], z * wip);
                }
            }
        }
        fwdfftany(k, m, k*s, !eo, y, x, W);
    }
}

///////////////////////////////////////////////////////////////////////////////
// Mixed Radix FFT
///////////////////////////////////////////////////////////////////////////////

void fwdfft(const int n, const int s, const bool eo,
        complex_vector x, complex_vector y, const_complex_vector W) noexcept
{
    const int N = n*s;
    if (N < 2) return;
    if (n%8 == 0) {
        if (n == 8)
            fwdend8(s, eo, x, y);
        else {
            fwdcore8(n, s, x, y, W);
            fwdfft(n/8, 8*s, !eo, y, x, W);
        }
    }
    else if (n%4 == 0) {
        if (n == 4)
            fwdend4(s, eo, x, y);
        else {
            fwdcore4(n, s, x, y, W);
            fwdfft(n/4, 4*s, !eo, y, x, W);
        }
    }
    else if (n%2 == 0) {
        if (n == 2)
            fwdend2(s, eo, x, y);
        else {
            fwdcore2(n, s, x, y, W);
            fwdfft(n/2, 2*s, !eo, y, x, W);
        }
    }
    else if (n%5 == 0) {
        if (n == 5)
            fwdend5(s, eo, x, y, W);
        else {
            fwdcore5(n, s, x, y, W);
            fwdfft(n/5, 5*s, !eo, y, x, W);
        }
    }
    else if (n%3 == 0) {
        if (n == 3)
            fwdend3(s, eo, x, y, W);
        else {
            fwdcore3(n, s, x, y, W);
            fwdfft(n/3, 3*s, !eo, y, x, W);
        }
    }
    else fwdfftany(7, n, s, eo, x, y, W);
}

///////////////////////////////////////////////////////////////////////////////
// Inverse Butterfly Operation
///////////////////////////////////////////////////////////////////////////////

void invend2(const int s, const bool eo,
        complex_vector x, complex_vector y) noexcept
{
    complex_vector z = eo ? y : x;
    if (s >= 2) {
        for (int q = 0; q < s; q += 2) {
            complex_vector xq = x + q;
            complex_vector zq = z + q;
            const Vec4d a = getpz2(xq+0);
            const Vec4d b = getpz2(xq+s);
            setpz2(zq+0, addpz2(a, b));
            setpz2(zq+s, subpz2(a, b));
        }
    }
    else {
        const Vec2d a = getpz(x[0]);
        const Vec2d b = getpz(x[1]);
        setpz(z[0], addpz(a, b));
        setpz(z[1], subpz(a, b));
    }
}

void invcore2(const int n, const int s,
        complex_vector x, complex_vector y, const_complex_vector W) noexcept
{
    const int m  = n/2;
    const int N  = n*s;
    const int N0 = 0;
    const int N1 = N/2;
    if (s >= 2) {
        for (int p = 0; p < m; p++) {
            const int sp = s*p;
            const int s2p = 2*sp;
            const Vec4d wp = duppz3(W[N-sp]);
            for (int q = 0; q < s; q += 2) {
                complex_vector xq_sp  = x + q + sp;
                complex_vector yq_s2p = y + q + s2p;
                const Vec4d a = getpz2(xq_sp+N0);
                const Vec4d b = getpz2(xq_sp+N1);
                setpz2(yq_s2p+s*0,            addpz2(a, b));
                setpz2(yq_s2p+s*1, mulpz2(wp, subpz2(a, b)));
            }
        }
    }
    else {
        for (int p = 0; p < m; p++) {
            complex_vector x_p  = x + p;
            complex_vector y_2p = y + 2*p;
            const Vec2d wp = getpz(W[N-p]);
            const Vec2d a = getpz(x_p[N0]);
            const Vec2d b = getpz(x_p[N1]);
            setpz(y_2p[0],           addpz(a, b));
            setpz(y_2p[1], mulpz(wp, subpz(a, b)));
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

void invend4(const int s, const bool eo,
        complex_vector x, complex_vector y) noexcept
{
    complex_vector z = eo ? y : x;
    if (s >= 2) {
        for (int q = 0; q < s; q += 2) {
            complex_vector xq = x + q;
            complex_vector zq = z + q;
            const Vec4d a = getpz2(xq+s*0);
            const Vec4d b = getpz2(xq+s*1);
            const Vec4d c = getpz2(xq+s*2);
            const Vec4d d = getpz2(xq+s*3);
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
    else {
        const Vec2d a = getpz(x[0]);
        const Vec2d b = getpz(x[1]);
        const Vec2d c = getpz(x[2]);
        const Vec2d d = getpz(x[3]);
        const Vec2d  apc =      addpz(a, c);
        const Vec2d  amc =      subpz(a, c);
        const Vec2d  bpd =      addpz(b, d);
        const Vec2d jbmd = jxpz(subpz(b, d));
        setpz(z[0], addpz(apc,  bpd));
        setpz(z[1], addpz(amc, jbmd));
        setpz(z[2], subpz(apc,  bpd));
        setpz(z[3], subpz(amc, jbmd));
    }
}

void invcore4(const int n, const int s,
        complex_vector x, complex_vector y, const_complex_vector W) noexcept
{
    const int m  = n/4;
    const int N  = n*s;
    const int N0 = 0;
    const int N1 = N/4;
    const int N2 = N1*2;
    const int N3 = N1*3;
    if (s >= 2) {
        for (int p = 0; p < m; p++) {
            const int sp = s*p;
            const int s4p = 4*sp;
            const Vec4d w1p = duppz3(W[N-1*sp]);
            const Vec4d w2p = duppz3(W[N-2*sp]);
            const Vec4d w3p = duppz3(W[N-3*sp]);
            for (int q = 0; q < s; q += 2) {
                complex_vector xq_sp  = x + q + sp;
                complex_vector yq_s4p = y + q + s4p;
                const Vec4d a = getpz2(xq_sp+N0);
                const Vec4d b = getpz2(xq_sp+N1);
                const Vec4d c = getpz2(xq_sp+N2);
                const Vec4d d = getpz2(xq_sp+N3);
                const Vec4d  apc =       addpz2(a, c);
                const Vec4d  amc =       subpz2(a, c);
                const Vec4d  bpd =       addpz2(b, d);
                const Vec4d jbmd = jxpz2(subpz2(b, d));
                setpz2(yq_s4p+s*0,             addpz2(apc,  bpd));
                setpz2(yq_s4p+s*1, mulpz2(w1p, addpz2(amc, jbmd)));
                setpz2(yq_s4p+s*2, mulpz2(w2p, subpz2(apc,  bpd)));
                setpz2(yq_s4p+s*3, mulpz2(w3p, subpz2(amc, jbmd)));
            }
        }
    }
    else {
        for (int p = 0; p < m; p++) {
            complex_vector x_p  = x + p;
            complex_vector y_4p = y + 4*p;
            const Vec2d w1p = cnjpz(getpz(W[p]));
            const Vec2d w2p = mulpz(w1p,w1p);
            const Vec2d w3p = mulpz(w1p,w2p);
            const Vec2d a = getpz(x_p[N0]);
            const Vec2d b = getpz(x_p[N1]);
            const Vec2d c = getpz(x_p[N2]);
            const Vec2d d = getpz(x_p[N3]);
            const Vec2d  apc =      addpz(a, c);
            const Vec2d  amc =      subpz(a, c);
            const Vec2d  bpd =      addpz(b, d);
            const Vec2d jbmd = jxpz(subpz(b, d));
            setpz(y_4p[0],            addpz(apc,  bpd));
            setpz(y_4p[1], mulpz(w1p, addpz(amc, jbmd)));
            setpz(y_4p[2], mulpz(w2p, subpz(apc,  bpd)));
            setpz(y_4p[3], mulpz(w3p, subpz(amc, jbmd)));
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

void invend8(const int s, const bool eo,
        complex_vector x, complex_vector y) noexcept
{
    complex_vector z = eo ? y : x;
    if (s >= 2) {
        for (int q = 0; q < s; q += 2) {
            complex_vector xq = x + q;
            complex_vector zq = z + q;
            const Vec4d x0 = getpz2(xq+s*0);
            const Vec4d x1 = getpz2(xq+s*1);
            const Vec4d x2 = getpz2(xq+s*2);
            const Vec4d x3 = getpz2(xq+s*3);
            const Vec4d x4 = getpz2(xq+s*4);
            const Vec4d x5 = getpz2(xq+s*5);
            const Vec4d x6 = getpz2(xq+s*6);
            const Vec4d x7 = getpz2(xq+s*7);
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
            setpz2(zq+s*0, addpz2(a04_p1_a26,    a15_p1_a37));
            setpz2(zq+s*1, addpz2(s04_pj_s26, v8_s15_pj_s37));
            setpz2(zq+s*2, addpz2(a04_m1_a26,  j_a15_m1_a37));
            setpz2(zq+s*3, subpz2(s04_mj_s26, w8_s15_mj_s37));
            setpz2(zq+s*4, subpz2(a04_p1_a26,    a15_p1_a37));
            setpz2(zq+s*5, subpz2(s04_pj_s26, v8_s15_pj_s37));
            setpz2(zq+s*6, subpz2(a04_m1_a26,  j_a15_m1_a37));
            setpz2(zq+s*7, addpz2(s04_mj_s26, w8_s15_mj_s37));
        }
    }
    else {
        const Vec2d x0 = getpz(x[0]);
        const Vec2d x1 = getpz(x[1]);
        const Vec2d x2 = getpz(x[2]);
        const Vec2d x3 = getpz(x[3]);
        const Vec2d x4 = getpz(x[4]);
        const Vec2d x5 = getpz(x[5]);
        const Vec2d x6 = getpz(x[6]);
        const Vec2d x7 = getpz(x[7]);
        const Vec2d  a04 =      addpz(x0, x4);
        const Vec2d  s04 =      subpz(x0, x4);
        const Vec2d  a26 =      addpz(x2, x6);
        const Vec2d js26 = jxpz(subpz(x2, x6));
        const Vec2d  a15 =      addpz(x1, x5);
        const Vec2d  s15 =      subpz(x1, x5);
        const Vec2d  a37 =      addpz(x3, x7);
        const Vec2d js37 = jxpz(subpz(x3, x7));
        const Vec2d    a04_p1_a26 =       addpz(a04,  a26);
        const Vec2d    s04_pj_s26 =       addpz(s04, js26);
        const Vec2d    a04_m1_a26 =       subpz(a04,  a26);
        const Vec2d    s04_mj_s26 =       subpz(s04, js26);
        const Vec2d    a15_p1_a37 =       addpz(a15,  a37);
        const Vec2d v8_s15_pj_s37 = v8xpz(addpz(s15, js37));
        const Vec2d  j_a15_m1_a37 =  jxpz(subpz(a15,  a37));
        const Vec2d w8_s15_mj_s37 = w8xpz(subpz(s15, js37));
        setpz(z[0], addpz(a04_p1_a26,    a15_p1_a37));
        setpz(z[1], addpz(s04_pj_s26, v8_s15_pj_s37));
        setpz(z[2], addpz(a04_m1_a26,  j_a15_m1_a37));
        setpz(z[3], subpz(s04_mj_s26, w8_s15_mj_s37));
        setpz(z[4], subpz(a04_p1_a26,    a15_p1_a37));
        setpz(z[5], subpz(s04_pj_s26, v8_s15_pj_s37));
        setpz(z[6], subpz(a04_m1_a26,  j_a15_m1_a37));
        setpz(z[7], addpz(s04_mj_s26, w8_s15_mj_s37));
    }
}

void invcore8(const int n, const int s,
        complex_vector x, complex_vector y, const_complex_vector W) noexcept
{
    const int m  = n/8;
    const int N  = n*s;
    const int N0 = 0;
    const int N1 = N/8;
    const int N2 = N1*2;
    const int N3 = N1*3;
    const int N4 = N1*4;
    const int N5 = N1*5;
    const int N6 = N1*6;
    const int N7 = N1*7;
    if (s >= 2) {
        for (int p = 0; p < m; p++) {
            const int sp = s*p;
            const int s8p = 8*sp;
            const Vec4d w1p = duppz3(W[N-1*sp]);
            const Vec4d w2p = duppz3(W[N-2*sp]);
            const Vec4d w3p = duppz3(W[N-3*sp]);
            const Vec4d w4p = mulpz2(w2p,w2p);
            const Vec4d w5p = mulpz2(w2p,w3p);
            const Vec4d w6p = mulpz2(w3p,w3p);
            const Vec4d w7p = mulpz2(w3p,w4p);
            for (int q = 0; q < s; q += 2) {
                complex_vector xq_sp  = x + q + sp;
                complex_vector yq_s8p = y + q + s8p;
                const Vec4d x0 = getpz2(xq_sp+N0);
                const Vec4d x1 = getpz2(xq_sp+N1);
                const Vec4d x2 = getpz2(xq_sp+N2);
                const Vec4d x3 = getpz2(xq_sp+N3);
                const Vec4d x4 = getpz2(xq_sp+N4);
                const Vec4d x5 = getpz2(xq_sp+N5);
                const Vec4d x6 = getpz2(xq_sp+N6);
                const Vec4d x7 = getpz2(xq_sp+N7);
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
                setpz2(yq_s8p+s*0,             addpz2(a04_p1_a26,    a15_p1_a37));
                setpz2(yq_s8p+s*1, mulpz2(w1p, addpz2(s04_pj_s26, v8_s15_pj_s37)));
                setpz2(yq_s8p+s*2, mulpz2(w2p, addpz2(a04_m1_a26,  j_a15_m1_a37)));
                setpz2(yq_s8p+s*3, mulpz2(w3p, subpz2(s04_mj_s26, w8_s15_mj_s37)));
                setpz2(yq_s8p+s*4, mulpz2(w4p, subpz2(a04_p1_a26,    a15_p1_a37)));
                setpz2(yq_s8p+s*5, mulpz2(w5p, subpz2(s04_pj_s26, v8_s15_pj_s37)));
                setpz2(yq_s8p+s*6, mulpz2(w6p, subpz2(a04_m1_a26,  j_a15_m1_a37)));
                setpz2(yq_s8p+s*7, mulpz2(w7p, addpz2(s04_mj_s26, w8_s15_mj_s37)));
            }
        }
    }
    else {
        for (int p = 0; p < m; p++) {
            complex_vector x_p  = x + p;
            complex_vector y_8p = y + 8*p;
            const Vec2d w1p = cnjpz(getpz(W[p]));
            const Vec2d w2p = mulpz(w1p,w1p);
            const Vec2d w3p = mulpz(w1p,w2p);
            const Vec2d w4p = mulpz(w2p,w2p);
            const Vec2d w5p = mulpz(w2p,w3p);
            const Vec2d w6p = mulpz(w3p,w3p);
            const Vec2d w7p = mulpz(w3p,w4p);
            const Vec2d x0 = getpz(x_p[N0]);
            const Vec2d x1 = getpz(x_p[N1]);
            const Vec2d x2 = getpz(x_p[N2]);
            const Vec2d x3 = getpz(x_p[N3]);
            const Vec2d x4 = getpz(x_p[N4]);
            const Vec2d x5 = getpz(x_p[N5]);
            const Vec2d x6 = getpz(x_p[N6]);
            const Vec2d x7 = getpz(x_p[N7]);
            const Vec2d  a04 =      addpz(x0, x4);
            const Vec2d  s04 =      subpz(x0, x4);
            const Vec2d  a26 =      addpz(x2, x6);
            const Vec2d js26 = jxpz(subpz(x2, x6));
            const Vec2d  a15 =      addpz(x1, x5);
            const Vec2d  s15 =      subpz(x1, x5);
            const Vec2d  a37 =      addpz(x3, x7);
            const Vec2d js37 = jxpz(subpz(x3, x7));
            const Vec2d    a04_p1_a26 =       addpz(a04,  a26);
            const Vec2d    s04_pj_s26 =       addpz(s04, js26);
            const Vec2d    a04_m1_a26 =       subpz(a04,  a26);
            const Vec2d    s04_mj_s26 =       subpz(s04, js26);
            const Vec2d    a15_p1_a37 =       addpz(a15,  a37);
            const Vec2d v8_s15_pj_s37 = v8xpz(addpz(s15, js37));
            const Vec2d  j_a15_m1_a37 =  jxpz(subpz(a15,  a37));
            const Vec2d w8_s15_mj_s37 = w8xpz(subpz(s15, js37));
            setpz(y_8p[0],            addpz(a04_p1_a26,    a15_p1_a37));
            setpz(y_8p[1], mulpz(w1p, addpz(s04_pj_s26, v8_s15_pj_s37)));
            setpz(y_8p[2], mulpz(w2p, addpz(a04_m1_a26,  j_a15_m1_a37)));
            setpz(y_8p[3], mulpz(w3p, subpz(s04_mj_s26, w8_s15_mj_s37)));
            setpz(y_8p[4], mulpz(w4p, subpz(a04_p1_a26,    a15_p1_a37)));
            setpz(y_8p[5], mulpz(w5p, subpz(s04_pj_s26, v8_s15_pj_s37)));
            setpz(y_8p[6], mulpz(w6p, subpz(a04_m1_a26,  j_a15_m1_a37)));
            setpz(y_8p[7], mulpz(w7p, addpz(s04_mj_s26, w8_s15_mj_s37)));
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

void invend5(const int s, const bool eo,
        complex_vector x, complex_vector y, const_complex_vector W) noexcept
{
    const cpx w1 = getpz(W[4*s]);
    const cpx w2 = mulpz(w1,w1);
    const cpx w3 = mulpz(w1,w2);
    const cpx w4 = mulpz(w2,w2);
    complex_vector z = eo ? y : x;
    for (int q = 0; q < s; q++) {
        const cpx a = getpz(x[q+s*0]);
        const cpx b = getpz(x[q+s*1]);
        const cpx c = getpz(x[q+s*2]);
        const cpx d = getpz(x[q+s*3]);
        const cpx e = getpz(x[q+s*4]);
        setpz(z[q+s*0], a + b + c + d + e);
        setpz(z[q+s*1], a + w1*b + w2*c + w3*d + w4*e);
        setpz(z[q+s*2], a + w2*b + w4*c + w1*d + w3*e);
        setpz(z[q+s*3], a + w3*b + w1*c + w4*d + w2*e);
        setpz(z[q+s*4], a + w4*b + w3*c + w2*d + w1*e);
    }
}

void invcore5(const int n, const int s,
        complex_vector x, complex_vector y, const_complex_vector W) noexcept
{
    const int N  = n*s;
    const int m  = n/5;
    const int N0 = 0;
    const int N1 = N/5;
    const int N2 = N1*2;
    const int N3 = N1*3;
    const int N4 = N1*4;
    const cpx w1 = getpz(W[N4]);
    const cpx w2 = mulpz(w1,w1);
    const cpx w3 = mulpz(w1,w2);
    const cpx w4 = mulpz(w2,w2);
    for (int p = 0; p < m; p++) {
        const int sp = s*p;
        const cpx w1p = getpz(W[N-sp]);
        const cpx w2p = mulpz(w1p,w1p);
        const cpx w3p = mulpz(w1p,w2p);
        const cpx w4p = mulpz(w2p,w2p);
        for (int q = 0; q < s; q++) {
            const int q_sp = q + sp;
            const cpx a = getpz(x[q_sp+N0]);
            const cpx b = getpz(x[q_sp+N1]);
            const cpx c = getpz(x[q_sp+N2]);
            const cpx d = getpz(x[q_sp+N3]);
            const cpx e = getpz(x[q_sp+N4]);
            const int q_s5p = q + sp*5;
            setpz(y[q_s5p+s*0],  a + b + c + d + e);
            setpz(y[q_s5p+s*1], (a + w1*b + w2*c + w3*d + w4*e)*w1p);
            setpz(y[q_s5p+s*2], (a + w2*b + w4*c + w1*d + w3*e)*w2p);
            setpz(y[q_s5p+s*3], (a + w3*b + w1*c + w4*d + w2*e)*w3p);
            setpz(y[q_s5p+s*4], (a + w4*b + w3*c + w2*d + w1*e)*w4p);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

void invend3(const int s, const bool eo,
        complex_vector x, complex_vector y, const_complex_vector W) noexcept
{
    const cpx w1 = getpz(W[2*s]);
    const cpx w2 = mulpz(w1,w1);
    complex_vector z = eo ? y : x;
    for (int q = 0; q < s; q++) {
        const cpx a = getpz(x[q+s*0]);
        const cpx b = getpz(x[q+s*1]);
        const cpx c = getpz(x[q+s*2]);
        setpz(z[q+s*0], a + b + c);
        setpz(z[q+s*1], a + w1*b + w2*c);
        setpz(z[q+s*2], a + w2*b + w1*c);
    }
}

void invcore3(const int n, const int s,
        complex_vector x, complex_vector y, const_complex_vector W) noexcept
{
    const int N  = n*s;
    const int m  = n/3;
    const int N0 = 0;
    const int N1 = N/3;
    const int N2 = N1*2;
    const cpx w1 = getpz(W[N2]);
    const cpx w2 = mulpz(w1,w1);
    for (int p = 0; p < m; p++) {
        const int sp = s*p;
        const cpx w1p = getpz(W[N-sp]);
        const cpx w2p = mulpz(w1p,w1p);
        for (int q = 0; q < s; q++) {
            const int q_sp = q + sp;
            const cpx a = getpz(x[q_sp+N0]);
            const cpx b = getpz(x[q_sp+N1]);
            const cpx c = getpz(x[q_sp+N2]);
            const int q_s3p = q + sp*3;
            setpz(y[q_s3p+s*0],  a + b + c);
            setpz(y[q_s3p+s*1], (a + w1*b + w2*c)*w1p);
            setpz(y[q_s3p+s*2], (a + w2*b + w1*c)*w2p);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
// Any Size IFFT except Radix-2,3,5
///////////////////////////////////////////////////////////////////////////////

void invfftany(const int r, const int n, const int s, const bool eo,
        complex_vector x, complex_vector y, const_complex_vector W) noexcept
{
    const Vec2d zero { 0, 0 };
    const int N = n*s;
    int k = r;
    while (n%k != 0) {
        if (k*k > n) { k = n; break; }
        k += 2;
    }
    if (k == n) {
        for (int q = 0; q < s; q++) {
            for (int i = 0; i < k; i++) {
                cpx z = zero;
                for (int j = 0; j < k; j++) {
                    const cpx a   = getpz(x[q+s*j]);
                    const cpx wij = getpz(W[N-s*((i*j)%k)]);
                    z = z + a*wij;
                }
                setpz(y[q+s*i], z);
            }
        }
        if (!eo) for (int p = 0; p < N; p++) setpz(x[p], getpz(y[p]));
    }
    else {
        const int m  = n/k;
        const int ms = m*s;
        for (int p = 0; p < m; p++) {
            const int sp = s*p;
            for (int q = 0; q < s; q++) {
                const int q_sp  = q + sp;
                const int q_spk = q + sp*k;
                for (int i = 0; i < k; i++) {
                    cpx z = zero;
                    for (int j = 0; j < k; j++) {
                        const cpx a   = getpz(x[q_sp+ms*j]);
                        const cpx wij = getpz(W[N-ms*((i*j)%k)]);
                        z = z + a*wij;
                    }
                    const cpx wip = getpz(W[N-i*sp]);
                    setpz(y[q_spk+s*i], z * wip);
                }
            }
        }
        invfftany(k, m, k*s, !eo, y, x, W);
    }
}

///////////////////////////////////////////////////////////////////////////////
// Mixed Radix IFFT
///////////////////////////////////////////////////////////////////////////////

void invfft(const int n, const int s, const bool eo,
        complex_vector x, complex_vector y, const_complex_vector W) noexcept
{
    const int N = n*s;
    if (N < 2) return;
    if (n%8 == 0) {
        if (n == 8)
            invend8(s, eo, x, y);
        else {
            invcore8(n, s, x, y, W);
            invfft(n/8, 8*s, !eo, y, x, W);
        }
    }
    else if (n%4 == 0) {
        if (n == 4)
            invend4(s, eo, x, y);
        else {
            invcore4(n, s, x, y, W);
            invfft(n/4, 4*s, !eo, y, x, W);
        }
    }
    else if (n%2 == 0) {
        if (n == 2)
            invend2(s, eo, x, y);
        else {
            invcore2(n, s, x, y, W);
            invfft(n/2, 2*s, !eo, y, x, W);
        }
    }
    else if (n%5 == 0) {
        if (n == 5)
            invend5(s, eo, x, y, W);
        else {
            invcore5(n, s, x, y, W);
            invfft(n/5, 5*s, !eo, y, x, W);
        }
    }
    else if (n%3 == 0) {
        if (n == 3)
            invend3(s, eo, x, y, W);
        else {
            invcore3(n, s, x, y, W);
            invfft(n/3, 3*s, !eo, y, x, W);
        }
    }
    else invfftany(7, n, s, eo, x, y, W);
}

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

    void setup(const int n)
    {
        const double theta0 = 2*M_PI/n;
        N = n;
        weight.setup(n+1); W = weight.get();

        for (int p = 0; p <= n; p++) {
            W[p] = complex_t(cos(p*theta0), -sin(p*theta0));
        }
    }

    inline void setup2(const int n) { setup(1 << n); }

    ///////////////////////////////////////////////////////////////////////////

    void fwd(complex_vector x, complex_vector y) const noexcept
    {
        const Vec2d rN = cmplx(1.0/N, 1.0/N);

        fwdfft(N, 1, 0, x, y, W);
        for (int k = 0; k < N; k++) setpz(x[k], mulpd(rN, getpz(x[k])));
    }

    void fwd0(complex_vector x, complex_vector y) const noexcept
    {
        fwdfft(N, 1, 0, x, y, W);
    }

    void fwdu(complex_vector x, complex_vector y) const noexcept
    {
        const double ssrN = sqrt(1.0/N);
        const Vec2d srN = cmplx(ssrN, ssrN);

        fwdfft(N, 1, 0, x, y, W);
        for (int k = 0; k < N; k++) setpz(x[k], mulpd(srN, getpz(x[k])));
    }

    inline void fwdn(complex_vector x, complex_vector y) const noexcept
    {
        fwd(x, y);
    }

    ///////////////////////////////////////////////////////////////////////////

    void inv(complex_vector x, complex_vector y) const noexcept
    {
        invfft(N, 1, 0, x, y, W);
    }

    inline void inv0(complex_vector x, complex_vector y) const noexcept
    {
        inv(x, y);
    }

    void invu(complex_vector x, complex_vector y) const noexcept
    {
        const double ssrN = sqrt(1.0/N);
        const Vec2d srN = cmplx(ssrN, ssrN);

        invfft(N, 1, 0, x, y, W);
        for (int p = 0; p < N; p++) setpz(x[p], mulpd(srN, getpz(x[p])));
    }

    void invn(complex_vector x, complex_vector y) const noexcept
    {
        const Vec2d rN = cmplx(1.0/N, 1.0/N);

        invfft(N, 1, 0, x, y, W);
        for (int p = 0; p < N; p++) setpz(x[p], mulpd(rN, getpz(x[p])));
    }
};

} /////////////////////////////////////////////////////////////////////////////

#endif // otfft_mixedradix_h
