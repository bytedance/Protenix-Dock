/*
* Copyright (C) 2025 ByteDance and/or its affiliates
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#pragma once

#include <cstdint>
#include <arm_neon.h>

namespace bytedock {

/**
 * ARM NEON SIMD for double precision (128-bit, 2 doubles per register).
 * Mirrors the AVX2 API in avx2_double.h so that pair.h code works unchanged.
 */

#define BDOCK_SIMD_DOUBLE_WIDTH 2

struct simd_double {
    simd_double() {}
    simd_double(double v) : internal_(vdupq_n_f64(v)) {}
    simd_double(float64x2_t v) : internal_(v) {}

    float64x2_t internal_;
};

static inline simd_double simd_load(const double* src) {
    return { vld1q_f64(src) };
}

static inline simd_double simd_loadu(const double* src) {
    return { vld1q_f64(src) };
}

static inline void simd_store(double* dst, simd_double& src) {
    vst1q_f64(dst, src.internal_);
}

static inline void simd_storeu(double* dst, simd_double& src) {
    vst1q_f64(dst, src.internal_);
}

static inline simd_double operator+(simd_double a, simd_double b) {
    return { vaddq_f64(a.internal_, b.internal_) };
}

static inline simd_double operator-(simd_double a, simd_double b) {
    return { vsubq_f64(a.internal_, b.internal_) };
}

static inline simd_double operator*(simd_double a, simd_double b) {
    return { vmulq_f64(a.internal_, b.internal_) };
}

static inline simd_double simd_max(simd_double a, simd_double b) {
    return { vmaxq_f64(a.internal_, b.internal_) };
}

static inline simd_double simd_min(simd_double a, simd_double b) {
    return { vminq_f64(a.internal_, b.internal_) };
}

// a*b+c
static inline simd_double simd_fma(simd_double a, simd_double b, simd_double c) {
    return { vfmaq_f64(c.internal_, a.internal_, b.internal_) };
}

// a*b-c
static inline simd_double simd_fms(simd_double a, simd_double b, simd_double c) {
    return { vfmaq_f64(vnegq_f64(c.internal_), a.internal_, b.internal_) };
}

// Approximate 1/x with Newton-Raphson refinement
static inline simd_double simd_rcp(simd_double x) {
    float64x2_t est = vrecpeq_f64(x.internal_);
    // 2 Newton iterations: est *= vrecpsq_f64(x, est) computes (2 - x*est)
    est = vmulq_f64(est, vrecpsq_f64(x.internal_, est));
    est = vmulq_f64(est, vrecpsq_f64(x.internal_, est));
    return { est };
}

// Approximate 1/sqrt(x)
static inline simd_double simd_rsqrt(simd_double x) {
    return { vrsqrteq_f64(x.internal_) };
}

// Newton-Raphson refinement for rsqrt using NEON's dedicated instruction
static inline simd_double simd_rsqrt_iter(simd_double lu, simd_double x) {
    // vrsqrtsq_f64(a, b) computes (3 - a*b) / 2
    float64x2_t step = vrsqrtsq_f64(vmulq_f64(x.internal_, lu.internal_), lu.internal_);
    return { vmulq_f64(lu.internal_, step) };
}

// 1/sqrt(x) with higher precision (2 Newton iterations)
static inline simd_double simd_invsqrt(simd_double x) {
    simd_double lu = simd_rsqrt(x);
    lu = simd_rsqrt_iter(lu, x);
    lu = simd_rsqrt_iter(lu, x);
    return lu;
}

static inline double simd_reduce(simd_double a) {
    return vaddvq_f64(a.internal_);
}

// Gather xyz: AoS [x0,y0,z0], [x1,y1,z1] => SoA [x0,x1], [y0,y1], [z0,z1]
static inline void simd_loadu_xyz(const double* base, const std::int32_t* offset,
                                  simd_double& rx, simd_double& ry, simd_double& rz) {
    // Load two consecutive pairs: [x0, y0] and [x1, y1]
    float64x2_t xy0 = vld1q_f64(base + 3 * offset[0]);
    float64x2_t xy1 = vld1q_f64(base + 3 * offset[1]);

    // Unzip to separate x and y channels
    rx.internal_ = vuzp1q_f64(xy0, xy1);  // [x0, x1]
    ry.internal_ = vuzp2q_f64(xy0, xy1);  // [y0, y1]

    // Z values loaded individually
    rz.internal_ = vsetq_lane_f64(
        base[3 * offset[1] + 2],
        vdupq_n_f64(base[3 * offset[0] + 2]),
        1
    );
}

// Scatter add: SoA [gx0,gx1], [gy0,gy1], [gz0,gz1] => AoS base[offset[i]] += {gx,gy,gz}
static inline void simd_incru_xyz(double* base, const std::int32_t* offset,
                                  simd_double rx, simd_double ry, simd_double rz) {
    base[3*offset[0]]   += vgetq_lane_f64(rx.internal_, 0);
    base[3*offset[0]+1] += vgetq_lane_f64(ry.internal_, 0);
    base[3*offset[0]+2] += vgetq_lane_f64(rz.internal_, 0);
    base[3*offset[1]]   += vgetq_lane_f64(rx.internal_, 1);
    base[3*offset[1]+1] += vgetq_lane_f64(ry.internal_, 1);
    base[3*offset[1]+2] += vgetq_lane_f64(rz.internal_, 1);
}

// Scatter sub: SoA [gx0,gx1], [gy0,gy1], [gz0,gz1] => AoS base[offset[i]] -= {gx,gy,gz}
static inline void simd_decru_xyz(double* base, const std::int32_t* offset,
                                  simd_double rx, simd_double ry, simd_double rz) {
    base[3*offset[0]]   -= vgetq_lane_f64(rx.internal_, 0);
    base[3*offset[0]+1] -= vgetq_lane_f64(ry.internal_, 0);
    base[3*offset[0]+2] -= vgetq_lane_f64(rz.internal_, 0);
    base[3*offset[1]]   -= vgetq_lane_f64(rx.internal_, 1);
    base[3*offset[1]+1] -= vgetq_lane_f64(ry.internal_, 1);
    base[3*offset[1]+2] -= vgetq_lane_f64(rz.internal_, 1);
}

}
