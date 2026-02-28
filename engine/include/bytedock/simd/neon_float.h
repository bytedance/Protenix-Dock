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

#define BDOCK_SIMD_FLOAT_WIDTH 4

struct simd_float {
    simd_float() {}
    simd_float(float v) : internal_(vdupq_n_f32(v)) {}
    simd_float(float32x4_t v) : internal_(v) {}

    float32x4_t internal_;
};

static inline simd_float simd_load(const float* src) {
    return { vld1q_f32(src) };
}

static inline simd_float simd_loadu(const float* src) {
    return { vld1q_f32(src) };
}

static inline void simd_store(float* dst, simd_float& src) {
    vst1q_f32(dst, src.internal_);
}

static inline void simd_storeu(float* dst, simd_float& src) {
    vst1q_f32(dst, src.internal_);
}

static inline simd_float operator+(simd_float a, simd_float b) {
    return { vaddq_f32(a.internal_, b.internal_) };
}

static inline simd_float operator-(simd_float a, simd_float b) {
    return { vsubq_f32(a.internal_, b.internal_) };
}

static inline simd_float operator*(simd_float a, simd_float b) {
    return { vmulq_f32(a.internal_, b.internal_) };
}

static inline simd_float simd_max(simd_float a, simd_float b) {
    return { vmaxq_f32(a.internal_, b.internal_) };
}

static inline simd_float simd_min(simd_float a, simd_float b) {
    return { vminq_f32(a.internal_, b.internal_) };
}

// a*b+c
static inline simd_float simd_fma(simd_float a, simd_float b, simd_float c) {
    return { vfmaq_f32(c.internal_, a.internal_, b.internal_) };
}

// a*b-c
static inline simd_float simd_fms(simd_float a, simd_float b, simd_float c) {
    return { vfmaq_f32(vnegq_f32(c.internal_), a.internal_, b.internal_) };
}

// Approximate 1/x (two Newton-Raphson steps for float precision)
static inline simd_float simd_rcp(simd_float x) {
    float32x4_t est = vrecpeq_f32(x.internal_);
    est = vmulq_f32(est, vrecpsq_f32(x.internal_, est));
    est = vmulq_f32(est, vrecpsq_f32(x.internal_, est));
    return { est };
}

// native 1/sqrt(x)
static inline simd_float simd_rsqrt(simd_float x) {
    return { vrsqrteq_f32(x.internal_) };
}

static inline simd_float simd_rsqrt_iter(simd_float lu, simd_float x) {
    // vrsqrtsq_f32(x*lu, lu) = (3 - (x*lu)*lu) / 2
    // result = lu * (3 - x*lu^2) / 2
    return { vmulq_f32(lu.internal_,
             vrsqrtsq_f32(vmulq_f32(x.internal_, lu.internal_), lu.internal_)) };
}

// 1/sqrt(x) with higher precision (two Newton-Raphson steps)
static inline simd_float simd_invsqrt(simd_float x) {
    simd_float lu = simd_rsqrt(x);
    lu = simd_rsqrt_iter(lu, x);
    lu = simd_rsqrt_iter(lu, x);
    return lu;
}

static inline float simd_reduce(simd_float a) {
    float32x2_t sum = vadd_f32(vget_low_f32(a.internal_), vget_high_f32(a.internal_));
    sum = vpadd_f32(sum, sum);
    return vget_lane_f32(sum, 0);
}

// xyz, xyz, ... => xx..., yy..., zz...  (4-wide gather)
static inline void simd_loadu_xyz(const float* base, const std::int32_t* offset,
                                  simd_float& rx, simd_float& ry, simd_float& rz) {
    float32x4_t t1 = vld1q_f32(base + 3*offset[0]);  // x0 y0 z0 ?
    float32x4_t t2 = vld1q_f32(base + 3*offset[1]);  // x1 y1 z1 ?
    float32x4_t t3 = vld1q_f32(base + 3*offset[2]);  // x2 y2 z2 ?
    float32x4_t t4 = vld1q_f32(base + 3*offset[3]);  // x3 y3 z3 ?

    // vtrnq_f32(a, b):
    //   val[0] = {a[0], b[0], a[2], b[2]}
    //   val[1] = {a[1], b[1], a[3], b[3]}
    float32x4x2_t trn12 = vtrnq_f32(t1, t2);
    // trn12.val[0] = {x0, x1, z0, z1}
    // trn12.val[1] = {y0, y1, ?,  ?}
    float32x4x2_t trn34 = vtrnq_f32(t3, t4);
    // trn34.val[0] = {x2, x3, z2, z3}
    // trn34.val[1] = {y2, y3, ?,  ?}

    rx.internal_ = vcombine_f32(vget_low_f32(trn12.val[0]),
                                vget_low_f32(trn34.val[0]));   // {x0,x1,x2,x3}
    ry.internal_ = vcombine_f32(vget_low_f32(trn12.val[1]),
                                vget_low_f32(trn34.val[1]));   // {y0,y1,y2,y3}
    rz.internal_ = vcombine_f32(vget_high_f32(trn12.val[0]),
                                vget_high_f32(trn34.val[0]));  // {z0,z1,z2,z3}
}

// xx..., yy..., zz... => +xyz, +xyz, ...  (4-wide scatter add)
static inline void simd_incru_xyz(float* base, const std::int32_t* offset,
                                  simd_float rx, simd_float ry, simd_float rz) {
    alignas(16) float gx[4], gy[4], gz[4];
    vst1q_f32(gx, rx.internal_);
    vst1q_f32(gy, ry.internal_);
    vst1q_f32(gz, rz.internal_);
    for (int i = 0; i < 4; ++i) {
        base[3*offset[i] + 0] += gx[i];
        base[3*offset[i] + 1] += gy[i];
        base[3*offset[i] + 2] += gz[i];
    }
}

// xx..., yy..., zz... => -xyz, -xyz, ...  (4-wide scatter subtract)
static inline void simd_decru_xyz(float* base, const std::int32_t* offset,
                                  simd_float rx, simd_float ry, simd_float rz) {
    alignas(16) float gx[4], gy[4], gz[4];
    vst1q_f32(gx, rx.internal_);
    vst1q_f32(gy, ry.internal_);
    vst1q_f32(gz, rz.internal_);
    for (int i = 0; i < 4; ++i) {
        base[3*offset[i] + 0] -= gx[i];
        base[3*offset[i] + 1] -= gy[i];
        base[3*offset[i] + 2] -= gz[i];
    }
}

}
