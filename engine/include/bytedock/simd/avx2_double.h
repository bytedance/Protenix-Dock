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

#include "immintrin.h"

namespace bytedock {

/**
 * Reference: https://gitlab.com/gromacs/gromacs/-/blob/release-2024/
 *            src/gromacs/simd/include/gromacs/simd
 */

#define BDOCK_SIMD_DOUBLE_WIDTH 4

struct simd_double {
    simd_double() {}
    simd_double(double v) : internal_(_mm256_set1_pd(v)) {}
    simd_double(__m256d v) : internal_(v) {}

    __m256d internal_;
};

static inline simd_double simd_load(const double* src) {
    // assert(std::size_t(src) % 32 == 0);
    return { _mm256_load_pd(src) };
}

static inline simd_double simd_loadu(const double* src) {
    return { _mm256_loadu_pd(src) };
}

static inline void simd_store(double* dst, simd_double& src) {
    // assert(std::size_t(dst) % 32 == 0);
    _mm256_store_pd(dst, src.internal_);
}

static inline void simd_storeu(double* dst, simd_double& src) {
    _mm256_storeu_pd(dst, src.internal_);
}

static inline simd_double operator+(simd_double a, simd_double b) {
    return { _mm256_add_pd(a.internal_, b.internal_) };
}

static inline simd_double operator-(simd_double a, simd_double b) {
    return { _mm256_sub_pd(a.internal_, b.internal_) };
}

static inline simd_double operator*(simd_double a, simd_double b) {
    return { _mm256_mul_pd(a.internal_, b.internal_) };
}

// a*b+c
static inline simd_double simd_fma(simd_double a, simd_double b, simd_double c) {
    return { _mm256_fmadd_pd(a.internal_, b.internal_, c.internal_) };
}

// a*b-c
static inline simd_double simd_fms(simd_double a, simd_double b, simd_double c) {
    return { _mm256_fmsub_pd(a.internal_, b.internal_, c.internal_) };
}

// native 1/sqrt(x)
static inline simd_double simd_rsqrt(simd_double x) {
    return { _mm256_cvtps_pd(_mm_rsqrt_ps(_mm256_cvtpd_ps(x.internal_))) };
}

static inline simd_double simd_rsqrt_iter(simd_double lu, simd_double x) {
    simd_double tmp1 = x * lu;
    simd_double tmp2 = simd_double(-0.5) * lu;
    tmp1 = simd_fma(tmp1, lu, simd_double(-3.0));
    return tmp1 * tmp2;
}

// 1/sqrt(x) with higher precision
static inline simd_double simd_invsqrt(simd_double x) {
    simd_double lu = simd_rsqrt(x);
    lu = simd_rsqrt_iter(lu, x);
    lu = simd_rsqrt_iter(lu, x);
    return lu;
}

static inline double simd_reduce(simd_double a) {
    __m128d a0, a1;
    a.internal_ = _mm256_add_pd(a.internal_, _mm256_permute_pd(a.internal_, 0b0101));
    a0 = _mm256_castpd256_pd128(a.internal_);
    a1 = _mm256_extractf128_pd(a.internal_, 0x1);
    a0 = _mm_add_sd(a0, a1);
    return *reinterpret_cast<double*>(&a0);
}

// xyz, xyz, ... => xx..., yy..., zz...
static inline void simd_loadu_xyz(const double* base, const std::int32_t* offset,
                                  simd_double& rx, simd_double& ry, simd_double& rz) {
    __m256d t1 = _mm256_loadu_pd(base + 3*offset[0]);
    __m256d t2 = _mm256_loadu_pd(base + 3*offset[1]);
    __m256d t3 = _mm256_loadu_pd(base + 3*offset[2]);
    __m256d t4 = _mm256_loadu_pd(base + 3*offset[3]);
    __m256d t5 = _mm256_unpacklo_pd(t1, t2);
    __m256d t6 = _mm256_unpackhi_pd(t1, t2);
    __m256d t7 = _mm256_unpacklo_pd(t3, t4);
    __m256d t8 = _mm256_unpackhi_pd(t3, t4);
    rx.internal_ = _mm256_permute2f128_pd(t5, t7, 0x20);
    ry.internal_ = _mm256_permute2f128_pd(t6, t8, 0x20);
    rz.internal_ = _mm256_permute2f128_pd(t5, t7, 0x31);
}

// xx..., yy..., zz... => +xyz, +xyz, ...
static inline void simd_incru_xyz(double* base, const std::int32_t* offset,
                                  simd_double rx, simd_double ry, simd_double rz) {
    __m256d t0 = _mm256_unpacklo_pd(rx.internal_, ry.internal_);
    __m256d t1 = _mm256_unpackhi_pd(rx.internal_, ry.internal_);
    __m128d t2 = _mm256_extractf128_pd(rz.internal_, 0x1);

    __m128d tA = _mm_loadu_pd(base + 3*offset[0]);
    __m128d tB = _mm_load_sd(base + 3*offset[0] + 2);
    tA = _mm_add_pd(tA, _mm256_castpd256_pd128(t0));
    tB = _mm_add_pd(tB, _mm256_castpd256_pd128(rz.internal_));
    _mm_storeu_pd(base + 3*offset[0], tA);
    _mm_store_sd(base + 3*offset[0] + 2, tB);

    tA = _mm_loadu_pd(base + 3*offset[1]);
    tB = _mm_loadh_pd(_mm_setzero_pd(), base + 3*offset[1] + 2);
    tA = _mm_add_pd(tA, _mm256_castpd256_pd128(t1));
    tB = _mm_add_pd(tB, _mm256_castpd256_pd128(rz.internal_));
    _mm_storeu_pd(base + 3*offset[1], tA);
    _mm_storeh_pd(base + 3*offset[1] + 2, tB);

    tA = _mm_loadu_pd(base + 3*offset[2]);
    tB = _mm_load_sd(base + 3*offset[2] + 2);
    tA = _mm_add_pd(tA, _mm256_extractf128_pd(t0, 0x1));
    tB = _mm_add_pd(tB, t2);
    _mm_storeu_pd(base + 3*offset[2], tA);
    _mm_store_sd(base + 3*offset[2] + 2, tB);

    tA = _mm_loadu_pd(base + 3*offset[3]);
    tB = _mm_loadh_pd(_mm_setzero_pd(), base + 3*offset[3] + 2);
    tA = _mm_add_pd(tA, _mm256_extractf128_pd(t1, 0x1));
    tB = _mm_add_pd(tB, t2);
    _mm_storeu_pd(base + 3*offset[3], tA);
    _mm_storeh_pd(base + 3*offset[3] + 2, tB);
}

// xx..., yy..., zz... => -xyz, -xyz, ...
static inline void simd_decru_xyz(double* base, const std::int32_t* offset,
                                  simd_double rx, simd_double ry, simd_double rz) {
    __m256d t0 = _mm256_unpacklo_pd(rx.internal_, ry.internal_);
    __m256d t1 = _mm256_unpackhi_pd(rx.internal_, ry.internal_);
    __m128d t2 = _mm256_extractf128_pd(rz.internal_, 0x1);

    __m128d tA = _mm_loadu_pd(base + 3*offset[0]);
    __m128d tB = _mm_load_sd(base + 3*offset[0] + 2);
    tA = _mm_sub_pd(tA, _mm256_castpd256_pd128(t0));
    tB = _mm_sub_pd(tB, _mm256_castpd256_pd128(rz.internal_));
    _mm_storeu_pd(base + 3*offset[0], tA);
    _mm_store_sd(base + 3*offset[0] + 2, tB);

    tA = _mm_loadu_pd(base + 3*offset[1]);
    tB = _mm_loadh_pd(_mm_setzero_pd(), base + 3*offset[1] + 2);
    tA = _mm_sub_pd(tA, _mm256_castpd256_pd128(t1));
    tB = _mm_sub_pd(tB, _mm256_castpd256_pd128(rz.internal_));
    _mm_storeu_pd(base + 3*offset[1], tA);
    _mm_storeh_pd(base + 3*offset[1] + 2, tB);

    tA = _mm_loadu_pd(base + 3*offset[2]);
    tB = _mm_load_sd(base + 3*offset[2] + 2);
    tA = _mm_sub_pd(tA, _mm256_extractf128_pd(t0, 0x1));
    tB = _mm_sub_pd(tB, t2);
    _mm_storeu_pd(base + 3*offset[2], tA);
    _mm_store_sd(base + 3*offset[2] + 2, tB);

    tA = _mm_loadu_pd(base + 3*offset[3]);
    tB = _mm_loadh_pd(_mm_setzero_pd(), base + 3*offset[3] + 2);
    tA = _mm_sub_pd(tA, _mm256_extractf128_pd(t1, 0x1));
    tB = _mm_sub_pd(tB, t2);
    _mm_storeu_pd(base + 3*offset[3], tA);
    _mm_storeh_pd(base + 3*offset[3] + 2, tB);
}

}
