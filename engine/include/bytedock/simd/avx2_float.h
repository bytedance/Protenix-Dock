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

#define BDOCK_SIMD_FLOAT_WIDTH 8

struct simd_float {
    simd_float() {}
    simd_float(float v) : internal_(_mm256_set1_ps(v)) {}
    simd_float(__m256 v) : internal_(v) {}

    __m256 internal_;
};

static inline simd_float simd_load(const float* src) {
    // assert(std::size_t(src) % 32 == 0);
    return { _mm256_load_ps(src) };
}

static inline simd_float simd_loadu(const float* src) {
    return { _mm256_loadu_ps(src) };
}

static inline void simd_store(float* dst, simd_float& src) {
    // assert(std::size_t(dst) % 32 == 0);
    _mm256_store_ps(dst, src.internal_);
}

static inline void simd_storeu(float* dst, simd_float& src) {
    _mm256_storeu_ps(dst, src.internal_);
}

static inline simd_float operator+(simd_float a, simd_float b) {
    return { _mm256_add_ps(a.internal_, b.internal_) };
}

static inline simd_float operator-(simd_float a, simd_float b) {
    return { _mm256_sub_ps(a.internal_, b.internal_) };
}

static inline simd_float operator*(simd_float a, simd_float b) {
    return { _mm256_mul_ps(a.internal_, b.internal_) };
}

// a*b+c
static inline simd_float simd_fma(simd_float a, simd_float b, simd_float c) {
    return { _mm256_fmadd_ps(a.internal_, b.internal_, c.internal_) };
}

// a*b-c
static inline simd_float simd_fms(simd_float a, simd_float b, simd_float c) {
    return { _mm256_fmsub_ps(a.internal_, b.internal_, c.internal_) };
}

// native 1/sqrt(x)
static inline simd_float simd_rsqrt(simd_float x) {
    return { _mm256_rsqrt_ps(x.internal_) };
}

static inline simd_float simd_rsqrt_iter(simd_float lu, simd_float x) {
    simd_float tmp1 = x * lu;
    simd_float tmp2 = simd_float(-0.5f) * lu;
    tmp1 = simd_fma(tmp1, lu, simd_float(-3.0f));
    return tmp1 * tmp2;
}

// 1/sqrt(x) with higher precision
static inline simd_float simd_invsqrt(simd_float x) {
    simd_float lu = simd_rsqrt(x);
    lu = simd_rsqrt_iter(lu, x);
    return lu;
}

static inline float simd_reduce(simd_float a) {
    __m128 t0;
    t0 = _mm_add_ps(_mm256_castps256_ps128(a.internal_),
                    _mm256_extractf128_ps(a.internal_, 0x1));
    t0 = _mm_add_ps(t0, _mm_permute_ps(t0, _MM_SHUFFLE(1, 0, 3, 2)));
    t0 = _mm_add_ss(t0, _mm_permute_ps(t0, _MM_SHUFFLE(0, 3, 2, 1)));
    return *reinterpret_cast<float*>(&t0);
}

// xyz, xyz, ... => xx..., yy..., zz...
static inline void simd_loadu_xyz(const float* base, const std::int32_t* offset,
                                  simd_float& rx, simd_float& ry, simd_float& rz) {
    __m256 t1 = _mm256_insertf128_ps(
        _mm256_castps128_ps256(_mm_loadu_ps(base + 3*offset[0])),
        _mm_loadu_ps(base + 3*offset[4]),
        0x1
    );
    __m256 t2 = _mm256_insertf128_ps(
        _mm256_castps128_ps256(_mm_loadu_ps(base + 3*offset[1])),
        _mm_loadu_ps(base + 3*offset[5]),
        0x1
    );
    __m256 t3 = _mm256_insertf128_ps(
        _mm256_castps128_ps256(_mm_loadu_ps(base + 3*offset[2])),
        _mm_loadu_ps(base + 3*offset[6]),
        0x1
    );
    __m256 t4 = _mm256_insertf128_ps(
        _mm256_castps128_ps256(_mm_loadu_ps(base + 3*offset[3])),
        _mm_loadu_ps(base + 3*offset[7]),
        0x1
    );
    __m256 t5 = _mm256_unpacklo_ps(t1, t2);
    __m256 t6 = _mm256_unpacklo_ps(t3, t4);
    __m256 t7 = _mm256_unpackhi_ps(t1, t2);
    __m256 t8 = _mm256_unpackhi_ps(t3, t4);
    rx.internal_ = _mm256_shuffle_ps(t5, t6, _MM_SHUFFLE(1, 0, 1, 0));
    ry.internal_ = _mm256_shuffle_ps(t5, t6, _MM_SHUFFLE(3, 2, 3, 2));
    rz.internal_ = _mm256_shuffle_ps(t7, t8, _MM_SHUFFLE(1, 0, 1, 0));
}

// xx..., yy..., zz... => +xyz, +xyz, ...
static inline void simd_incru_xyz(float* base, const std::int32_t* offset,
                                  simd_float rx, simd_float ry, simd_float rz) {
    __m256 t1  = _mm256_unpacklo_ps(ry.internal_, rz.internal_);
    __m256 t2  = _mm256_unpackhi_ps(ry.internal_, rz.internal_);
    __m256 t3  = _mm256_shuffle_ps(rx.internal_, t1, _MM_SHUFFLE(1, 0, 0, 0));
    __m256 t4  = _mm256_shuffle_ps(rx.internal_, t1, _MM_SHUFFLE(3, 2, 0, 1));
    __m256 t5  = _mm256_shuffle_ps(rx.internal_, t2, _MM_SHUFFLE(1, 0, 0, 2));
    __m256 t6 = _mm256_shuffle_ps(rx.internal_, t2, _MM_SHUFFLE(3, 2, 0, 3));

    __m128 tA = _mm256_castps256_ps128(t3);
    __m128 tB = _mm256_castps256_ps128(t4);
    __m128 tC = _mm256_castps256_ps128(t5);
    __m128 tD = _mm256_castps256_ps128(t6);
    __m128 tE = _mm256_extractf128_ps(t3, 0x1);
    __m128 tF = _mm256_extractf128_ps(t4, 0x1);
    __m128 tG = _mm256_extractf128_ps(t5, 0x1);
    __m128 tH = _mm256_extractf128_ps(t6, 0x1);

    __m128 tX = _mm_load_ss(base + 3*offset[0]);
    tX = _mm_loadh_pi(tX, reinterpret_cast<__m64*>(base + 3*offset[0] + 1));
    tX = _mm_add_ps(tX, tA);
    _mm_store_ss(base + 3*offset[0], tX);
    _mm_storeh_pi(reinterpret_cast<__m64*>(base + 3*offset[0] + 1), tX);

    tX = _mm_load_ss(base + 3*offset[1]);
    tX = _mm_loadh_pi(tX, reinterpret_cast<__m64*>(base + 3*offset[1] + 1));
    tX = _mm_add_ps(tX, tB);
    _mm_store_ss(base + 3*offset[1], tX);
    _mm_storeh_pi(reinterpret_cast<__m64*>(base + 3*offset[1] + 1), tX);

    tX = _mm_load_ss(base + 3*offset[2]);
    tX = _mm_loadh_pi(tX, reinterpret_cast<__m64*>(base + 3*offset[2] + 1));
    tX = _mm_add_ps(tX, tC);
    _mm_store_ss(base + 3*offset[2], tX);
    _mm_storeh_pi(reinterpret_cast<__m64*>(base + 3*offset[2] + 1), tX);

    tX = _mm_load_ss(base + 3*offset[3]);
    tX = _mm_loadh_pi(tX, reinterpret_cast<__m64*>(base + 3*offset[3] + 1));
    tX = _mm_add_ps(tX, tD);
    _mm_store_ss(base + 3*offset[3], tX);
    _mm_storeh_pi(reinterpret_cast<__m64*>(base + 3*offset[3] + 1), tX);

    tX = _mm_load_ss(base + 3*offset[4]);
    tX = _mm_loadh_pi(tX, reinterpret_cast<__m64*>(base + 3*offset[4] + 1));
    tX = _mm_add_ps(tX, tE);
    _mm_store_ss(base + 3*offset[4], tX);
    _mm_storeh_pi(reinterpret_cast<__m64*>(base + 3*offset[4] + 1), tX);

    tX = _mm_load_ss(base + 3*offset[5]);
    tX = _mm_loadh_pi(tX, reinterpret_cast<__m64*>(base + 3*offset[5] + 1));
    tX = _mm_add_ps(tX, tF);
    _mm_store_ss(base + 3*offset[5], tX);
    _mm_storeh_pi(reinterpret_cast<__m64*>(base + 3*offset[5] + 1), tX);

    tX = _mm_load_ss(base + 3*offset[6]);
    tX = _mm_loadh_pi(tX, reinterpret_cast<__m64*>(base + 3*offset[6] + 1));
    tX = _mm_add_ps(tX, tG);
    _mm_store_ss(base + 3*offset[6], tX);
    _mm_storeh_pi(reinterpret_cast<__m64*>(base + 3*offset[6] + 1), tX);

    tX = _mm_load_ss(base + 3*offset[7]);
    tX = _mm_loadh_pi(tX, reinterpret_cast<__m64*>(base + 3*offset[7] + 1));
    tX = _mm_add_ps(tX, tH);
    _mm_store_ss(base + 3*offset[7], tX);
    _mm_storeh_pi(reinterpret_cast<__m64*>(base + 3*offset[7] + 1), tX);
}

// xx..., yy..., zz... => -xyz, -xyz, ...
static inline void simd_decru_xyz(float* base, const std::int32_t* offset,
                                  simd_float rx, simd_float ry, simd_float rz) {
    __m256 t1  = _mm256_unpacklo_ps(ry.internal_, rz.internal_);
    __m256 t2  = _mm256_unpackhi_ps(ry.internal_, rz.internal_);
    __m256 t3  = _mm256_shuffle_ps(rx.internal_, t1, _MM_SHUFFLE(1, 0, 0, 0));
    __m256 t4  = _mm256_shuffle_ps(rx.internal_, t1, _MM_SHUFFLE(3, 2, 0, 1));
    __m256 t5  = _mm256_shuffle_ps(rx.internal_, t2, _MM_SHUFFLE(1, 0, 0, 2));
    __m256 t6 = _mm256_shuffle_ps(rx.internal_, t2, _MM_SHUFFLE(3, 2, 0, 3));

    __m128 tA = _mm256_castps256_ps128(t3);
    __m128 tB = _mm256_castps256_ps128(t4);
    __m128 tC = _mm256_castps256_ps128(t5);
    __m128 tD = _mm256_castps256_ps128(t6);
    __m128 tE = _mm256_extractf128_ps(t3, 0x1);
    __m128 tF = _mm256_extractf128_ps(t4, 0x1);
    __m128 tG = _mm256_extractf128_ps(t5, 0x1);
    __m128 tH = _mm256_extractf128_ps(t6, 0x1);

    __m128 tX = _mm_load_ss(base + 3*offset[0]);
    tX = _mm_loadh_pi(tX, reinterpret_cast<__m64*>(base + 3*offset[0] + 1));
    tX = _mm_sub_ps(tX, tA);
    _mm_store_ss(base + 3*offset[0], tX);
    _mm_storeh_pi(reinterpret_cast<__m64*>(base + 3*offset[0] + 1), tX);

    tX = _mm_load_ss(base + 3*offset[1]);
    tX = _mm_loadh_pi(tX, reinterpret_cast<__m64*>(base + 3*offset[1] + 1));
    tX = _mm_sub_ps(tX, tB);
    _mm_store_ss(base + 3*offset[1], tX);
    _mm_storeh_pi(reinterpret_cast<__m64*>(base + 3*offset[1] + 1), tX);

    tX = _mm_load_ss(base + 3*offset[2]);
    tX = _mm_loadh_pi(tX, reinterpret_cast<__m64*>(base + 3*offset[2] + 1));
    tX = _mm_sub_ps(tX, tC);
    _mm_store_ss(base + 3*offset[2], tX);
    _mm_storeh_pi(reinterpret_cast<__m64*>(base + 3*offset[2] + 1), tX);

    tX = _mm_load_ss(base + 3*offset[3]);
    tX = _mm_loadh_pi(tX, reinterpret_cast<__m64*>(base + 3*offset[3] + 1));
    tX = _mm_sub_ps(tX, tD);
    _mm_store_ss(base + 3*offset[3], tX);
    _mm_storeh_pi(reinterpret_cast<__m64*>(base + 3*offset[3] + 1), tX);

    tX = _mm_load_ss(base + 3*offset[4]);
    tX = _mm_loadh_pi(tX, reinterpret_cast<__m64*>(base + 3*offset[4] + 1));
    tX = _mm_sub_ps(tX, tE);
    _mm_store_ss(base + 3*offset[4], tX);
    _mm_storeh_pi(reinterpret_cast<__m64*>(base + 3*offset[4] + 1), tX);

    tX = _mm_load_ss(base + 3*offset[5]);
    tX = _mm_loadh_pi(tX, reinterpret_cast<__m64*>(base + 3*offset[5] + 1));
    tX = _mm_sub_ps(tX, tF);
    _mm_store_ss(base + 3*offset[5], tX);
    _mm_storeh_pi(reinterpret_cast<__m64*>(base + 3*offset[5] + 1), tX);

    tX = _mm_load_ss(base + 3*offset[6]);
    tX = _mm_loadh_pi(tX, reinterpret_cast<__m64*>(base + 3*offset[6] + 1));
    tX = _mm_sub_ps(tX, tG);
    _mm_store_ss(base + 3*offset[6], tX);
    _mm_storeh_pi(reinterpret_cast<__m64*>(base + 3*offset[6] + 1), tX);

    tX = _mm_load_ss(base + 3*offset[7]);
    tX = _mm_loadh_pi(tX, reinterpret_cast<__m64*>(base + 3*offset[7] + 1));
    tX = _mm_sub_ps(tX, tH);
    _mm_store_ss(base + 3*offset[7], tX);
    _mm_storeh_pi(reinterpret_cast<__m64*>(base + 3*offset[7] + 1), tX);
}

}
