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

#include <cmath>

#include "bytedock/lib/dtype.h"

namespace bytedock {

struct vector_3d {
    param_t xyz[3];
};

struct vector_4d {
    param_t abcd[4];
};

struct matrix_3x3 {
    param_t data[9];
};

inline matrix_3x3 multiply_matrix_3x3(const matrix_3x3& r, const matrix_3x3& s) {
    matrix_3x3 out;
    out.data[0] = r.data[0]*s.data[0] + r.data[1]*s.data[3] + r.data[2]*s.data[6];
    out.data[1] = r.data[0]*s.data[1] + r.data[1]*s.data[4] + r.data[2]*s.data[7];
    out.data[2] = r.data[0]*s.data[2] + r.data[1]*s.data[5] + r.data[2]*s.data[8];
    out.data[3] = r.data[3]*s.data[0] + r.data[4]*s.data[3] + r.data[5]*s.data[6];
    out.data[4] = r.data[3]*s.data[1] + r.data[4]*s.data[4] + r.data[5]*s.data[7];
    out.data[5] = r.data[3]*s.data[2] + r.data[4]*s.data[5] + r.data[5]*s.data[8];
    out.data[6] = r.data[6]*s.data[0] + r.data[7]*s.data[3] + r.data[8]*s.data[6];
    out.data[7] = r.data[6]*s.data[1] + r.data[7]*s.data[4] + r.data[8]*s.data[7];
    out.data[8] = r.data[6]*s.data[2] + r.data[7]*s.data[5] + r.data[8]*s.data[8];
    return out;
}

// a-b => c
inline void minus_3d(const param_t* a, const param_t* b, param_t* c) {
    c[0] = a[0] - b[0];
    c[1] = a[1] - b[1];
    c[2] = a[2] - b[2];
}

// a+b => c
inline void add_3d(const param_t* a, const param_t* b, param_t* c) {
    c[0] = a[0] + b[0];
    c[1] = a[1] + b[1];
    c[2] = a[2] + b[2];
}

inline void add_inplace_3d(param_t* vector, const param_t* other) {
    vector[0] += other[0];
    vector[1] += other[1];
    vector[2] += other[2];
}

inline void minus_inplace_3d(param_t* vector, const param_t* other) {
    vector[0] -= other[0];
    vector[1] -= other[1];
    vector[2] -= other[2];
}

inline void multiply_3d(const param_t* vin, const param_t* other, param_t* vout) {
    vout[0] = vin[0] * other[0];
    vout[1] = vin[1] * other[1];
    vout[2] = vin[2] * other[2];
}

inline void scale_3d(const param_t* vin, const param_t factor, param_t* vout) {
    vout[0] = vin[0] * factor;
    vout[1] = vin[1] * factor;
    vout[2] = vin[2] * factor;
}

inline void scale_inplace_3d(param_t* vector, const param_t factor) {
    vector[0] *= factor;
    vector[1] *= factor;
    vector[2] *= factor;
}

// Dot product
inline param_t dot_product_3d(const param_t* a, const param_t* b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

inline void multiply_3x3_3d(const param_t* matrix, const param_t* vin, param_t* vout) {
    vout[0] = dot_product_3d(matrix, vin);
    vout[1] = dot_product_3d(matrix+3, vin);
    vout[2] = dot_product_3d(matrix+6, vin);
}

#define MATH_SQUARE_SCALAR(x) ((x)*(x))
#define MATH_CUBE_SCALAR(x) ((x)*(x)*(x))

inline param_t square(const param_t x) {
    return x * x;
}

inline param_t square_sum_3d(const param_t* vector) {
    param_t sum = MATH_SQUARE_SCALAR(vector[0]);
    sum += MATH_SQUARE_SCALAR(vector[1]);
    sum += MATH_SQUARE_SCALAR(vector[2]);
    return sum;
}

inline param_t get_norm_3d(const param_t* vector) {
    param_t sum = MATH_SQUARE_SCALAR(vector[0]);
    sum += MATH_SQUARE_SCALAR(vector[1]);
    sum += MATH_SQUARE_SCALAR(vector[2]);
    return std::sqrt(sum);
}

const param_t kParamEpsilon = 1e-20;

inline void normalize_3d(param_t* vector) {
    param_t sum = get_norm_3d(vector);
    vector[0] /= sum;
    vector[1] /= sum;
    vector[2] /= sum;
}

void get_rotation_matrix(const vector_4d& qt, matrix_3x3& mat);

void get_rotation_matrix(const vector_3d& axis, const param_t theta, matrix_3x3& out);

/**
 * Given one 3D vector to rotate and two 3x3 matrix:
 * - `out @ vec` => rotated vector
 * - `jbf @ vec` => jacobian over theta
 */
void get_rotation_matrix_and_gradient(const vector_3d& axis, const param_t theta,
                                      matrix_3x3& out, matrix_3x3& jbf);

#define MATH_PAIR_MAX(x, y) (x) > (y) ? (x) : (y)
#define MATH_PAIR_MIN(x, y) (x) < (y) ? (x) : (y)

inline param_t get_distance_square_3d(const param_t* xyz1, const param_t* xyz2) {
    param_t tmp = 0_r;
    tmp += MATH_SQUARE_SCALAR(xyz1[0] - xyz2[0]);
    tmp += MATH_SQUARE_SCALAR(xyz1[1] - xyz2[1]);
    tmp += MATH_SQUARE_SCALAR(xyz1[2] - xyz2[2]);
    return tmp;
}

inline param_t get_distance_3d(const param_t* xyz1, const param_t* xyz2) {
    param_t tmp = 0_r;
    tmp += MATH_SQUARE_SCALAR(xyz1[0] - xyz2[0]);
    tmp += MATH_SQUARE_SCALAR(xyz1[1] - xyz2[1]);
    tmp += MATH_SQUARE_SCALAR(xyz1[2] - xyz2[2]);
    return std::sqrt(tmp);
}

inline param_t safe_sqrt(const param_t x) {
    return std::sqrt(x + kParamEpsilon);
}

inline param_t safe_ln(const param_t x) {
    return std::log(x + kParamEpsilon);
}

inline param_t safe_divide(const param_t dividend, const param_t divisor) {
    return dividend / (divisor + kParamEpsilon);
}

inline param_t clamp(const param_t value, const param_t min, const param_t max) {
    return std::min(MATH_PAIR_MAX(value, min), max);
}

inline param_t get_linear_percent(
    const param_t value, const param_t minv, const param_t maxv
) {
    if (value > maxv) {
        return 1_r;
    } else if (value >= minv) {
        return (value - minv) / (maxv - minv);
    } else {
        return 0_r;
    }
}

inline param_t get_reversed_linear_percent(
    const param_t value, const param_t minv, const param_t maxv
) {
    if (value < minv) {
        return 1_r;
    } else if (value <= maxv) {
        return 1_r - (value - minv) / (maxv - minv);
    } else {
        return 0_r;
    }
}

inline param_t get_reversed_linear_percent(const param_t value, const param_t center,
                                           const param_t minv, const param_t maxv) {
    if (value < center) {
        return get_reversed_linear_percent(center - value, minv, maxv);
    } else {
        return get_reversed_linear_percent(value - center, minv, maxv);
    }
}

inline void cross_product_3d(const param_t* v1, const param_t* v2, param_t* out) {
    out[0] = v1[1] * v2[2] - v1[2] * v2[1];
    out[1] = v1[2] * v2[0] - v1[0] * v2[2];
    out[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

/**
 * Copy definition of macro `M_PI` from `math.h`.
 * Its value is the same as Python constant `math.pi`.
 */
const param_t kMathPi = 3.141592653589793;
const param_t kMathTwoPi = 2 * kMathPi;

#define MATH_TO_DEGREE(x) ((x) / kMathPi * 180_r)
#define MATH_TO_RADIAN(x) ((x) / 180_r * kMathPi)

inline param_t get_normal_angle(const vector_3d& nv1, const vector_3d& nv2) {
    param_t tmp = dot_product_3d(nv1.xyz, nv2.xyz);
    tmp = clamp(tmp, -1_r, 1_r);
    return std::acos(tmp);
}

// bytedock.utils.geometry.vecangle
inline param_t get_vector_angle(const vector_3d& v1, const vector_3d& v2) {
    param_t dp = dot_product_3d(v1.xyz, v2.xyz);
    param_t ll = get_norm_3d(v1.xyz) * get_norm_3d(v2.xyz);
    dp = clamp(safe_divide(dp, ll), -1_r, 1_r);
    return std::acos(dp);
}

// bytedock.utils.geometry.get_angle_vec(with_vec=False)
inline param_t get_angle_3d(const param_t* pos_i,
                            const param_t* pos_j,
                            const param_t* pos_k) {
    param_t v0[3], v1[3], cross[3];
    minus_3d(pos_j, pos_i, v0);
    minus_3d(pos_j, pos_k, v1);
    cross_product_3d(v0, v1, cross);
    return std::atan2(get_norm_3d(cross), dot_product_3d(v0, v1) + kParamEpsilon);
}

// bytedock.utils.geometry.get_angle_vec(with_vec=True)
param_t get_angle_and_gradient_3d(
    const param_t* pos_i, const param_t* pos_j, const param_t* pos_k,
    param_t* grad_i, param_t* grad_k
);

// bytedock.utils.geometry.get_dihedral_angle_vec(with_vec=False)
param_t get_dihedral_angle_3d(const param_t* pos_i, const param_t* pos_j,
                              const param_t* pos_k, const param_t* pos_l);

// bytedock.utils.geometry.get_dihedral_angle_vec(with_vec=True)
param_t get_dihedral_angle_and_gradient_3d(const param_t* pos_i, const param_t* pos_j,
                                           const param_t* pos_k, const param_t* pos_l,
                                           param_t* grad_i, param_t* grad_j,
                                           param_t* grad_k, param_t* grad_l);

}
