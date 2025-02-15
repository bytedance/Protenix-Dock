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

#include "bytedock/lib/math.h"

namespace bytedock {

void get_rotation_matrix(const vector_4d& qt, matrix_3x3& mat) {
    const param_t aa = MATH_SQUARE_SCALAR(qt.abcd[0]);
    const param_t ab = qt.abcd[0] * qt.abcd[1];
    const param_t ac = qt.abcd[0] * qt.abcd[2];
    const param_t ad = qt.abcd[0] * qt.abcd[3];
    const param_t bb = MATH_SQUARE_SCALAR(qt.abcd[1]);
    const param_t bc = qt.abcd[1] * qt.abcd[2];
    const param_t bd = qt.abcd[1] * qt.abcd[3];
    const param_t cc = MATH_SQUARE_SCALAR(qt.abcd[2]);
    const param_t cd = qt.abcd[2] * qt.abcd[3];
    const param_t dd = MATH_SQUARE_SCALAR(qt.abcd[3]);
    mat.data[0] = aa + bb - cc - dd;  // (0, 0)
    mat.data[1] = 2 * (bc - ad);      // (0, 1)
    mat.data[2] = 2 * (ac + bd);      // (0, 2)
    mat.data[3] = 2 * (ad + bc);      // (1, 0)
    mat.data[4] = aa - bb + cc - dd;  // (1, 1)
    mat.data[5] = 2 * (cd - ab);      // (1, 2)
    mat.data[6] = 2 * (bd - ac);      // (2, 0)
    mat.data[7] = 2 * (ab + cd);      // (2, 1)
    mat.data[8] = aa - bb - cc + dd;  // (2, 2)
}

void get_rotation_matrix(const vector_3d& axis, const param_t theta, matrix_3x3& out) {
    const param_t& kx = axis.xyz[0];
    const param_t& ky = axis.xyz[1];
    const param_t& kz = axis.xyz[2];
    param_t st = std::sin(theta);
    param_t ct = std::cos(theta);
    param_t one_minus_ct = 1 - ct;
    out.data[0] = ct + MATH_SQUARE_SCALAR(kx) * one_minus_ct;
    out.data[1] = kx * ky * one_minus_ct - kz * st;
    out.data[2] = ky * st + kx * kz * one_minus_ct;
    out.data[3] = kz * st + kx * ky * one_minus_ct;
    out.data[4] = ct + MATH_SQUARE_SCALAR(ky) * one_minus_ct;
    out.data[5] = ky * kz * one_minus_ct - kx * st;
    out.data[6] = kx * kz * one_minus_ct - ky * st;
    out.data[7] = kx * st + ky * kz * one_minus_ct;
    out.data[8] = ct + MATH_SQUARE_SCALAR(kz) * one_minus_ct;
}

void get_rotation_matrix_and_gradient(const vector_3d& axis, const param_t theta,
                                      matrix_3x3& out, matrix_3x3& jbf) {
    const param_t& kx = axis.xyz[0];
    const param_t& ky = axis.xyz[1];
    const param_t& kz = axis.xyz[2];
    param_t st = std::sin(theta);
    param_t ct = std::cos(theta);
    param_t one_minus_ct = 1 - ct;
    param_t kx_omc = kx * one_minus_ct;
    param_t ky_omc = ky * one_minus_ct;
    param_t kz_omc = kz * one_minus_ct;
    param_t kx_st = kx * st;
    param_t ky_st = ky * st;
    param_t kz_st = kz * st;
    out.data[0] = ct + kx * kx_omc;
    out.data[1] = ky * kx_omc - kz_st;
    out.data[2] = ky_st + kz * kx_omc;
    out.data[3] = kz_st + kx * ky_omc;
    out.data[4] = ct + ky * ky_omc;
    out.data[5] = kz * ky_omc - kx_st;
    out.data[6] = kx * kz_omc - ky_st;
    out.data[7] = kx_st + ky * kz_omc;
    out.data[8] = ct + kz * kz_omc;
    jbf.data[0] = -st + kx * kx_st;
    jbf.data[1] = -kz * ct + kx * ky_st;
    jbf.data[2] = ky * ct + kx * kz_st;
    jbf.data[3] = kz * ct + ky * kx_st;
    jbf.data[4] = -st + ky * ky_st;
    jbf.data[5] = -kx * ct + ky * kz_st;
    jbf.data[6] = -ky * ct + kz * kx_st;
    jbf.data[7] = kx * ct + kz * ky_st;
    jbf.data[8] = -st + kz * kz_st;
}

param_t get_angle_and_gradient_3d(
    const param_t* pos_i, const param_t* pos_j, const param_t* pos_k,
    param_t* grad_i, param_t* grad_k
) {
    param_t v0[3], v1[3], cross[3];
    minus_3d(pos_j, pos_i, v0);
    minus_3d(pos_j, pos_k, v1);
    cross_product_3d(v0, v1, cross);
    param_t rp = get_norm_3d(cross);
    cross_product_3d(v0, cross, grad_i);
    cross_product_3d(cross, v1, grad_k);
    scale_inplace_3d(grad_i, -safe_divide(1_r, square_sum_3d(v0)*rp));
    scale_inplace_3d(grad_k, -safe_divide(1_r, square_sum_3d(v1)*rp));
    return std::atan2(rp, dot_product_3d(v0, v1) + kParamEpsilon);
}

param_t get_dihedral_angle_3d(const param_t* pos_i, const param_t* pos_j,
                             const param_t* pos_k, const param_t* pos_l) {
    param_t r_ij[3], r_kj[3], w[3];  // w => r_kl
    minus_3d(pos_j, pos_i, r_ij);
    minus_3d(pos_j, pos_k, r_kj);
    minus_3d(pos_l, pos_k, w);
    param_t m[3], n[3];
    cross_product_3d(r_ij, r_kj, m);
    cross_product_3d(r_kj, w, n);
    cross_product_3d(m, n, w);
    param_t wlen = get_norm_3d(w);
    // w[3] => {s, phi, ipr}
    w[0] = dot_product_3d(m, n);
    w[1] = std::atan2(wlen, w[0] + kParamEpsilon);
    w[2] = dot_product_3d(r_ij, n);
    return w[2] < 0_r ? w[1] : (-w[1]);
}

param_t get_dihedral_angle_and_gradient_3d(
    const param_t* pos_i, const param_t* pos_j,
    const param_t* pos_k, const param_t* pos_l,
    param_t* grad_i, param_t* grad_j, param_t* grad_k, param_t* grad_l
) {
    param_t r_ij[3], r_kj[3], r_kl[3];
    minus_3d(pos_j, pos_i, r_ij);
    minus_3d(pos_j, pos_k, r_kj);
    minus_3d(pos_l, pos_k, r_kl);
    param_t m[3], n[3], w[3];
    cross_product_3d(r_ij, r_kj, m);
    cross_product_3d(r_kj, r_kl, n);
    cross_product_3d(m, n, w);
    param_t wlen = get_norm_3d(w);
    // w[3] => {s, phi, ipr}
    w[0] = dot_product_3d(m, n);
    w[1] = std::atan2(wlen, w[0] + kParamEpsilon);
    w[2] = dot_product_3d(r_ij, n);
    wlen = w[2] < 0_r ? w[1] : (-w[1]);  // wlen => theta
    param_t a = square_sum_3d(r_kj);  // a => nrkj2
    param_t b = safe_divide(1_r, std::sqrt(a));  // b => nrkj_1
    param_t nrkj_2 = MATH_SQUARE_SCALAR(b);
    param_t nrkj = a * b;
    a = -nrkj / square_sum_3d(m);
    b = nrkj / square_sum_3d(n);
    scale_3d(m, -a, grad_i);
    scale_3d(n, -b, grad_l);
    a = dot_product_3d(r_ij, r_kj) * nrkj_2;  // a => p
    b = dot_product_3d(r_kj, r_kl) * nrkj_2;  // b => q
    scale_3d(grad_i, a, m);  // m => uvec
    scale_3d(grad_l, b, n);  // n => vvec
    minus_3d(m, n, w);  // w => svec
    minus_3d(grad_i, w, grad_j);
    add_3d(grad_l, w, grad_k);
    scale_inplace_3d(grad_j, -1_r);
    scale_inplace_3d(grad_k, -1_r);
    return wlen;
}

}
