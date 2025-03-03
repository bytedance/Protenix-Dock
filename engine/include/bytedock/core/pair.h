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

#include "bytedock/core/system.h"
#include "bytedock/lib/math.h"
#include "bytedock/simd/avx2_def.h"

namespace bytedock {

inline void precalculate_coul_pair(const std::vector<param_t>& charges,
                                   nonbonded_atom_pair& item) {
    item.qq = charges[item.ids[0]] * charges[item.ids[1]];
}

inline void precalculate_vdw_pair(const lj_vdw* type_i, const lj_vdw* type_j,
                                  nonbonded_atom_pair& item) {
    param_t sigma6 = std::pow((type_i->sigma + type_j->sigma) * 0.5_r, 6_r);
    item.c6 = 4_r * std::sqrt(type_i->epsilon * type_j->epsilon) * sigma6;
    item.c12 = item.c6 * sigma6;
}

inline param_t calculate_coul_pair(const param_t qq,
                                   const param_t* r_ij,
                                   const param_t ir,  // 1/distance
                                   const param_t scale,
                                   atom_position& grad_i,
                                   atom_position& grad_j) {
    param_t qqk = qq * kCoulombFactor * scale;
    param_t de_dr_r = -qqk * MATH_CUBE_SCALAR(ir);
    param_t grad[3];
    scale_3d(r_ij, de_dr_r, grad);
    minus_inplace_3d(grad_i.xyz, grad);
    add_inplace_3d(grad_j.xyz, grad);
    return qqk * ir;
}

inline param_t calculate_coul_pair(const param_t qi,
                                   const param_t qj,
                                   const param_t* r_ij,
                                   const param_t distance,
                                   const param_t scale,
                                   atom_position& grad_i,
                                   atom_position& grad_j) {
    param_t qi_qj = qi * qj * scale;
    param_t de_dr_r = -safe_divide(qi_qj, MATH_CUBE_SCALAR(distance)) * kCoulombFactor;
    param_t grad[3];
    scale_3d(r_ij, de_dr_r, grad);
    minus_inplace_3d(grad_i.xyz, grad);
    add_inplace_3d(grad_j.xyz, grad);
    return safe_divide(qi_qj, distance) * kCoulombFactor;
}

inline param_t calculate_vdw_pair(const param_t c6,
                                  const param_t c12,
                                  const param_t* r_ij,
                                  const param_t ir,  // 1/distance
                                  const param_t scale,
                                  atom_position& grad_i,
                                  atom_position& grad_j) {
    param_t ir2 = MATH_SQUARE_SCALAR(ir);
    param_t ir6 = MATH_CUBE_SCALAR(ir2);
    param_t u6 = c6 * ir6 * scale;
    param_t u12 = c12 * MATH_SQUARE_SCALAR(ir6) * scale;
    param_t de_dr_r = (-12_r * u12 + 6_r * u6) * ir2;
    param_t grad[3];
    scale_3d(r_ij, de_dr_r, grad);
    minus_inplace_3d(grad_i.xyz, grad);
    add_inplace_3d(grad_j.xyz, grad);
    return u12 - u6;
}

inline param_t calculate_vdw_pair(const lj_vdw& type_i,
                                  const lj_vdw& type_j,
                                  const param_t* r_ij,
                                  const param_t distance,
                                  const param_t scale,
                                  atom_position& grad_i,
                                  atom_position& grad_j) {
    param_t fused = (type_i.sigma + type_j.sigma) * 0.5_r;
    param_t u6 = std::pow(safe_divide(fused, distance), 6_r);
    param_t u12 = MATH_SQUARE_SCALAR(u6);
    fused = std::sqrt(type_i.epsilon * type_j.epsilon) * scale * 4_r;
    param_t irsqr = safe_divide(1_r, MATH_SQUARE_SCALAR(distance));
    param_t de_dr_r = (-12_r * u12 + 6_r * u6) * irsqr * fused;
    param_t grad[3];
    scale_3d(r_ij, de_dr_r, grad);
    minus_inplace_3d(grad_i.xyz, grad);
    add_inplace_3d(grad_j.xyz, grad);
    return (u12 - u6) * fused;
}

#ifdef ENABLE_SIMD_AVX2
static void calculate_nb_pairs(const std::vector<nonbonded_atom_pair>& pairs,
                               const param_t coul_scale, param_t vdw_scale,
                               const molecule_pose& mol_xyz,
                               param_t& m_coul_e, param_t& m_vdw_e,
                               molecule_pose& xyz_gradient) {
    auto coords = reinterpret_cast<const param_t*>(mol_xyz.data());
    auto grads = reinterpret_cast<param_t*>(xyz_gradient.data());

    // Constant values
    simd_real cf(kCoulombFactor * coul_scale);
    simd_real vf(vdw_scale);
    simd_real six(6_r);
    simd_real twelve(12_r);

    // Variables for SIMD
    alignas(BDOCK_SIMD_ALIGNMENT) index_t m_ai[BDOCK_SIMD_REAL_WIDTH];
    alignas(BDOCK_SIMD_ALIGNMENT) index_t m_aj[BDOCK_SIMD_REAL_WIDTH];
    alignas(BDOCK_SIMD_ALIGNMENT) param_t m_qq[BDOCK_SIMD_REAL_WIDTH];
    alignas(BDOCK_SIMD_ALIGNMENT) param_t m_c6[BDOCK_SIMD_REAL_WIDTH];
    alignas(BDOCK_SIMD_ALIGNMENT) param_t m_c12[BDOCK_SIMD_REAL_WIDTH];

    // Calculate by pack
    const size_t npairs = pairs.size();
    simd_real coul_e(0_r), vdw_e(0_r);
    for (size_t i = 0; i < npairs; i += BDOCK_SIMD_REAL_WIDTH) {
        size_t s, iu;

        // Assemble inputs
        for (s = 0, iu = i; s < BDOCK_SIMD_REAL_WIDTH; ++s, ++iu) {
            if (iu == npairs) break;
            auto& item = pairs[iu];
            m_ai[s] = item.ids[0];
            m_aj[s] = item.ids[1];
            m_qq[s] = item.qq;
            m_c6[s] = item.c6;
            m_c12[s] = item.c12;
        }

        // Apply padding
        iu = s - 1;  // It is ensured that s>=1
        for (; s < BDOCK_SIMD_REAL_WIDTH; ++s) {
            m_ai[s] = m_ai[iu];
            m_aj[s] = m_aj[iu];
            m_qq[s] = m_c6[s] = m_c12[s] = 0_r;
        }

        // Load SIMD variables
        simd_real xi, yi, zi, xj, yj, zj;
        simd_loadu_xyz(coords, m_ai, xi, yi, zi);
        simd_loadu_xyz(coords, m_aj, xj, yj, zj);
        auto qqk = simd_load(m_qq) * cf;
        auto c6 = simd_load(m_c6);
        auto c12 = simd_load(m_c12);

        // Calculate distance & its inverse
        auto dx = xj - xi;
        auto dy = yj - yi;
        auto dz = zj - zi;
        auto rinv = simd_invsqrt(dx*dx + dy*dy + dz*dz);
        auto rinv2 = rinv * rinv;
        auto rinv6 = rinv2 * rinv2 * rinv2;

        // Calculate the Coulomb part
        coul_e = simd_fma(qqk, rinv, coul_e);
        auto de_dr_r = qqk * rinv2 * rinv;

        // Calculate the VDW part
        rinv = rinv6 * vf;  // Use rinv as tmp
        auto u6 = c6 * rinv;
        auto u12 = c12 * rinv6 * rinv;
        vdw_e = u12 - u6 + vdw_e;
        rinv6 = simd_fms(twelve, u12, six * u6);  // Use rinv6 as tmp
        de_dr_r = simd_fma(rinv6, rinv2, de_dr_r);
        auto gx = de_dr_r * dx;  // atom i
        auto gy = de_dr_r * dy;  // atom i
        auto gz = de_dr_r * dz;  // atom i

        // Update gradients
        simd_incru_xyz(grads, m_ai, gx, gy, gz);
        simd_decru_xyz(grads, m_aj, gx, gy, gz);
    }
    m_coul_e += simd_reduce(coul_e);
    m_vdw_e += simd_reduce(vdw_e);
}
#endif

}
