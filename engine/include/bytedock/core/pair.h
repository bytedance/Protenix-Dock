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

namespace bytedock {

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
    param_t energy = safe_divide(qi_qj, distance) * kCoulombFactor;
    return energy;
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

}
