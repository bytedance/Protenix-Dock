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

#include "bytedock/core/pair.h"

#include "test_lib.h"

namespace bytedock {

TEST(PairPrecalculationTest, CoulumbEnergy) {
    const index_t ii = 1, jj = 3;
    std::vector<param_t> charges = {2_r, 3_r, 4_r, 5_r};
    nonbonded_atom_pair item = {ii, jj};
    precalculate_coul_pair(charges, item);
    param_t r_ij[] = {3_r, -4_r, 5_r};
    param_t distance = get_norm_3d(r_ij);
    param_t ir = safe_divide(1_r, distance);

    const param_t scale = 0.8_r;
    atom_position grad_i = {0_r}, grad_j = {0_r};
    param_t energy = calculate_coul_pair(charges[ii], charges[jj], r_ij, distance,
                                         scale, grad_i, grad_j);
    atom_position new_grad_i = {0_r}, new_grad_j = {0_r};
    param_t new_energy = calculate_coul_pair(item.qq, r_ij, ir, scale,
                                             new_grad_i, new_grad_j);
    EXPECT_NEAR(energy, new_energy, 1e-7);
    check_3d(grad_i.xyz, new_grad_i.xyz, 1e-7);
    check_3d(grad_j.xyz, new_grad_j.xyz, 1e-7);
}

TEST(PairPrecalculationTest, VdwEnergy) {
    lj_vdw type_i = {2_r, 0.4_r}; 
    lj_vdw type_j = {16_r, 0.8_r};
    nonbonded_atom_pair item;
    precalculate_vdw_pair(&type_i, &type_j, item);
    param_t r_ij[] = {-5_r, 7_r, -11_r};
    param_t distance = get_norm_3d(r_ij);
    param_t ir = safe_divide(1_r, distance);

    const param_t scale = 1.6_r;
    atom_position grad_i = {0_r}, grad_j = {0_r};
    param_t energy = calculate_vdw_pair(type_i, type_j, r_ij, distance, scale,
                                        grad_i, grad_j);
    atom_position new_grad_i = {0_r}, new_grad_j = {0_r};
    param_t new_energy = calculate_vdw_pair(item.c6, item.c12, r_ij, ir, scale,
                                            new_grad_i, new_grad_j);
    EXPECT_NEAR(energy, new_energy, 1e-7);
    check_3d(grad_i.xyz, new_grad_i.xyz, 1e-7);
    check_3d(grad_j.xyz, new_grad_j.xyz, 1e-7);
}

}
