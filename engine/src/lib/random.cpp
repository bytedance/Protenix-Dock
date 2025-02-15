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

#include "bytedock/lib/random.h"

namespace bytedock {

void random_generator::random_in_box(const param_t* init_xyz, const param_t* oppo_xyz,
                                     param_t* rand_xyz) {
    rand_xyz[0] = uniform_real(init_xyz[0], oppo_xyz[0]);
    rand_xyz[1] = uniform_real(init_xyz[1], oppo_xyz[1]);
    rand_xyz[2] = uniform_real(init_xyz[2], oppo_xyz[2]);
}

void random_generator::random_orientation(param_t* quaternion) {
    // Reference: https://stackoverflow.com/a/44031492/28276974
    param_t u1 = urd_pos_one_(engine_);
    param_t u2 = urd_pos_one_(engine_) * kMathTwoPi;
    param_t u3 = urd_pos_one_(engine_) * kMathTwoPi;
    param_t sqrt1u1 = std::sqrt(1_r - u1);
    param_t sqrtu1 = std::sqrt(u1);
    quaternion[0] = sqrt1u1 * std::sin(u2);
    quaternion[1] = sqrt1u1 * std::cos(u2);
    quaternion[2] = sqrtu1 * std::sin(u3);
    quaternion[3] = sqrtu1 * std::cos(u3);
}

void random_generator::random_in_sphere(param_t* xyz) {
    while (true) {
        xyz[0] = urd_one_(engine_);
        xyz[1] = urd_one_(engine_);
        xyz[2] = urd_one_(engine_);
        if (square_sum_3d(xyz) < 1_r) break;
    }
}

}
