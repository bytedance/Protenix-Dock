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

#include <random>

#include "bytedock/lib/math.h"

namespace bytedock {

inline int auto_seed(int min_val) {
    std::random_device rd;
    std::uniform_int_distribution<> distr(min_val);
    return distr(rd);
}

// Generate a sequence of random numbers with the same engine
class random_generator {
public:
    random_generator(unsigned int seed) : engine_(seed), urd_pi_(-kMathPi, kMathPi),
                                          urd_pos_one_(0_r, 1_r),
                                          urd_one_(-1_r, 1_r) {}

    int uniform_int(int min_val, int max_val) {
        std::uniform_int_distribution<> uid(min_val, max_val);
        return uid(engine_);
    }

    param_t uniform_real(const param_t& min_val, const param_t& max_val) {
        std::uniform_real_distribution<> urd(min_val, max_val);
        return urd(engine_);
    }

    param_t uniform_torsion() {
        return urd_pi_(engine_);
    }

    void random_in_box(const param_t* init_xyz, const param_t* oppo_xyz,
                       param_t* rand_xyz);
    void random_orientation(param_t* quaternion);
    void random_in_sphere(param_t* xyz);

private:
    std::mt19937 engine_;
    std::uniform_real_distribution<param_t> urd_pi_;  // [-pi, pi)
    std::uniform_real_distribution<param_t> urd_pos_one_;  // [0, 1)
    std::uniform_real_distribution<param_t> urd_one_;  // [-1, 1)
};

}
