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

#include "test_lib.h"

namespace bytedock {

TEST(RandomTest, GenerateSeed) {
    int min_val = 1000000;
    int seed = auto_seed(min_val);
    EXPECT_GE(seed, min_val);
    min_val = -10;
    seed = auto_seed(min_val);
    EXPECT_GE(seed, min_val);
}

class RandomGeneratorTest : public testing::Test {
protected:
    void SetUp() override {
        rg = std::make_unique<random_generator>(59);
    }

    std::unique_ptr<random_generator> rg;
};

TEST_F(RandomGeneratorTest, GenerateInBox) {
    param_t init_xyz[3], oppo_xyz[3], rand_xyz[3];
    for (size_t i = 0; i < 5; ++i) {
        rg->random_in_sphere(init_xyz);
        scale_inplace_3d(init_xyz, i+1);
        for (size_t j = 0; j < 3; ++j) {
            oppo_xyz[j] = init_xyz[j] + rg->uniform_int(1, 10);
        }
        rg->random_in_box(init_xyz, oppo_xyz, rand_xyz);
        for (size_t j = 0; j < 3; ++j) {
            EXPECT_GT(rand_xyz[j], init_xyz[j]);
            EXPECT_LT(rand_xyz[j], oppo_xyz[j]);
        }
    }
}

TEST_F(RandomGeneratorTest, GenerateOrientation) {
    param_t qt[4], l2_norm;
    for (size_t i = 0; i < 5; ++i) {
        rg->random_orientation(qt);
        l2_norm = MATH_SQUARE_SCALAR(qt[0]) + MATH_SQUARE_SCALAR(qt[1])
                + MATH_SQUARE_SCALAR(qt[2]) + MATH_SQUARE_SCALAR(qt[3]);
        EXPECT_NEAR(1_r, l2_norm, 1e-7);
    }
}

TEST_F(RandomGeneratorTest, GenerateInSphere) {
    param_t xyz[3], l2_norm, prev = -1_r;
    for (size_t i = 0; i < 5; ++i) {
        rg->random_in_sphere(xyz);
        l2_norm = get_norm_3d(xyz);
        EXPECT_LT(l2_norm, 1_r);
        EXPECT_NE(l2_norm, prev);
        prev = l2_norm;
    }
}

}
