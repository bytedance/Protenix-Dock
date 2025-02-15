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

#include "test_lib.h"

namespace bytedock {

TEST(AlgebraTest, MultiplyBetweenMatrix3D) {
    matrix_3x3 lhs, rhs;
    for (size_t i = 0; i < 9; ++i) {
        lhs.data[i] = i + 1;
        rhs.data[i] = 20_r - i;
    }
    auto out = multiply_matrix_3x3(lhs, rhs);
    EXPECT_NEAR(96_r, out.data[0], 1e-7);
    EXPECT_NEAR(90_r, out.data[1], 1e-7);
    EXPECT_NEAR(84_r, out.data[2], 1e-7);
    EXPECT_NEAR(249_r, out.data[3], 1e-7);
    EXPECT_NEAR(234_r, out.data[4], 1e-7);
    EXPECT_NEAR(219_r, out.data[5], 1e-7);
    EXPECT_NEAR(402_r, out.data[6], 1e-7);
    EXPECT_NEAR(378_r, out.data[7], 1e-7);
    EXPECT_NEAR(354_r, out.data[8], 1e-7);
}

TEST(GeometryTest, GetDihedralAngle) {
    param_t pos_i[3] = {7.908, 27.763, 9.369};
    param_t pos_j[3] = {8.890, 26.662, 8.959};
    param_t pos_k[3] = {9.472, 26.977, 7.565};
    param_t pos_l[3] = {8.806, 26.234, 6.393};
    param_t theta = get_dihedral_angle_3d(pos_i, pos_j, pos_k, pos_l);
    EXPECT_NEAR(-1.7239, theta, 3e-5);
}

TEST(GeometryTest, GetDihedralAngleAndGradient) {
    param_t pos_i[3] = {7.908, 27.763, 9.369};
    param_t pos_j[3] = {8.890, 26.662, 8.959};
    param_t pos_k[3] = {9.472, 26.977, 7.565};
    param_t pos_l[3] = {8.806, 26.234, 6.393};
    param_t grad_i[3], grad_j[3], grad_k[3], grad_l[3];
    param_t theta = get_dihedral_angle_and_gradient_3d(pos_i, pos_j, pos_k, pos_l,
                                                       grad_i, grad_j, grad_k, grad_l);
    EXPECT_NEAR(-1.7239, theta, 3e-5);
    param_t ref_i[3] = {-0.5188, -0.3524, -0.2962};
    param_t ref_j[3] = { 0.4927,  0.6991,  0.3637};
    param_t ref_k[3] = { 0.4956, -0.8849,  0.0070};
    param_t ref_l[3] = {-0.4695,  0.5383, -0.0744};
    check_3d(ref_i, grad_i, 1e-4);
    check_3d(ref_j, grad_j, 1e-4);
    check_3d(ref_k, grad_k, 1e-4);
    check_3d(ref_l, grad_l, 1e-4);
}

TEST(GeometryTest, GetRotationMatrix) {
    vector_3d axis = {1., 2., 3.};
    normalize_3d(axis.xyz);
    param_t vin[3] = {1., 1., 1.};
    param_t theta = MATH_TO_RADIAN(45.);
    matrix_3x3 rot;
    param_t tmp[3];
    get_rotation_matrix(axis, theta, rot);
    multiply_3x3_3d(rot.data, vin, tmp);
    param_t ref_out[3] = {0.6437, 1.3361, 0.8947};
    check_3d(ref_out, tmp, 5e-5);
}

TEST(GeometryTest, GetRotationMatrixAndGradient) {
    vector_3d axis = {0.26726124, 0.5345225, 0.8017837};
    param_t vin[3] = {1., 1., 1.};
    param_t theta = M_PI / 4.;
    matrix_3x3 rot, jbf;
    param_t tmp[3];
    get_rotation_matrix_and_gradient(axis, theta, rot, jbf);

    param_t ref_out[3] = {0.6437, 1.3361, 0.8947};
    multiply_3x3_3d(rot.data, vin, tmp);  // Rotated vector
    check_3d(ref_out, tmp, 5e-5);
    param_t ref_jcb[3] = {-0.5930, 0.2769, 0.0130};
    multiply_3x3_3d(jbf.data, vin, tmp);  // Jacobian over theta
    check_3d(ref_jcb, tmp, 5e-5);
}

}
