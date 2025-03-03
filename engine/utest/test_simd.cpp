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

#ifdef ENABLE_SIMD_AVX2
#include "bytedock/simd/avx2_def.h"

#include "test_lib.h"

namespace bytedock {

TEST(SimdFloatTest, Intialization) {
    simd_float a(1.5f);
    alignas(BDOCK_SIMD_ALIGNMENT) float buffer[BDOCK_SIMD_FLOAT_WIDTH];
    simd_store(buffer, a);
    for (int i = 0; i < BDOCK_SIMD_FLOAT_WIDTH; ++i) {
        EXPECT_NEAR(buffer[i], 1.5, 1e-6);
    }
}

TEST(SimdDoubleTest, Intialization) {
    simd_double a(-1.5);
    alignas(BDOCK_SIMD_ALIGNMENT) double buffer[BDOCK_SIMD_DOUBLE_WIDTH];
    simd_store(buffer, a);
    for (int i = 0; i < BDOCK_SIMD_DOUBLE_WIDTH; ++i) {
        EXPECT_NEAR(buffer[i], -1.5, 1e-15);
    }
}

TEST(SimdFloatTest, Addition) {
    auto a = simd_float(-3.6f) + simd_float(2.9f);
    float buffer[BDOCK_SIMD_FLOAT_WIDTH];
    simd_storeu(buffer, a);
    for (int i = 0; i < BDOCK_SIMD_FLOAT_WIDTH; ++i) {
        EXPECT_NEAR(buffer[i], -0.7, 1e-6);
    }
}

TEST(SimdDoubleTest, Addition) {
    auto a = simd_double(-4.2) + simd_double(11.5);
    double buffer[BDOCK_SIMD_DOUBLE_WIDTH];
    simd_storeu(buffer, a);
    for (int i = 0; i < BDOCK_SIMD_DOUBLE_WIDTH; ++i) {
        EXPECT_NEAR(buffer[i], 7.3, 1e-15);
    }
}

TEST(SimdFloatTest, Subtraction) {
    auto a = simd_float(1.5f) - simd_float(3.2f);
    alignas(BDOCK_SIMD_ALIGNMENT) float buffer[BDOCK_SIMD_FLOAT_WIDTH];
    simd_store(buffer, a);
    for (int i = 0; i < BDOCK_SIMD_FLOAT_WIDTH; ++i) {
        EXPECT_NEAR(buffer[i], -1.7, 1e-6);
    }
}

TEST(SimdDoubleTest, Subtraction) {
    auto a = simd_double(9.5) - simd_double(-2.4);
    alignas(BDOCK_SIMD_ALIGNMENT) double buffer[BDOCK_SIMD_DOUBLE_WIDTH];
    simd_store(buffer, a);
    for (int i = 0; i < BDOCK_SIMD_DOUBLE_WIDTH; ++i) {
        EXPECT_NEAR(buffer[i], 11.9, 1e-15);
    }
}

TEST(SimdFloatTest, Multiplication) {
    auto a = simd_float(1.5f) * simd_float(3.2f);
    float buffer[BDOCK_SIMD_FLOAT_WIDTH];
    simd_storeu(buffer, a);
    for (int i = 0; i < BDOCK_SIMD_FLOAT_WIDTH; ++i) {
        EXPECT_NEAR(buffer[i], 4.8, 1e-6);
    }
}

TEST(SimdDoubleTest, Multiplication) {
    auto a = simd_double(9.5) * simd_double(-2.4);
    double buffer[BDOCK_SIMD_DOUBLE_WIDTH];
    simd_storeu(buffer, a);
    for (int i = 0; i < BDOCK_SIMD_DOUBLE_WIDTH; ++i) {
        EXPECT_NEAR(buffer[i], -22.8, 1e-15);
    }
}

TEST(SimdFloatTest, FusedMultiplyAdd) {
    auto a = simd_fma(simd_float(5.3f), simd_float(0.7f), simd_float(-9.4f));
    alignas(BDOCK_SIMD_ALIGNMENT) float buffer[BDOCK_SIMD_FLOAT_WIDTH];
    simd_store(buffer, a);
    for (int i = 0; i < BDOCK_SIMD_FLOAT_WIDTH; ++i) {
        EXPECT_NEAR(buffer[i], -5.69, 1e-6);
    }
}

TEST(SimdDoubleTest, FusedMultiplyAdd) {
    auto a = simd_fma(simd_double(7.2), simd_double(-1.4), simd_double(3.5));
    alignas(BDOCK_SIMD_ALIGNMENT) double buffer[BDOCK_SIMD_DOUBLE_WIDTH];
    simd_store(buffer, a);
    for (int i = 0; i < BDOCK_SIMD_DOUBLE_WIDTH; ++i) {
        EXPECT_NEAR(buffer[i], -6.58, 1e-15);
    }
}

TEST(SimdFloatTest, FusedMultiplySubtract) {
    auto a = simd_fms(simd_float(0.7f), simd_float(5.3f), simd_float(9.4f));
    float buffer[BDOCK_SIMD_FLOAT_WIDTH];
    simd_storeu(buffer, a);
    for (int i = 0; i < BDOCK_SIMD_FLOAT_WIDTH; ++i) {
        EXPECT_NEAR(buffer[i], -5.69, 1e-6);
    }
}

TEST(SimdDoubleTest, FusedMultiplySubtract) {
    auto a = simd_fms(simd_double(-1.4), simd_double(7.2), simd_double(-3.5));
    double buffer[BDOCK_SIMD_DOUBLE_WIDTH];
    simd_storeu(buffer, a);
    for (int i = 0; i < BDOCK_SIMD_DOUBLE_WIDTH; ++i) {
        EXPECT_NEAR(buffer[i], -6.58, 1e-15);
    }
}

TEST(SimdFloatTest, InvSqrt) {
    auto a = simd_invsqrt(simd_float(6.54321f));
    float buffer[BDOCK_SIMD_FLOAT_WIDTH];
    simd_storeu(buffer, a);
    for (int i = 0; i < BDOCK_SIMD_FLOAT_WIDTH; ++i) {
        EXPECT_NEAR(buffer[i], 0.390935, 1e-6);
    }
}

TEST(SimdDoubleTest, InvSqrt) {
    auto a = simd_invsqrt(simd_double(9.7654321));
    double buffer[BDOCK_SIMD_DOUBLE_WIDTH];
    simd_storeu(buffer, a);
    for (int i = 0; i < BDOCK_SIMD_DOUBLE_WIDTH; ++i) {
        EXPECT_NEAR(buffer[i], 0.320003160520422, 2e-15);
    }
}

TEST(SimdFloatTest, SumElements) {
    alignas(BDOCK_SIMD_ALIGNMENT) float buffer[] = {
        1.f, -2.f, 3.f, 4.f, -5.f, 6.f, 7.f, 8.f
    };
    auto a = simd_reduce(simd_load(buffer));
    EXPECT_NEAR(a, 22, 1e-6);
}

TEST(SimdDoubleTest, SumElements) {
    alignas(BDOCK_SIMD_ALIGNMENT) double buffer[] = {1., -12., 23., 74.};
    auto a = simd_reduce(simd_load(buffer));
    EXPECT_NEAR(a, 86, 1e-15);
}

TEST(SimdFloatTest, LoadXYZ) {
    float coords[BDOCK_SIMD_ALIGNMENT];  // Larger than 3*BDOCK_SIMD_FLOAT_WIDTH
    for (int i = 0; i < BDOCK_SIMD_ALIGNMENT; ++i) coords[i] = i;
    std::int32_t offsets[] = {9, 7, 6, 5, 3, 2, 1, 0};
    simd_float rx, ry, rz;
    simd_loadu_xyz(coords, offsets, rx, ry, rz);
    float m_rx[BDOCK_SIMD_FLOAT_WIDTH], m_ry[BDOCK_SIMD_FLOAT_WIDTH],
          m_rz[BDOCK_SIMD_FLOAT_WIDTH];
    simd_storeu(m_rx, rx);
    simd_storeu(m_ry, ry);
    simd_storeu(m_rz, rz);
    for (int i = 0; i < BDOCK_SIMD_FLOAT_WIDTH; ++i) {
        auto ref = coords + 3*offsets[i];
        EXPECT_NEAR(m_rx[i], ref[0], 1e-6);
        EXPECT_NEAR(m_ry[i], ref[1], 1e-6);
        EXPECT_NEAR(m_rz[i], ref[2], 1e-6);
    }
}

TEST(SimdDoubleTest, LoadXYZ) {
    double coords[BDOCK_SIMD_ALIGNMENT];  // Larger than 3*BDOCK_SIMD_DOUBLE_WIDTH
    for (int i = 0; i < BDOCK_SIMD_ALIGNMENT; ++i) coords[i] = i;
    std::int32_t offsets[] = {8, 5, 4, 1};
    simd_double rx, ry, rz;
    simd_loadu_xyz(coords, offsets, rx, ry, rz);
    double m_rx[BDOCK_SIMD_DOUBLE_WIDTH], m_ry[BDOCK_SIMD_DOUBLE_WIDTH],
           m_rz[BDOCK_SIMD_DOUBLE_WIDTH];
    simd_storeu(m_rx, rx);
    simd_storeu(m_ry, ry);
    simd_storeu(m_rz, rz);
    for (int i = 0; i < BDOCK_SIMD_DOUBLE_WIDTH; ++i) {
        auto ref = coords + 3*offsets[i];
        EXPECT_NEAR(m_rx[i], ref[0], 1e-15);
        EXPECT_NEAR(m_ry[i], ref[1], 1e-15);
        EXPECT_NEAR(m_rz[i], ref[2], 1e-15);
    }
}

TEST(SimdFloatTest, IncrXYZ) {
    float grads[BDOCK_SIMD_ALIGNMENT];  // Larger than 3*BDOCK_SIMD_FLOAT_WIDTH
    for (int i = 0; i < BDOCK_SIMD_ALIGNMENT; ++i) grads[i] = -12.f;
    std::int32_t offsets[] = {9, 7, 6, 5, 3, 2, 1, 0};
    alignas(BDOCK_SIMD_ALIGNMENT) float buffer[] = {
        1.f, 4.f, 7.f, 10.f, 13.f, 16.f, 19.f, 22.f,
        2.f, 5.f, 8.f, 11.f, 14.f, 17.f, 20.f, 23.f,
        3.f, 6.f, 9.f, 12.f, 15.f, 18.f, 21.f, 24.f
    };
    simd_float rx = simd_load(buffer);
    simd_float ry = simd_load(buffer + BDOCK_SIMD_FLOAT_WIDTH);
    simd_float rz = simd_load(buffer + 2*BDOCK_SIMD_FLOAT_WIDTH);
    simd_incru_xyz(grads, offsets, rx, ry, rz);
    float ref[] = {
         10.f,  11.f,  12.f,
          7.f,   8.f,   9.f,
          4.f,   5.f,   6.f,
          1.f,   2.f,   3.f,
        -12.f, -12.f, -12.f,
         -2.f,  -1.f,   0.f,
         -5.f,  -4.f,  -3.f,
         -8.f,  -7.f,  -6.f,
        -12.f, -12.f, -12.f,
        -11.f, -10.f,  -9.f,
        -12.f, -12.f
    };
    for (int i = 0; i < BDOCK_SIMD_ALIGNMENT; ++i) {
        EXPECT_NEAR(ref[i], grads[i], 1e-6);
    }
}

TEST(SimdDoubleTest, IncrXYZ) {
    double grads[BDOCK_SIMD_ALIGNMENT];  // Larger than 3*BDOCK_SIMD_DOUBLE_WIDTH
    for (int i = 0; i < BDOCK_SIMD_ALIGNMENT; ++i) grads[i] = 5.;
    std::int32_t offsets[] = {8, 5, 4, 1};
    alignas(BDOCK_SIMD_ALIGNMENT) double buffer[] = {
        -1.,  4.,  7., -10.,
         2., -5.,  8., -11.,
         3.,  6., -9.,  12.,
    };
    simd_double rx = simd_load(buffer);
    simd_double ry = simd_load(buffer + BDOCK_SIMD_DOUBLE_WIDTH);
    simd_double rz = simd_load(buffer + 2*BDOCK_SIMD_DOUBLE_WIDTH);
    simd_incru_xyz(grads, offsets, rx, ry, rz);
    double ref[] = {
         5.,  5.,  5.,
        -5., -6., 17.,
         5.,  5.,  5.,
         5.,  5.,  5.,
        12., 13., -4.,
         9.,  0., 11.,
         5.,  5.,  5.,
         5.,  5.,  5.,
         4.,  7.,  8.,
         5.,  5.,  5.,
         5.,  5.
    };
    for (int i = 0; i < BDOCK_SIMD_ALIGNMENT; ++i) {
        EXPECT_NEAR(ref[i], grads[i], 1e-15);
    }
}

TEST(SimdFloatTest, DecrXYZ) {
    float grads[BDOCK_SIMD_ALIGNMENT];  // Larger than 3*BDOCK_SIMD_FLOAT_WIDTH
    for (int i = 0; i < BDOCK_SIMD_ALIGNMENT; ++i) grads[i] = -9.f;
    std::int32_t offsets[] = {9, 8, 7, 5, 4, 3, 2, 1};
    float buffer[] = {
         1.f,  4.f, -7.f,  10.f, -13.f, -16.f, -19.f, 22.f,
        -2.f,  5.f,  8.f, -11.f, -14.f,  17.f, -20.f, 23.f,
         3.f, -6.f,  9.f, -12.f,  15.f, -18.f, -21.f, 24.f
    };
    simd_float rx = simd_loadu(buffer);
    simd_float ry = simd_loadu(buffer + BDOCK_SIMD_FLOAT_WIDTH);
    simd_float rz = simd_loadu(buffer + 2*BDOCK_SIMD_FLOAT_WIDTH);
    simd_decru_xyz(grads, offsets, rx, ry, rz);
    float ref[] = {
         -9.f,  -9.f,  -9.f,
        -31.f, -32.f, -33.f,
         10.f,  11.f,  12.f,
          7.f, -26.f,   9.f,
          4.f,   5.f, -24.f,
        -19.f,   2.f,   3.f,
         -9.f,  -9.f,  -9.f,
         -2.f, -17.f, -18.f,
        -13.f, -14.f,  -3.f,
        -10.f,  -7.f, -12.f,
         -9.f,  -9.f
    };
    for (int i = 0; i < BDOCK_SIMD_ALIGNMENT; ++i) {
        EXPECT_NEAR(ref[i], grads[i], 1e-6);
    }
}

TEST(SimdDoubleTest, DecrXYZ) {
    double grads[BDOCK_SIMD_ALIGNMENT];  // Larger than 3*BDOCK_SIMD_DOUBLE_WIDTH
    for (int i = 0; i < BDOCK_SIMD_ALIGNMENT; ++i) grads[i] = 17.;
    std::int32_t offsets[] = {7, 6, 3, 0};
    alignas(BDOCK_SIMD_ALIGNMENT) double buffer[] = {
         1.,  4., -7.,  10.,
        -2.,  5., -8.,  11.,
         3., -6.,  9., -12.,
    };
    simd_double rx = simd_loadu(buffer);
    simd_double ry = simd_loadu(buffer + BDOCK_SIMD_DOUBLE_WIDTH);
    simd_double rz = simd_loadu(buffer + 2*BDOCK_SIMD_DOUBLE_WIDTH);
    simd_decru_xyz(grads, offsets, rx, ry, rz);
    double ref[] = {
         7.,  6., 29.,
        17., 17., 17.,
        17., 17., 17.,
        24., 25.,  8.,
        17., 17., 17.,
        17., 17., 17.,
        13., 12., 23.,
        16., 19., 14.,
        17., 17., 17.,
        17., 17., 17.,
        17., 17.
    };
    for (int i = 0; i < BDOCK_SIMD_ALIGNMENT; ++i) {
        EXPECT_NEAR(ref[i], grads[i], 1e-15);
    }
}

}
#endif
