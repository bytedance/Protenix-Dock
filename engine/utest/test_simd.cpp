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

#ifdef BDOCK_HAS_SIMD
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

TEST(SimdFloatTest, Reciprocal) {
    auto a = simd_rcp(simd_float(6.54321f));
    alignas(BDOCK_SIMD_ALIGNMENT) float buffer[BDOCK_SIMD_FLOAT_WIDTH];
    simd_store(buffer, a);
    for (int i = 0; i < BDOCK_SIMD_FLOAT_WIDTH; ++i) {
        EXPECT_NEAR(buffer[i], 0.15283, 1e-6);
    }
}

TEST(SimdDoubleTest, Reciprocal) {
    auto a = simd_rcp(simd_double(9.7654321));
    alignas(BDOCK_SIMD_ALIGNMENT) double buffer[BDOCK_SIMD_DOUBLE_WIDTH];
    simd_store(buffer, a);
    for (int i = 0; i < BDOCK_SIMD_DOUBLE_WIDTH; ++i) {
        EXPECT_NEAR(buffer[i], 0.10240202274305, 1e-14);
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
    alignas(BDOCK_SIMD_ALIGNMENT) float buffer[BDOCK_SIMD_FLOAT_WIDTH];
    float expected = 0;
    for (int i = 0; i < BDOCK_SIMD_FLOAT_WIDTH; ++i) {
        buffer[i] = static_cast<float>(i * 3 - 5);
        expected += buffer[i];
    }
    auto a = simd_reduce(simd_load(buffer));
    EXPECT_NEAR(a, expected, 1e-6);
}

TEST(SimdDoubleTest, SumElements) {
    alignas(BDOCK_SIMD_ALIGNMENT) double buffer[BDOCK_SIMD_DOUBLE_WIDTH];
    double expected = 0;
    for (int i = 0; i < BDOCK_SIMD_DOUBLE_WIDTH; ++i) {
        buffer[i] = static_cast<double>(i * 7 - 12);
        expected += buffer[i];
    }
    auto a = simd_reduce(simd_load(buffer));
    EXPECT_NEAR(a, expected, 1e-15);
}

TEST(SimdFloatTest, LoadXYZ) {
    static constexpr int N = 48;  // Enough for max offset * 3
    float coords[N];
    for (int i = 0; i < N; ++i) coords[i] = static_cast<float>(i);
    std::int32_t offsets[] = {9, 7, 6, 5, 3, 2, 1, 0};  // First WIDTH used
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
    static constexpr int N = 48;
    double coords[N];
    for (int i = 0; i < N; ++i) coords[i] = static_cast<double>(i);
    std::int32_t offsets[] = {8, 5, 4, 1};  // First WIDTH used
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
    static constexpr int N = 48;
    const float init_val = -12.f;
    float grads[N], ref[N];
    for (int i = 0; i < N; ++i) grads[i] = ref[i] = init_val;
    std::int32_t offsets[] = {9, 7, 6, 5, 3, 2, 1, 0};  // First WIDTH used
    alignas(BDOCK_SIMD_ALIGNMENT) float buffer[3 * BDOCK_SIMD_FLOAT_WIDTH];
    for (int i = 0; i < 3 * BDOCK_SIMD_FLOAT_WIDTH; ++i)
        buffer[i] = static_cast<float>(i + 1);
    simd_float rx = simd_load(buffer);
    simd_float ry = simd_load(buffer + BDOCK_SIMD_FLOAT_WIDTH);
    simd_float rz = simd_load(buffer + 2*BDOCK_SIMD_FLOAT_WIDTH);
    simd_incru_xyz(grads, offsets, rx, ry, rz);
    for (int i = 0; i < BDOCK_SIMD_FLOAT_WIDTH; ++i) {
        ref[3*offsets[i] + 0] += buffer[i];
        ref[3*offsets[i] + 1] += buffer[BDOCK_SIMD_FLOAT_WIDTH + i];
        ref[3*offsets[i] + 2] += buffer[2*BDOCK_SIMD_FLOAT_WIDTH + i];
    }
    for (int i = 0; i < N; ++i) {
        EXPECT_NEAR(ref[i], grads[i], 1e-6);
    }
}

TEST(SimdDoubleTest, IncrXYZ) {
    static constexpr int N = 48;
    const double init_val = 5.;
    double grads[N], ref[N];
    for (int i = 0; i < N; ++i) grads[i] = ref[i] = init_val;
    std::int32_t offsets[] = {8, 5, 4, 1};  // First WIDTH used
    alignas(BDOCK_SIMD_ALIGNMENT) double buffer[3 * BDOCK_SIMD_DOUBLE_WIDTH];
    for (int i = 0; i < 3 * BDOCK_SIMD_DOUBLE_WIDTH; ++i)
        buffer[i] = static_cast<double>(i * 2 - 3);
    simd_double rx = simd_load(buffer);
    simd_double ry = simd_load(buffer + BDOCK_SIMD_DOUBLE_WIDTH);
    simd_double rz = simd_load(buffer + 2*BDOCK_SIMD_DOUBLE_WIDTH);
    simd_incru_xyz(grads, offsets, rx, ry, rz);
    for (int i = 0; i < BDOCK_SIMD_DOUBLE_WIDTH; ++i) {
        ref[3*offsets[i] + 0] += buffer[i];
        ref[3*offsets[i] + 1] += buffer[BDOCK_SIMD_DOUBLE_WIDTH + i];
        ref[3*offsets[i] + 2] += buffer[2*BDOCK_SIMD_DOUBLE_WIDTH + i];
    }
    for (int i = 0; i < N; ++i) {
        EXPECT_NEAR(ref[i], grads[i], 1e-15);
    }
}

TEST(SimdFloatTest, DecrXYZ) {
    static constexpr int N = 48;
    const float init_val = -9.f;
    float grads[N], ref[N];
    for (int i = 0; i < N; ++i) grads[i] = ref[i] = init_val;
    std::int32_t offsets[] = {9, 8, 7, 5, 4, 3, 2, 1};  // First WIDTH used
    float buffer[3 * BDOCK_SIMD_FLOAT_WIDTH];
    for (int i = 0; i < 3 * BDOCK_SIMD_FLOAT_WIDTH; ++i)
        buffer[i] = static_cast<float>(i * 3 - 7);
    simd_float rx = simd_loadu(buffer);
    simd_float ry = simd_loadu(buffer + BDOCK_SIMD_FLOAT_WIDTH);
    simd_float rz = simd_loadu(buffer + 2*BDOCK_SIMD_FLOAT_WIDTH);
    simd_decru_xyz(grads, offsets, rx, ry, rz);
    for (int i = 0; i < BDOCK_SIMD_FLOAT_WIDTH; ++i) {
        ref[3*offsets[i] + 0] -= buffer[i];
        ref[3*offsets[i] + 1] -= buffer[BDOCK_SIMD_FLOAT_WIDTH + i];
        ref[3*offsets[i] + 2] -= buffer[2*BDOCK_SIMD_FLOAT_WIDTH + i];
    }
    for (int i = 0; i < N; ++i) {
        EXPECT_NEAR(ref[i], grads[i], 1e-6);
    }
}

TEST(SimdDoubleTest, DecrXYZ) {
    static constexpr int N = 48;
    const double init_val = 17.;
    double grads[N], ref[N];
    for (int i = 0; i < N; ++i) grads[i] = ref[i] = init_val;
    std::int32_t offsets[] = {7, 6, 3, 0};  // First WIDTH used
    alignas(BDOCK_SIMD_ALIGNMENT) double buffer[3 * BDOCK_SIMD_DOUBLE_WIDTH];
    for (int i = 0; i < 3 * BDOCK_SIMD_DOUBLE_WIDTH; ++i)
        buffer[i] = static_cast<double>(i * 4 - 1);
    simd_double rx = simd_loadu(buffer);
    simd_double ry = simd_loadu(buffer + BDOCK_SIMD_DOUBLE_WIDTH);
    simd_double rz = simd_loadu(buffer + 2*BDOCK_SIMD_DOUBLE_WIDTH);
    simd_decru_xyz(grads, offsets, rx, ry, rz);
    for (int i = 0; i < BDOCK_SIMD_DOUBLE_WIDTH; ++i) {
        ref[3*offsets[i] + 0] -= buffer[i];
        ref[3*offsets[i] + 1] -= buffer[BDOCK_SIMD_DOUBLE_WIDTH + i];
        ref[3*offsets[i] + 2] -= buffer[2*BDOCK_SIMD_DOUBLE_WIDTH + i];
    }
    for (int i = 0; i < N; ++i) {
        EXPECT_NEAR(ref[i], grads[i], 1e-15);
    }
}

}
#endif
