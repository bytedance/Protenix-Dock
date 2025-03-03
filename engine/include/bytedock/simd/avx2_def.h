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

#if ENABLE_SIMD_AVX2

#define BDOCK_SIMD_ALIGNMENT 32  // Bytes

#include "bytedock/simd/avx2_double.h"
#include "bytedock/simd/avx2_float.h"

namespace bytedock {
#if ENABLE_DOUBLE_PRECISION
#define BDOCK_SIMD_REAL_WIDTH BDOCK_SIMD_DOUBLE_WIDTH
typedef simd_double simd_real;
#else
#define BDOCK_SIMD_REAL_WIDTH BDOCK_SIMD_FLOAT_WIDTH
typedef simd_float simd_real;
#endif
}

#endif
