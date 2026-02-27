/*
* Unified SIMD dispatch header.
* Includes the appropriate backend (AVX2 or NEON) and defines BDOCK_HAS_SIMD
* so that calling code can use a single guard.
*/

#pragma once

#if ENABLE_SIMD_AVX2
#include "bytedock/simd/avx2_def.h"
#define BDOCK_HAS_SIMD 1
#elif ENABLE_SIMD_NEON
#include "bytedock/simd/neon_def.h"
#define BDOCK_HAS_SIMD 1
#endif
