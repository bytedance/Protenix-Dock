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

#include "bytedock/simd/allocator.h"

#include <cstdlib>

namespace bytedock {

// For most Intel processors, the cache line size is 64 bytes
static const size_t kSimdAlignmentInBytes = 64;

void* simd_allocation_policy::alloc(std::size_t bytes) {
    bytes += kSimdAlignmentInBytes;  // Pad memory at the end to avoid false sharing
    return std::aligned_alloc(kSimdAlignmentInBytes, bytes);
}

void simd_allocation_policy::free(void* p) {
    if (p) std::free(p);
}

}
