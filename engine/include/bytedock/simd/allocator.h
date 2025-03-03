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

#include <memory>
#include <vector>

namespace bytedock {

class simd_allocation_policy {
public:
    static void* alloc(std::size_t bytes);
    static void free(void* p);
};

template<class T>
class simd_allocator {
public:
    typedef T value_type;

    value_type* allocate(std::size_t n) {
        void* p = simd_allocation_policy::alloc(n * sizeof(T));
        if (p == nullptr) throw std::bad_alloc();
        return static_cast<value_type*>(p);
    }

    void deallocate(value_type* p, std::size_t n) {
        simd_allocation_policy::free(p);
    }

    template<class T2>
    bool operator==(const simd_allocator<T2>& /*unused*/) const {
        return true;
    }

    template<class T2>
    bool operator!=(const simd_allocator<T2>& rhs) const {
        return !(*this == rhs);
    }
};

template<typename T>
using simd_vector = std::vector<T, simd_allocator<T>>;

}
