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

namespace bytedock {

#if ENABLE_DOUBLE_PRECISION
typedef double param_t;
#else
typedef float param_t;
#endif

// Refer to https://stackoverflow.com/a/36594669/28276974
using size_t = decltype(sizeof(int));
typedef size_t index_t;

constexpr param_t operator ""_r(long double value) {
    return (param_t)value;
}

constexpr param_t operator ""_r(unsigned long long int value) {
    return (param_t)value;
}

}
