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

#include <gtest/gtest.h>

#include "bytedock/lib/dtype.h"

namespace bytedock {

inline void check_3d(const param_t* ref, const param_t* pred, const param_t atol) {
    EXPECT_NEAR(ref[0], pred[0], atol);
    EXPECT_NEAR(ref[1], pred[1], atol);
    EXPECT_NEAR(ref[2], pred[2], atol);
}

}
