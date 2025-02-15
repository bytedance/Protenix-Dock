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

#include <sstream>
#include <vector>

namespace bytedock {

#define NELEMENTS_PER_LINE 5

template<typename T, size_t NCOLS = NELEMENTS_PER_LINE>
std::string format(const std::vector<T> values, const std::string& prefix = "") {
    std::ostringstream oss;
    const std::string indent(prefix + "  ");
    size_t nrows = values.size() / NCOLS;
    size_t remained = values.size() % NCOLS;
    if (remained == 0 && nrows > 0) {
        nrows--;
        remained = NCOLS;
    }
    for (size_t i = 0; i < nrows; i++) {
        oss << indent;
        for (size_t j = 0; j < NCOLS; j++) {
            oss << values[i*NCOLS+j] << ", ";
        }
        oss << std::endl;
    }
    if (remained > 0) {
        nrows *= NCOLS;  // Used as `offset`
        oss << indent;
        oss << values[nrows];
        for (size_t j = 1; j < remained; j++) oss << ", " << values[nrows+j];
        oss << std::endl;
    }
    return oss.str();
}

}
