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

#include "bytedock/ext/progress.h"

namespace bytedock {

progress_bar::progress_bar(size_t expected_count) {
    pd_.reset(new boost::progress_display(expected_count));
}

size_t progress_bar::operator++() {
    std::unique_lock<std::mutex> lock(mtx_);
    return ++(*pd_);
}

size_t progress_bar::operator+=(size_t increment) {
    std::unique_lock<std::mutex> lock(mtx_);
    return (*pd_) += increment;
}

}
