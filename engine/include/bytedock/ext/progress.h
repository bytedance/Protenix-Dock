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
#include <mutex>

#include <boost/progress.hpp>

#include "bytedock/lib/utility.h"

namespace bytedock {

// A thread-safe wrapper of boost::timer::progress_display
class progress_bar {
public:
    progress_bar(size_t expected_count);
    DISABLE_COPY_AND_ASSIGN(progress_bar);

    size_t operator++();
    size_t operator+=(size_t increment);

private:
	std::unique_ptr<boost::progress_display> pd_;
    mutable std::mutex mtx_;
};

}
