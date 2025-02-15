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

#include <iostream>
#include <memory>
#include <vector>

namespace bytedock {

#if ENABLE_VETOS_ACCESS
inline bool is_tos_path(const std::string& path) {
    return path.compare(0, 4, "tos:") == 0;
}
#endif

uintmax_t get_file_length(const std::string& path);

std::shared_ptr<std::istream> open_for_read(const std::string& path);

std::shared_ptr<std::ostream> open_for_write(const std::string& path);

/**
 * The only difference between "write" and "append" is whether characters in
 * buffer are auto flushed before destructed.
 */
std::unique_ptr<std::ostream> open_for_append(const std::string& path);

bool is_file(const std::string& path);

void create_directory(const std::string& path);

std::string get_executable_path();

}
