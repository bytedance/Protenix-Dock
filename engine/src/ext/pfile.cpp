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

#include "bytedock/ext/pfile.h"

#include <cassert>
#include <fstream>
#include <limits.h>
#include <unistd.h>
#include <utility>

#include <boost/filesystem/operations.hpp>

#if ENABLE_VETOS_ACCESS
#include "bytedock/ext/tos.h"
#endif
#include "bytedock/lib/error.h"

namespace bytedock {

namespace fs = boost::filesystem;

#if ENABLE_VETOS_ACCESS
inline std::pair<std::string, std::string> split_tos_path(const std::string& inp) {
    size_t offset = 4;  // Length of "tos:"
    while (inp[offset] == '/') offset++;
    size_t dot = inp.find("/", offset);
    if (dot == std::string::npos) {
        throw file_system_error("No bucket field found in path: " + inp);
    }
    return std::pair(inp.substr(offset, dot-offset), inp.substr(dot+1));
}
#endif

uintmax_t get_file_length(const std::string& path) {
#if ENABLE_VETOS_ACCESS
    if (is_tos_path(path)) {
        auto bucket_and_key = split_tos_path(path);
        int64_t length = toutiao_object_storage::singleton().get_file_length(
            bucket_and_key.first, bucket_and_key.second
        );
        assert(length > 0);
        return static_cast<uintmax_t>(length);
    } else {
#endif
        return fs::file_size(path);
#if ENABLE_VETOS_ACCESS
    }
#endif
}

std::shared_ptr<std::istream> open_for_read(const std::string& path) {
#if ENABLE_VETOS_ACCESS
    if (is_tos_path(path)) {
        auto bucket_and_key = split_tos_path(path);
        return toutiao_object_storage::singleton().open_for_read(
            bucket_and_key.first, bucket_and_key.second
        );
    } else {
#endif
        auto in = std::make_shared<std::ifstream>(path);
        if (!in->good()) {
            throw file_system_error("Failed to open file for read: " + path);
        }
        return std::static_pointer_cast<std::istream>(in);
#if ENABLE_VETOS_ACCESS
    }
#endif
}

std::shared_ptr<std::ostream> open_for_write(const std::string& path) {
#if ENABLE_VETOS_ACCESS
    if (is_tos_path(path)) {
        auto bucket_and_key = split_tos_path(path);
        return toutiao_object_storage::singleton().open_for_write(
            bucket_and_key.first, bucket_and_key.second
        );
    } else {
#endif
        auto out = std::make_shared<std::ofstream>(path);
        if (!out->good()) {
            throw file_system_error("Failed to open file for write: " + path);
        }
        return std::static_pointer_cast<std::ostream>(out);
#if ENABLE_VETOS_ACCESS
    }
#endif
}

std::unique_ptr<std::ostream> open_for_append(const std::string& path) {
#if ENABLE_VETOS_ACCESS
    if (is_tos_path(path)) {
        auto bucket_and_key = split_tos_path(path);
        return toutiao_object_storage::singleton().open_for_append(
            bucket_and_key.first, bucket_and_key.second
        );
    } else {
#endif
        auto raw = new std::ofstream(path, std::ios_base::app);
        return std::unique_ptr<std::ostream>(raw);
#if ENABLE_VETOS_ACCESS
    }
#endif
}

bool is_file(const std::string& path) {
#if ENABLE_VETOS_ACCESS
    if (is_tos_path(path)) {
        auto bucket_and_key = split_tos_path(path);
        return toutiao_object_storage::singleton().is_file(
            bucket_and_key.first, bucket_and_key.second
        );
    } else {
#endif
        return fs::exists(path);
#if ENABLE_VETOS_ACCESS
    }
#endif
}

void create_directory(const std::string& path) {
#if ENABLE_VETOS_ACCESS
    if (!is_tos_path(path)) {
#endif
        /**
         * If failed, it will throw an exception.
         * If existed, it returns false without exception.
         */
        fs::create_directory(path);
#if ENABLE_VETOS_ACCESS
    }
#endif
}

std::string get_executable_path() {
    char buffer[PATH_MAX];
    ssize_t count = readlink( "/proc/self/exe", buffer, PATH_MAX);
    if (count < 0 || count >= PATH_MAX) return "";
    buffer[count] = '\0';
    return std::string(buffer);
}

}
