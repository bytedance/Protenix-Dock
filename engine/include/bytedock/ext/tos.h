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

#include "TosClientV2.h"

#include "bytedock/lib/utility.h"

namespace bytedock {

class toutiao_object_storage {
public:
    static toutiao_object_storage& singleton();
    ~toutiao_object_storage();

    int64_t get_file_length(const std::string& bucket, const std::string& key);
    std::shared_ptr<std::istream> open_for_read(const std::string& bucket,
                                                const std::string& key);
    std::shared_ptr<std::ostream> open_for_write(const std::string& bucket,
                                                 const std::string& key);
    std::unique_ptr<std::ostream> open_for_append(const std::string& bucket,
                                                  const std::string& key);
    void upload_file(const std::string& bucket, const std::string& keyconst,
                     const std::string& local_path, size_t num_threads = 1);
    bool is_file(const std::string& bucket, const std::string& key);

private:
    toutiao_object_storage();  // Forbid class inheritence as well
    DISABLE_COPY_AND_ASSIGN(toutiao_object_storage);

    std::unique_ptr<VolcengineTos::TosClientV2> client_;
};

}
