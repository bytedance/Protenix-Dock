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

#include <string>

#include "bytedock/core/task.h"
#include "bytedock/lib/queue.h"

namespace bytedock {

class single_pose_reader {
public:
    single_pose_reader(
        const docking_task& origin, size_t nposes_to_generate,
        blocking_queue<std::string>& file_queue
    ): origin_(origin), gen_factor_(nposes_to_generate), file_queue_(file_queue) {}

    void fill(blocking_queue<name_and_task>& parsed_queue);

private:
    void handle_json_file(const std::string& path,
                          blocking_queue<name_and_task>& parsed_queue);

    const docking_task origin_;
    size_t gen_factor_;
    blocking_queue<std::string>& file_queue_;
};

// Read multiple poses per ligand
class multi_pose_reader {
public:
    multi_pose_reader(
        const docking_task& origin,
        blocking_queue<std::string>& file_queue
    ): origin_(origin), file_queue_(file_queue) {}

    void fill(blocking_queue<name_and_task>& parsed_queue);

private:
    void handle_json_file(const std::string& path,
                          blocking_queue<name_and_task>& parsed_queue);

    const docking_task origin_;
    blocking_queue<std::string>& file_queue_;
};

}
