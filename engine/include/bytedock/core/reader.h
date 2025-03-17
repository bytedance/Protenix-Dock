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

#include "bytedock/core/factory.h"
#include "bytedock/core/task.h"
#include "bytedock/lib/queue.h"

namespace bytedock {

class ligand_pose_reader {
public:
    ligand_pose_reader(
        const scoring_function_factory& sf_manager, bool include_bscore,
        const docking_task& task_templ, blocking_queue<std::string>& file_queue,
        size_t max_ntasks = 0  // Zero means no limitation
    ): sf_mgr_(sf_manager), include_bscore_(include_bscore),
       task_templ_(task_templ), file_queue_(file_queue), max_ntasks_(max_ntasks) {}

    void fill(blocking_queue<name_and_task>& parsed_queue);

private:
    void handle_json_file(const std::string& path,
                          blocking_queue<name_and_task>& parsed_queue);

    const scoring_function_factory& sf_mgr_;
    const bool include_bscore_;
    const docking_task& task_templ_;
    const size_t max_ntasks_;
    blocking_queue<std::string>& file_queue_;
};

}
