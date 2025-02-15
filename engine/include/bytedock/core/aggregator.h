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

#include <unordered_map>

#include "bytedock/core/task.h"
#include "bytedock/lib/queue.h"

namespace bytedock {

class pose_aggregator {
public:
    pose_aggregator(
        blocking_queue<name_and_task>& pose_queue
    ): pose_queue_(pose_queue) {}

    void fill(blocking_queue<name_and_batch>& ligand_queue);

private:
    blocking_queue<name_and_task>& pose_queue_;

    mutable std::mutex mtx_;
    std::unordered_map<std::string, pose_batch> store_;
    std::unordered_map<std::string, size_t> remained_;
};

}
