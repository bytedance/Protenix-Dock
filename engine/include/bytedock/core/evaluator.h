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

#include "bytedock/core/factory.h"
#include "bytedock/core/task.h"
#include "bytedock/lib/queue.h"

namespace bytedock {

class scoring_evaluator {
public:
    explicit scoring_evaluator(
        const scoring_function_factory& sf_manager,
        blocking_queue<name_and_task>& in_queue
    ) : sf_mgr_(sf_manager), task_queue_(in_queue) {}
    void fill(blocking_queue<name_and_report>& out_queue);

private:
    void calculate_scores(std::string& name, binding_input& in,
                          blocking_queue<name_and_report>& out_queue);

    const scoring_function_factory& sf_mgr_;
    blocking_queue<name_and_task>& task_queue_;
};

}
