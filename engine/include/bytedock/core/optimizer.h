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

#include "bytedock/core/interaction.h"
#include "bytedock/core/task.h"
#include "bytedock/lib/queue.h"

namespace bytedock {

class optimization_step {
public:
    // Update out.{ligand_xyz, torsions, energy}
    virtual bool apply(binding_system_interactions& model, optimized_result& out) = 0;
};

// `Serial` means it only uses 1 CPU thread
class serial_optimizer {
public:
    serial_optimizer(
        optimization_step& step, blocking_queue<name_and_task>& in_queue
    ): step_(step), task_queue_(in_queue) {}

    void fill(blocking_queue<name_and_task>& out_queue);

private:
    // Return `false` if it failed to converge
    bool minimize_energy(binding_input& in, optimized_result& out);

    optimization_step& step_;
    blocking_queue<name_and_task>& task_queue_;
};

// bytedock.docking.optimizer.lbfgsgl.LBFGSGL
class lbfgs_step : public optimization_step {
public:
    lbfgs_step(
        size_t max_niters, size_t max_nevals = 0
    ) : max_iter_(max_niters),
        max_eval_(max_nevals > 0 ? max_nevals : max_niters * 5 / 4) {}

    // tolerance_type=max & line_search_fn=strong_wolfe
    bool apply(binding_system_interactions& model, optimized_result& out) override;

private:
    const size_t max_iter_;
    const size_t max_eval_;
    const param_t lr_ = 1.;
    const param_t tolerance_grad_ = 0.5;
    const param_t tolerance_change_ = 1e-4;
    const size_t history_size_ = 100;
};

}
