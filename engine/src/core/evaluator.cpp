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

#include "bytedock/core/evaluator.h"

#include "bytedock/ext/counter.h"

namespace bytedock {

void scoring_evaluator::calculate_scores(std::string& name, binding_input& in,
                                         blocking_queue<name_and_report>& out_queue) {
    evaluation_task task;
    instant_cache memo;
    auto& receptor_xyz = in.receptor->get_positions();
    auto& receptor_ffdata = in.receptor->get_ffdata();
    auto& ligand_xyz = in.ligand->get_pose(in.pose_id);
    auto& ligand_ffdata = in.ligand->get_ffdata();
    {
        step_timer t(EVALUATE_FAST_SCORE_STEP);
        auto& pp = sf_mgr_.get_pose_selection();
        memo = pp.precalculate(receptor_xyz, receptor_ffdata,
                               ligand_xyz, ligand_ffdata);  // Move
        task.fetch.pscore = pp.combine(receptor_xyz, receptor_ffdata,
                                       ligand_xyz, ligand_ffdata, memo);
    }
    {
        step_timer t(EVALUATE_AFFINITY_SCORE_STEP);
        auto& bb = sf_mgr_.get_affinity_ranking();
        task.fetch.bscore = bb.combine(receptor_xyz, receptor_ffdata,
                                       ligand_xyz, ligand_ffdata, memo);
    }
    task.feed = std::move(in);
    out_queue.push({std::move(name), std::move(task)});
}

void scoring_evaluator::fill(blocking_queue<name_and_report>& out_queue) {
    name_and_task pair = std::move(task_queue_.pop());
    while (!task_queue_.is_eoq(pair)) {
        calculate_scores(pair.first, pair.second.feed, out_queue);
        pair = std::move(task_queue_.pop());
    }
    task_queue_.close();  // Put back the EOF flag
    out_queue.close();
}

}
