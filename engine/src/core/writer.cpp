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

#include "bytedock/core/writer.h"

#include <iomanip>

#include "bytedock/ext/logging.h"
#include "bytedock/lib/printer.h"

namespace bytedock {

inline void fix_backslashes(std::ostream& out, const std::string& val) {
    const char* buffer = val.data();
    size_t start = 0, i;
    for (i = 0; i < val.size(); ++i) {
        if (buffer[i] == '\\') {
            out.write(buffer + start, i - start + 1);
            out.write(buffer + i, 1);  // Double backslash
            start = i + 1;
        }
    }
    out.write(buffer + start, i - start);  // i == val.size()
}

std::string pose_writer::consume(const std::string& name, pose_batch& batch) {
    auto& candidates = batch.candidates;
    auto output_file = output_dir_ / (name + "_out.json");
    auto path = output_file.string();
    auto content = open_for_write(path);
    *content << "{" << std::endl;
    const auto& smiles = batch.ligand->get_smiles();
    if (!smiles.empty()) {
        *content << "  \"mapped_smiles\": \"";
        fix_backslashes(*content, smiles);
        *content << "\"," << std::endl;
    }
    *content << "  \"best_pose\": {" << std::endl;
    *content << "    \"index\": " << batch.best_id << "," << std::endl;
    *content << "    \"bscore\": " << batch.bscore << std::endl;
    *content << "  }," << std::endl;
    *content << "  \"poses\": [" << std::endl;
    for (size_t i = 0; i < candidates.size(); ++i) {
        const auto& result = candidates[i];
        auto aligned = batch.receptor->restore_bond_orders(result.torsions);  // Move
        *content << "    {" << std::endl;
        *content << "      \"offset\": " << result.offset << "," << std::endl;
        *content << "      \"energy\": " << result.energy << "," << std::endl;
        *content << "      \"pscore\": " << result.pscore << "," << std::endl;
        *content << "      \"nevals\": " << result.nevals << "," << std::endl;
        *content << "      \"receptor\": {" << std::endl;
        *content << "        \"torsions\": [" << std::endl
                 << format<param_t, 10>(aligned, "          ")
                 << "        ]" << std::endl;
        *content << "      }," << std::endl;
        *content << "      \"ligand\": {" << std::endl;
        *content << "        \"xyz\": " << to_string(result.ligand_xyz, "        ")
                 << std::endl;
        *content << "      }" << std::endl;
        *content << "    }";
        if (i < candidates.size() - 1) *content << ",";
        *content << std::endl;
    }
    *content << "  ]" << std::endl;
    *content << "}" << std::endl;
    return path;
}

void pose_writer::fill(blocking_queue<std::string>& file_queue) {
    auto pair = std::move(batch_queue_.pop());
    while (!batch_queue_.is_eoq(pair)) {
        nposes_manager::singleton().erase(pair.first);
        if (preprocess(pair.second)) {
            LOG_DEBUG << "Tasks with prefix [" << pair.first << "] has finished.";
            file_queue.push(consume(pair.first, pair.second));
        } else {
            LOG_DEBUG << "Tasks with prefix [" << pair.first << "] has failed.";
        }
        pair = std::move(batch_queue_.pop());
    }
    batch_queue_.close();  // Put back the EOF flag
    file_queue.close();
}

inline param_t calculate_rmsd_upper_bound(const molecule_pose& lhs,
                                          const molecule_pose& rhs) {
    param_t sum = 0_r;
    for (size_t i = 0; i < lhs.size(); ++i) {
        sum += get_distance_square_3d(lhs[i].xyz, rhs[i].xyz);
    }
    return std::sqrt(sum / lhs.size());
}

static void find_closest_pose(
    const std::vector<optimized_result>& candidates, const optimized_result& target,
    size_t& offset, param_t& rmsd
) {
    // It is ensured that there is at least 1 candidate
    param_t tmp;
    for (size_t i = 0; i < candidates.size(); ++i) {
        auto& cdd = candidates[i];
        tmp = calculate_rmsd_upper_bound(cdd.ligand_xyz, target.ligand_xyz);
        if (i == 0 || tmp < rmsd) {
            offset = i;
            rmsd = tmp;
        }
    }
}

bool pose_cluster::preprocess(pose_batch& aggregated) {
    auto& candidates = aggregated.candidates;  // At least 1 pose

    // Cluster conformer candidates
    std::vector<optimized_result> sorted;
    sorted.reserve(num_modes_);
    size_t offset;
    param_t rmsd;
    bool updated;
    for (size_t i = 0; i < candidates.size(); ++i) {
        optimized_result& cdd = candidates[i];
        if (cdd.nevals == 0) continue;  // Failure mark from searcher
        if (sorted.size() == 0) {
            sorted.push_back({});
            sorted.back() = std::move(cdd);
            continue;
        }
        find_closest_pose(sorted, cdd, offset, rmsd);
        updated = false;
        if (rmsd < min_rmsd_) {
            if (cdd.pscore < sorted[offset].pscore) {
                sorted[offset] = std::move(cdd);
                updated = true;
            }
        } else if (sorted.size() < num_modes_) {
            sorted.push_back({});
            sorted.back() = std::move(cdd);
            updated = true;
        } else if (cdd.pscore < sorted.back().pscore) {
            sorted.back() = std::move(cdd);
            updated = true;
        }
        if (updated) {
            std::sort(sorted.begin(), sorted.end(), [](auto& lhs, auto& rhs) {
                return lhs.pscore < rhs.pscore;
            });
        }
    }
    if (sorted.size() == 0) return false;
    candidates.swap(sorted);

    // Calculate affinity score
    aggregated.best_id = 0;
    {
        step_timer t(EVALUATE_AFFINITY_SCORE_STEP);
        auto& bb = sf_mgr_.get_affinity_ranking();
        auto& best_pose = candidates[0];
        aggregated.bscore = bb.combine(
            best_pose.receptor_xyz, aggregated.receptor->get_ffdata(),
            best_pose.ligand_xyz, aggregated.ligand->get_ffdata()
        );
    }
    return true;
}

bool pose_ranker::preprocess(pose_batch& aggregated) {
    auto& candidates = aggregated.candidates;  // At least 1 pose
    auto& pp = sf_mgr_.get_pose_selection();
    auto& receptor_ffdata = aggregated.receptor->get_ffdata();
    auto& ligand_ffdata = aggregated.ligand->get_ffdata();

    size_t best_id;
    instant_cache tmp, memo;
    for (size_t i = 0; i < candidates.size(); ++i) {
        auto& cdd = candidates[i];
        step_timer t(EVALUATE_FAST_SCORE_STEP);
        memo = pp.precalculate(cdd.receptor_xyz, receptor_ffdata,
                               cdd.ligand_xyz, ligand_ffdata);  // Swap
        cdd.pscore = pp.combine(cdd.receptor_xyz, receptor_ffdata,
                                cdd.ligand_xyz, ligand_ffdata, memo);
        if (i == 0 || cdd.pscore < candidates[best_id].pscore) {
            best_id = i;
            tmp = std::move(memo);
        }
    }
    aggregated.best_id = best_id;
    {
        step_timer t(EVALUATE_AFFINITY_SCORE_STEP);
        auto& bb = sf_mgr_.get_affinity_ranking();
        auto& best_pose = candidates[best_id];
        aggregated.bscore = bb.combine(best_pose.receptor_xyz, receptor_ffdata,
                                       best_pose.ligand_xyz, ligand_ffdata, tmp);
    }
    return true;
}

void report_writer::collect(blocking_queue<name_and_report>& in_queue) {
    auto content = open_for_write(output_file_.string());
    *content << "ligand_file,pose_id,pscore,bscore" << std::endl;
    auto pair = std::move(in_queue.pop());
    size_t num_poses = 0;
    while (!in_queue.is_eoq(pair)) {
        auto& task = pair.second;
        *content << remove_pose_suffix_value(pair.first) << "," << task.feed.pose_id << ","
#if ENABLE_DOUBLE_PRECISION
                 << std::setprecision(16) << task.fetch.pscore << ","
                 << std::setprecision(16) << task.fetch.bscore << std::endl;
#else
                 << std::setprecision(7) << task.fetch.pscore << ","
                 << std::setprecision(7) << task.fetch.bscore << std::endl;
#endif
        num_poses += 1;
        pair = std::move(in_queue.pop());
    }
    LOG_INFO << "Number of successful evaluations is: " << num_poses;
    in_queue.close();  // Put back the EOF flag
}

}
