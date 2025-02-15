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

#include "bytedock/core/grid.h"
#include "bytedock/lib/error.h"

namespace bytedock {

struct binding_input {
    std::shared_ptr<torsional_receptor> receptor;
    std::shared_ptr<receptor_cache> cache;  // Could be empty!!!
    std::shared_ptr<free_ligand> ligand;
    index_t pose_id;  // Slice of field `ligand.poses_`
};

struct optimized_result {
    std::vector<param_t> torsions;  // Torisons for receptor leaves
    molecule_pose receptor_xyz;  // Cache coordinates calculated from torsions
    molecule_pose ligand_xyz;
    param_t energy;  // Molcular mechanics energy of both receptor and ligand
    param_t pscore;  // Pose score correlated with interaction affinity
    index_t nevals;  // Count of energy & force evaluations
    index_t offset;  // Walker index in searching or pose index in optimization
    std::string name;  // name of this pose
};

// Suffix is defined in file "src/core/reader.cpp"
inline std::string remove_pose_suffix_value(const std::string& value) {
    auto pos = value.find_last_of('/');
    if (pos == std::string::npos) return value;
    return value.substr(0, pos);
}

extern std::unordered_map<std::string,
    std::unordered_map<std::string, param_t>> global_energy_info;

// receptor_variable=torsion & ligand_varaible=xyz
struct docking_task {
    binding_input feed;
    optimized_result fetch;
};

struct scoring_result {
    param_t pscore;  // Pose score
    param_t bscore;  // Affinity score
};

struct evaluation_task {
    binding_input feed;
    scoring_result fetch;
};

struct pose_batch {
    std::shared_ptr<torsional_receptor> receptor;
    std::shared_ptr<free_ligand> ligand;
    index_t best_id;  // Whose pscore is the lowest
    param_t bscore;  // For the best pose
    std::vector<optimized_result> candidates;
};

typedef std::pair<std::string, docking_task> name_and_task;
typedef std::pair<std::string, evaluation_task> name_and_report;
typedef std::pair<std::string, pose_batch> name_and_batch;

}
