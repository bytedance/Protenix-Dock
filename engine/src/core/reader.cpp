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

#include "bytedock/core/reader.h"

#include <boost/filesystem/path.hpp>

#include "bytedock/ext/pfile.h"
#include "bytedock/ext/logging.h"

namespace bytedock {

std::unordered_map<std::string,
    std::unordered_map<std::string, param_t>> global_energy_info = {};
inline std::string get_file_name(const std::string& path) {
    boost::filesystem::path inp(path);
    return inp.stem().string();
}

void ligand_pose_reader::handle_json_file(
    const std::string& path, blocking_queue<name_and_task>& parsed_queue
) {
    auto ligand = std::make_shared<free_ligand>();
    try {
        auto in = open_for_read(path);
        ligand->parse(in);
    } catch (std::exception& e) {
        LOG_WARNING << "Due to invalid format or content, skip ligand: " << path
                    << std::endl << " with reason: " << e.what();
        return;
    }

    // Precalculate coefficients for nonbonded terms in pscore
    auto pscorer = std::make_shared<root_scorer>(sf_mgr_.get_pose_selection());
    std::shared_ptr<root_scorer> bscorer;
    if (include_bscore_) bscorer = std::make_shared<root_scorer>(
        sf_mgr_.get_affinity_ranking()
    );

    // Used by `mc_prune_energy_threshold` in Monte-Carlo searching
    auto file_name = get_file_name(path);
    global_energy_info[file_name] = {
        {"minimum_pre_opt_energy", 1e8_r},
        {"minimum_accepted_energy", 1e8_r},
    };

    // Generate tasks per pose
    size_t ntasks_to_generate = MATH_PAIR_MAX(ligand->get_nposes(), max_ntasks_);
    nposes_manager::singleton().insert(file_name, ntasks_to_generate);
    auto prefix = file_name + "/pose-";
    for (size_t i = 0; i < ntasks_to_generate; ++i) {
        std::string name(prefix + std::to_string(i));
        docking_task task(task_templ_);  // Copy
        task.feed.ligand = ligand;
        task.feed.pose_id = i;
        task.feed.pscorer = pscorer;
        task.feed.bscorer = bscorer;
        parsed_queue.push({std::move(name), std::move(task)});
    }
}

void ligand_pose_reader::fill(blocking_queue<name_and_task>& parsed_queue) {
    std::string path = std::move(file_queue_.pop());
    while (!file_queue_.is_eoq(path)) {
        handle_json_file(path, parsed_queue);
        path = std::move(file_queue_.pop());
    }
    file_queue_.close();  // Put back the EOF flag
    parsed_queue.close();
}

}
