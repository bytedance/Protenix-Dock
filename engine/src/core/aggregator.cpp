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

#include "bytedock/core/aggregator.h"

#include "bytedock/ext/logging.h"

namespace bytedock {

// Suffix is defined in file "src/core/reader.cpp"
inline void remove_pose_suffix(std::string& pose_name) {
    auto pos = pose_name.rfind('/');
    if (pos != std::string::npos) pose_name.resize(pos);
}

void pose_aggregator::fill(blocking_queue<name_and_batch>& ligand_queue) {
    auto pair = std::move(pose_queue_.pop());
    size_t nposes;
    while (!pose_queue_.is_eoq(pair)) {
        auto& name = pair.first;
        LOG_INFO << "Task [" << name << "] has finished.";
        remove_pose_suffix(name);
        if (!nposes_manager::singleton().get(name, &nposes)) {
            LOG_WARNING << "An unexpected ligand (no pose count declared) is found "
                        << "during aggregation: " << name;
            continue;
        }
        auto& task = pair.second;
        auto iter1 = remained_.find(name);
        if (iter1 == remained_.end()) {
            remained_[name] = nposes;
            iter1 = remained_.find(name);
            auto& batch = store_[name];  // Create a new entry
            batch.receptor = task.feed.receptor;
            batch.ligand = task.feed.ligand;
            batch.candidates.resize(nposes);
        }
        store_[name].candidates[task.feed.pose_id] = std::move(task.fetch);
        if (iter1->second > 1) {
            iter1->second -= 1;
        } else {  // All poses of the same ligand are collected
            auto iter2 = store_.find(name);
            ligand_queue.push({std::move(name), std::move(iter2->second)});
            remained_.erase(iter1);
            store_.erase(iter2);
        }
        pair = std::move(pose_queue_.pop());
    }
    pose_queue_.close();  // Put back the EOF flag
    ligand_queue.close();
}

}
