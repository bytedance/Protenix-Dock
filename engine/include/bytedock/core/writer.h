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

#include <boost/filesystem/path.hpp>

#include "bytedock/core/factory.h"
#include "bytedock/core/task.h"
#include "bytedock/ext/pfile.h"
#include "bytedock/lib/queue.h"

namespace bytedock {

namespace fs = boost::filesystem;

// Both ligand and flexible leaves of receptor
class pose_writer {
public:
    pose_writer(
        const scoring_function_factory& sf_manager, const std::string& output_dir,
        blocking_queue<name_and_batch>& in_queue
    ) : sf_mgr_(sf_manager), output_dir_(output_dir), batch_queue_(in_queue) {}

    void fill(blocking_queue<std::string>& file_queue);

protected:
    virtual bool preprocess(pose_batch& aggregated) = 0;

    const scoring_function_factory& sf_mgr_;

private:
    std::string consume(const std::string& name, pose_batch& batch);

    const fs::path output_dir_;
    blocking_queue<name_and_batch>& batch_queue_;
};

class pose_cluster : public pose_writer {
public:
    pose_cluster(
        const scoring_function_factory& sf_manager, size_t num_modes, param_t min_rmsd,
        const std::string& output_dir, blocking_queue<name_and_batch>& in_queue
    ) : pose_writer(sf_manager, output_dir, in_queue),
        num_modes_(num_modes), min_rmsd_(min_rmsd) {}

protected:
    bool preprocess(pose_batch& aggregated) override;

private:
    const size_t num_modes_;
    const param_t min_rmsd_;
};

class pose_ranker : public pose_writer {
public:
    pose_ranker(
        const scoring_function_factory& sf_manager, const std::string& output_dir,
        blocking_queue<name_and_batch>& in_queue
    ) : pose_writer(sf_manager, output_dir, in_queue) {}

protected:
    bool preprocess(pose_batch& aggregated) override;
};

class report_writer {
public:
    report_writer(const std::string& output_file) : output_file_(output_file) {}
    void collect(blocking_queue<name_and_report>& in_queue);

private:
    const fs::path output_file_;
};

}
