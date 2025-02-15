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

#include "bytedock/core/scorer.h"

namespace bytedock {

class scoring_function_factory {
public:
    explicit scoring_function_factory(const std::string& config_path);
    explicit scoring_function_factory() : scoring_function_factory("") {}

    const ScoringFunctionConfig& get_config() const { return sf_config_; }
    const root_scorer& get_pose_selection() const { return *pose_scorer_; }
    const root_scorer& get_affinity_ranking() const { return *affinity_scorer_; }

private:
    ScoringFunctionConfig sf_config_;
    std::shared_ptr<root_scorer> pose_scorer_;
    std::shared_ptr<root_scorer> affinity_scorer_;
};

}
