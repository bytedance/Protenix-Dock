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

#include <string>
#include <vector>

#include "bytedock/lib/utility.h"

namespace bytedock {

size_t get_sequential_id(const std::string& category);

enum workflow_step {
    PARSE_RECEPTOR_DATA_STEP,
    CONSTRUCT_RECEPTOR_CACHE_STEP,
    RESTORE_RECEPTOR_CACHE_STEP,
    SAVE_RECEPTOR_CACHE_STEP,
    PARSE_LIGAND_DATA_STEP,
    SEARCH_BINDING_POSE_STEP,
    REFINE_BINDING_POSE_STEP,
    EVALUATE_FAST_SCORE_STEP,
    EVALUATE_AFFINITY_SCORE_STEP,
    WORKFLOW_STEP_COUNT  // Only used as the upper limit
};

const char* to_string(workflow_step step);

class perf_counter {
public:
    static perf_counter& singleton();
    DISABLE_COPY_AND_ASSIGN(perf_counter);

    // Reset timer value of each step to zero
    void reset();

    // Return average timer values of all steps with unit "us"
    std::vector<size_t> get_time_cost() const;

private:
    perf_counter() {}
};

class step_timer {
public:
    explicit step_timer(workflow_step step);
    ~step_timer();

private:
    size_t get_timestamp_in_us();

    workflow_step step_;
    size_t start_;
};

class cache_reporter {
public:
    cache_reporter() : hit_(0), miss_(0) {}

    void hit() { hit_++; }
    void miss() { miss_++; }

    void tell(size_t *hit, size_t *miss) const {
        *hit += hit_;
        *miss += miss_;
    }

    void reset() {
        hit_ = 0;
        miss_ = 0;
    }

private:
    size_t hit_;
    size_t miss_;
};

class cache_counter {
public:
    static cache_counter& singleton();
    DISABLE_COPY_AND_ASSIGN(cache_counter);

    cache_reporter* get_or_create();
    double summary(size_t* hit, size_t* miss, bool clear = true);

private:
    cache_counter() {}
};

class nposes_manager {
public:
    static nposes_manager& singleton();
    bool insert(const std::string& key, size_t value);
    bool get(const std::string& key, size_t* value);
    void erase(const std::string& key);
};

}
