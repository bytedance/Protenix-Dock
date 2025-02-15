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

#include "bytedock/ext/counter.h"

#include <atomic>
#include <mutex>
#include <thread>
#include <unordered_map>

#include <time.h>
#include <sys/time.h>

namespace bytedock {

static std::mutex mtx;
static std::unordered_map<std::string, size_t> id_counters;

size_t get_sequential_id(const std::string& category) {
    std::lock_guard<std::mutex> lock(mtx);
    if (id_counters.count(category) == 0) {
        id_counters.emplace(category, 0);
    }
    return ++(id_counters[category]);
}

const char* to_string(workflow_step step) {
    switch (step) {
        case PARSE_RECEPTOR_DATA_STEP: return "parse_receptor_data_step";
        case CONSTRUCT_RECEPTOR_CACHE_STEP: return "construct_receptor_cache_step";
        case RESTORE_RECEPTOR_CACHE_STEP: return "restore_receptor_cache_step";
        case SAVE_RECEPTOR_CACHE_STEP: return "save_receptor_cache_step";
        case PARSE_LIGAND_DATA_STEP:   return "parse_ligand_data_step";
        case SEARCH_BINDING_POSE_STEP: return "search_binding_pose_step";
        case REFINE_BINDING_POSE_STEP: return "refine_binding_pose_step";
        case EVALUATE_FAST_SCORE_STEP: return "evaluate_fast_score_step";
        case EVALUATE_AFFINITY_SCORE_STEP: return "evaluate_affinity_score_step";
        default:                       return "unknown_step";
    };
}

using counter_t = std::atomic<size_t>;
std::vector<counter_t> time_counters(WORKFLOW_STEP_COUNT);

perf_counter& perf_counter::singleton() {
    static perf_counter instance;
    return instance;
}

void perf_counter::reset() {
    for (size_t i = 0; i < WORKFLOW_STEP_COUNT; i++) {
        time_counters[i].store(0);
    }
}

std::vector<size_t> perf_counter::get_time_cost() const {
    std::vector<size_t> captured(WORKFLOW_STEP_COUNT);
    for (size_t i = 0; i < WORKFLOW_STEP_COUNT; i++) {
        captured[i] = time_counters[i].load();
    }
    return captured;
}

step_timer::step_timer(workflow_step step): step_(step) {
    start_ = get_timestamp_in_us();
}

step_timer::~step_timer() {
    size_t elapsed_us = get_timestamp_in_us() - start_;
    time_counters[step_] += elapsed_us;
}

size_t step_timer::get_timestamp_in_us() {
    timeval time;
    if (gettimeofday(&time, NULL)){
        return 0;
    }
    return time.tv_sec * 1000000 + time.tv_usec;
}

static std::unordered_map<std::thread::id, cache_reporter> reporters;

cache_counter& cache_counter::singleton() {
    static cache_counter instance;
    return instance;
}

cache_reporter* cache_counter::get_or_create() {
    std::lock_guard<std::mutex> lock(mtx);
    auto tid = std::this_thread::get_id();
    if (reporters.find(tid) == reporters.end()) {
        reporters[tid] = cache_reporter();
    }
    return &reporters[tid];
}

double cache_counter::summary(size_t* hit, size_t* miss, bool clear) {
    std::lock_guard<std::mutex> lock(mtx);
    size_t tmp_hit = 0, tmp_miss = 0;
    for (const auto& iter : reporters) iter.second.tell(&tmp_hit, &tmp_miss);
    if (clear) {
        for (auto& iter : reporters) iter.second.reset();
    }
    if (hit != nullptr) *hit = tmp_hit;
    if (miss != nullptr) *miss = tmp_miss;
    size_t total = tmp_hit + tmp_miss;
    return total > 0 ? (tmp_hit * 1. / total) : 1.;
}

nposes_manager& nposes_manager::singleton() {
    static nposes_manager instance;
    return instance;
}

static std::mutex npose_mtx;
static std::unordered_map<std::string, size_t> npose_store;

bool nposes_manager::insert(const std::string& key, size_t value) {
    std::lock_guard<std::mutex> lock(npose_mtx);
    auto iter = npose_store.find(key);
    if (iter == npose_store.end()) {
        npose_store[key] = value;
        return true;
    } else {
        return false;
    }
}

bool nposes_manager::get(const std::string& key, size_t* value) {
    std::lock_guard<std::mutex> lock(npose_mtx);
    auto iter = npose_store.find(key);
    if (iter == npose_store.end()) {
        return false;
    } else {
        *value = iter->second;
        return true;
    }
}

void nposes_manager::erase(const std::string& key) {
    std::lock_guard<std::mutex> lock(npose_mtx);
    npose_store.erase(key);
}

}
