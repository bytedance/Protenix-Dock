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

#include "bytedock/core/entrance.h"

#include <iomanip>
#include <thread>

#include "bytedock/core/aggregator.h"
#include "bytedock/core/evaluator.h"
#include "bytedock/core/optimizer.h"
#include "bytedock/core/reader.h"
#include "bytedock/core/searcher.h"
#include "bytedock/core/source.h"
#include "bytedock/core/writer.h"
#include "bytedock/ext/logging.h"

namespace bytedock {

class global_logging_manager {
public:
    ~global_logging_manager() {
        if (sink_) disable_sink(sink_, false);
    }

    void setup(const std::string& path, int verbosity) {
        if (sink_) disable_sink(sink_, true);
        sink_ = enable_global_logging(path, verbosity);
    }

private:
    boost::shared_ptr<sink_t> sink_;
};

void setup_global_logging(const std::string& path, int verbosity) {
    static global_logging_manager glm;  // Ensured to be destroyed at exit
    glm.setup(path, verbosity);
}

static void tell_metrics() {
    auto metrics = perf_counter::singleton().get_time_cost();
    std::ostringstream oss;
    oss << "Performance report is:" << std::endl
        << std::setw(32) << "workflow_step" << " " << "time_cost(us)" << std::endl;
    for (size_t i = 0; i < WORKFLOW_STEP_COUNT; i++) {
        if (metrics[i] > 0) {
            auto step = static_cast<workflow_step>(i);
            oss << std::setw(32) << to_string(step) << " " << metrics[i] << std::endl;
        }
    }
    LOG_INFO << oss.str();
}

ReusableEngine::ReusableEngine(
    const std::string& receptor_file, const std::string& sf_file, int cpu_nthreads
) : nreaders_(1), nwriters_(1), naggregators_(1), sf_mgr_(sf_file) {
    receptor_.reset(new torsional_receptor);
    receptor_->parse(open_for_read(receptor_file));
    if (cpu_nthreads < 1) {
        nworkers_ = std::thread::hardware_concurrency();
        LOG_INFO << "CPU threads is auto set to: " << nworkers_;
    } else {
        nworkers_ = cpu_nthreads;
    }
}

void ReusableEngine::require_box() {
    if (box_.size() == 0) {
        throw failed_conf_error("Please call `set_box(...)` first!");
    }
}

void ReusableEngine::set_box(double center_x, double center_y, double center_z,
                             double size_x, double size_y, double size_z) {
    box_.push_back(center_x);
    box_.push_back(center_y);
    box_.push_back(center_z);
    box_.push_back(size_x);
    box_.push_back(size_y);
    box_.push_back(size_z);
}

void ReusableEngine::require_maps() {
    if (!cache_) throw failed_conf_error(
        "Please call either `generate_cache_maps(...)` or `load_maps(...)` first!"
    );
}

void ReusableEngine::generate_maps(double grid_spacing, double max_size) {
    require_box();
    double size_x = std::max(box_[3], max_size);
    double size_y = std::max(box_[4], max_size);
    double size_z = std::max(box_[5], max_size);
    bc_ = decide_grid_dims(box_[0], box_[1], box_[2],
                           size_x, size_y, size_z, grid_spacing);  // Copy
    LOG_INFO << "Generating receptor cache with grid settings:" << std::endl
             << "- Front-left-bottom point: ["
             << std::fixed << std::setprecision(3) << bc_.init_xyz[0] << ", "
             << std::fixed << std::setprecision(3) << bc_.init_xyz[1] << ", "
             << std::fixed << std::setprecision(3) << bc_.init_xyz[2] << "]"
             << std::endl
             << "- Sample point count per dimension: ["
             << bc_.npoints[0] << ", " << bc_.npoints[1] << ", "
             << bc_.npoints[2] << "]";
    cache_.reset(new receptor_cache(bc_));
    cache_->populate(*receptor_, nworkers_);
}

void ReusableEngine::drop_maps() {
    cache_.reset();
}

void ReusableEngine::dump_maps(const std::string& cache_dir) {
    require_maps();
    create_directory(cache_dir);
    cache_->save(cache_dir, nworkers_);
}

void ReusableEngine::load_maps(const std::string& cache_dir) {
    cache_.reset(new receptor_cache());
    cache_->restore(cache_dir, nworkers_);
    bc_ = cache_->get_box_config();  // Copy
    LOG_INFO << "Previous grid settings:" << std::endl
             << "- Front-left-bottom point: ["
             << std::fixed << std::setprecision(3) << bc_.init_xyz[0] << ", "
             << std::fixed << std::setprecision(3) << bc_.init_xyz[1] << ", "
             << std::fixed << std::setprecision(3) << bc_.init_xyz[2] << "]"
             << std::endl
             << "- Spacing: " << bc_.spacing << std::endl
             << "- Sample point count per dimension: ["
             << bc_.npoints[0] << ", " << bc_.npoints[1] << ", "
             << bc_.npoints[2] << "]";
}

void ReusableEngine::evaluate(const std::string& ligand_index_file,
                              const std::string& output_file) {
    perf_counter::singleton().reset();

    // Prepare root task
    docking_task origin;
    origin.feed.receptor = receptor_;

    // Extract entries from source
    const std::string eoq_file = "";
    blocking_queue<std::string> inf_queue(nreaders_*2, eoq_file, 1);
    std::thread loader([&ligand_index_file, &inf_queue]() {
        LOG_TELL_THREAD_ID();
        LOG_DEBUG << "Loader has started...";
        index_file_source source(ligand_index_file);
        source.fill(inf_queue);
        LOG_DEBUG << "Loader has ended.";
    });

    // Start readers
    const name_and_task eoq_task("", {});
    blocking_queue<name_and_task> parsed_queue(nworkers_*2, eoq_task, nreaders_);
    std::vector<std::thread> readers;
    for (size_t i = 0; i < nreaders_; i++) {
        readers.emplace_back([&origin, &inf_queue, &parsed_queue]() {
            LOG_TELL_THREAD_ID();
            LOG_DEBUG << "Reader has started...";
            multi_pose_reader one(origin, inf_queue);
            one.fill(parsed_queue);
            LOG_DEBUG << "Reader has ended.";
        });
    }

    // Launch multiple evaluators to consume tasks
    const name_and_report eoq_report("", {});
    blocking_queue<name_and_report> pose_queue(nwriters_*2, eoq_report, nworkers_);
    std::vector<std::thread> evaluators;
    for (size_t i = 0; i < nworkers_; i++) {
        evaluators.emplace_back([&]() {
            LOG_TELL_THREAD_ID();
            LOG_INFO << "Evaluator has started...";
            scoring_evaluator one(sf_mgr_, parsed_queue);
            one.fill(pose_queue);
            LOG_INFO << "Evaluator has ended.";
        });
    }

    // Dump scores and wait for each thread's end
    report_writer writer(output_file);
    writer.collect(pose_queue);
    loader.join();
    for (size_t i = 0; i < readers.size(); i++) readers[i].join();
    for (size_t i = 0; i < evaluators.size(); i++) evaluators[i].join();
    tell_metrics();
}

void ReusableEngine::optimize(const std::string& ligand_index_file,
                              const std::string& output_dir,
                              int max_niters,
                              double slope) {
    perf_counter::singleton().reset();

    // Prepare root task
    docking_task origin;
    origin.feed.receptor = receptor_;
    if (!cache_) {
        LOG_INFO << "Receptor cache is disabled in this run.";
    } else {
        cache_->set_slope(slope);
        origin.feed.cache = cache_;
    }

    // Extract entries from source
    const std::string eoq_file = "";
    blocking_queue<std::string> inf_queue(nreaders_*2, eoq_file, 1);
    std::thread loader([&ligand_index_file, &inf_queue]() {
        LOG_TELL_THREAD_ID();
        LOG_DEBUG << "Loader has started...";
        index_file_source source(ligand_index_file);
        source.fill(inf_queue);
        LOG_DEBUG << "Loader has ended.";
    });

    // Start readers
    const name_and_task eoq_task("", {});
    blocking_queue<name_and_task> parsed_queue(nworkers_*2, eoq_task, nreaders_);
    std::vector<std::thread> readers;
    for (size_t i = 0; i < nreaders_; i++) {
        readers.emplace_back([&origin, &inf_queue, &parsed_queue]() {
            LOG_TELL_THREAD_ID();
            LOG_DEBUG << "Reader has started...";
            multi_pose_reader one(origin, inf_queue);
            one.fill(parsed_queue);
            LOG_DEBUG << "Reader has ended.";
        });
    }

    // Launch multiple optimizers to consume tasks
    blocking_queue<name_and_task> pose_queue(naggregators_*2, eoq_task, nworkers_);
    std::vector<std::thread> optimizers;
    lbfgs_step algo(max_niters);
    for (size_t i = 0; i < nworkers_; i++) {
        optimizers.emplace_back([&]() {
            LOG_TELL_THREAD_ID();
            LOG_INFO << "Optimizer has started...";
            serial_optimizer one(algo, parsed_queue);
            one.fill(pose_queue);
            LOG_INFO << "Optimizer has ended.";
        });
    }

    // Collect poses of the same ligand
    const name_and_batch eoq_batch("", {});
    blocking_queue<name_and_batch> ligand_queue(nwriters_*2,
                                                eoq_batch,
                                                naggregators_);
    std::vector<std::thread> aggregators;
    for (size_t i = 0; i < naggregators_; i++) {
        aggregators.emplace_back([&ligand_queue, &pose_queue]() {
            LOG_TELL_THREAD_ID();
            LOG_DEBUG << "Aggregator has started...";
            pose_aggregator one(pose_queue);
            one.fill(ligand_queue);
            LOG_DEBUG << "Aggregator has ended.";
        });
    }

    // Dump optimized poses to files
    create_directory(output_dir);
    blocking_queue<std::string> outf_queue(2, eoq_file, nwriters_);
    std::vector<std::thread> writers;
    for (size_t i = 0; i < nwriters_; i++) {
        writers.emplace_back([&]() {
            LOG_TELL_THREAD_ID();
            LOG_DEBUG << "Writer has started...";
            pose_ranker one(sf_mgr_, output_dir, ligand_queue);
            one.fill(outf_queue);
            LOG_DEBUG << "Writer has ended.";
        });
    }

    // Tell count of finished poses
    std::string entry = outf_queue.pop();
    size_t num_files = 0;
    while (!outf_queue.is_eoq(entry)) {
        num_files += 1;
        entry = std::move(outf_queue.pop());
    }
    LOG_INFO << "Number of successful optimizations is: " << num_files;
    outf_queue.close();

    // Wait for each thread's end
    loader.join();
    for (size_t i = 0; i < readers.size(); i++) readers[i].join();
    for (size_t i = 0; i < optimizers.size(); i++) optimizers[i].join();
    for (size_t i = 0; i < aggregators.size(); i++) aggregators[i].join();
    for (size_t i = 0; i < writers.size(); i++) writers[i].join();

    // Summary speed & hit
    tell_metrics();
    if (origin.feed.cache) {
        size_t hit, miss;
        double hit_rate = cache_counter::singleton().summary(&hit, &miss);
        LOG_INFO << "Receptor cache hit rate: "
                 << std::fixed << std::setprecision(2) << (hit_rate * 100.) << " % "
                 << "(" << hit << " / " << (hit + miss) << ")";
    }
}

void ReusableEngine::search(const std::string& ligand_index_file,
                            const std::string& output_dir,
                            int seed,
                            int exhaustiveness,
                            int max_nsteps,
                            int relax_nsteps,
                            int num_modes,
                            double min_rmsd,
                            double mc_mmenergy_threshold,
                            double slope) {
    perf_counter::singleton().reset();

    // Prepare root task
    docking_task origin;
    origin.feed.receptor = receptor_;
    if (!cache_) {
        require_box();
        bc_ = describe_boundary(box_[0], box_[1], box_[2],  // In case cache is dropped
                                box_[3], box_[4], box_[5]);
        LOG_INFO << "Receptor cache is disabled in this run.";
    } else {
        cache_->set_slope(slope);
        origin.feed.cache = cache_;
    }

    // Extract entries from source
    const std::string eoq_file = "";
    blocking_queue<std::string> inf_queue(nreaders_*2, eoq_file, 1);
    std::thread loader([&ligand_index_file, &inf_queue]() {
        LOG_TELL_THREAD_ID();
        LOG_DEBUG << "Loader has started...";
        index_file_source source(ligand_index_file);
        source.fill(inf_queue);
        LOG_DEBUG << "Loader has ended.";
    });

    // Start readers
    const name_and_task eoq_task("", {});
    blocking_queue<name_and_task> parsed_queue(nworkers_*2, eoq_task, nreaders_);
    std::vector<std::thread> readers;
    for (size_t i = 0; i < nreaders_; i++) {
        readers.emplace_back([&origin, exhaustiveness, &inf_queue, &parsed_queue]() {
            LOG_TELL_THREAD_ID();
            LOG_DEBUG << "Reader has started...";
            single_pose_reader one(origin, exhaustiveness, inf_queue);
            one.fill(parsed_queue);
            LOG_DEBUG << "Reader has ended.";
        });
    }

    // Launch multiple searchers to consume tasks
    blocking_queue<name_and_task> pose_queue(naggregators_*2, eoq_task, nworkers_);
    std::vector<std::thread> searchers;
    for (size_t i = 0; i < nworkers_; i++) {
        searchers.emplace_back([&]() {
            LOG_TELL_THREAD_ID();
            LOG_INFO << "Searcher has started...";
            monte_carlo_searcher one(sf_mgr_, bc_, seed, max_nsteps, relax_nsteps,
                                     mc_mmenergy_threshold, parsed_queue);
            one.fill(pose_queue);
            LOG_INFO << "Searcher has ended.";
        });
    }

    // Collect poses of the same ligand
    const name_and_batch eoq_batch("", {});
    blocking_queue<name_and_batch> ligand_queue(nwriters_*2, eoq_batch, naggregators_);
    std::vector<std::thread> aggregators;
    for (size_t i = 0; i < naggregators_; i++) {
        aggregators.emplace_back([&ligand_queue, &pose_queue]() {
            LOG_TELL_THREAD_ID();
            LOG_DEBUG << "Aggregator has started...";
            pose_aggregator one(pose_queue);
            one.fill(ligand_queue);
            LOG_DEBUG << "Aggregator has ended.";
        });
    }

    // Dump optimized poses to files
    create_directory(output_dir);
    blocking_queue<std::string> outf_queue(2, eoq_file, nwriters_);
    std::vector<std::thread> writers;
    for (size_t i = 0; i < nwriters_; i++) {
        writers.emplace_back([&]() {
            LOG_TELL_THREAD_ID();
            LOG_DEBUG << "Writer has started...";
            pose_cluster one(sf_mgr_, num_modes, min_rmsd, output_dir, ligand_queue);
            one.fill(outf_queue);
            LOG_DEBUG << "Writer has ended.";
        });
    }

    // Tell count of finished poses
    std::string entry = outf_queue.pop();
    size_t num_files = 0;
    while (!outf_queue.is_eoq(entry)) {
        num_files += 1;
        entry = std::move(outf_queue.pop());
    }
    LOG_INFO << "Number of successful dockings is: " << num_files;
    outf_queue.close();

    // Wait for each thread's end
    loader.join();
    for (size_t i = 0; i < readers.size(); i++) readers[i].join();
    for (size_t i = 0; i < searchers.size(); i++) searchers[i].join();
    for (size_t i = 0; i < aggregators.size(); i++) aggregators[i].join();
    for (size_t i = 0; i < writers.size(); i++) writers[i].join();

    // Summary speed & hit
    tell_metrics();
    if (origin.feed.cache) {
        size_t hit, miss;
        double hit_rate = cache_counter::singleton().summary(&hit, &miss);
        LOG_INFO << "Receptor cache hit rate: "
                 << std::fixed << std::setprecision(2) << (hit_rate * 100.) << " % "
                 << "(" << hit << " / " << (hit + miss) << ")";
    }
}

}
