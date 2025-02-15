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

#include "bytedock/core/factory.h"
#include "bytedock/core/grid.h"

namespace bytedock {

// No thread safety!
void setup_global_logging(const std::string& path, int verbosity);

class ReusableEngine {
public:
    explicit ReusableEngine(const std::string& receptor_file,
                            const std::string& sf_file,
                            int cpu_nthreads);

    void require_box();
    void set_box(double center_x, double center_y, double center_z,
                 double size_x, double size_y, double size_z);

    void require_maps();
    void generate_maps(double grid_spacing, double max_size = 0.);
    void drop_maps();
    void dump_maps(const std::string& cache_dir);
    void load_maps(const std::string& cache_dir);

    void evaluate(const std::string& ligand_index_file,
                  const std::string& output_file);
    void optimize(const std::string& ligand_index_file,
                  const std::string& output_dir,
                  int max_niters, double slope);
    void search(const std::string& ligand_index_file,
                const std::string& output_dir,
                int seed,
                int exhaustiveness,
                int max_nsteps,
                int relax_nsteps,
                int num_modes,
                double min_rmsd,
                double mc_mmenergy_threshold,
                double slope);

private:
    const size_t nreaders_;
    const size_t nwriters_;
    const size_t naggregators_;
    size_t nworkers_;

    std::vector<double> box_;  // center_{x,y,z}, size_{x,y,z}
    box_config bc_;

    scoring_function_factory sf_mgr_;
    std::shared_ptr<torsional_receptor> receptor_;
    std::shared_ptr<receptor_cache> cache_;
};

}
