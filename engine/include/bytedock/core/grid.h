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

#include <unordered_set>

#include <boost/filesystem/path.hpp>

#include "bytedock/core/data.h"
#include "bytedock/ext/progress.h"

namespace bytedock {

namespace fs = boost::filesystem;

struct box_config {
    param_t init_xyz[3];  // Init point is the front-left-bottom corner of the box
    index_t npoints[3];
    param_t spacing;

    // Derived fields
    param_t oppo_xyz[3];  // Opposite point is the back-right-top corner of the box
    index_t npoints_minus_one[3];  // npoints_minus_one = npoints - 1
    param_t spacing_inv;  // (atom_xyz - init_xyz) * spacing_inv => flb_ijk + flb_xyz
};

inline bool check_atom_in_box(const atom_position& atom, const box_config& bc) {
    return atom.xyz[0] > bc.init_xyz[0] && atom.xyz[0] < bc.oppo_xyz[0] &&
           atom.xyz[1] > bc.init_xyz[1] && atom.xyz[1] < bc.oppo_xyz[1] &&
           atom.xyz[2] > bc.init_xyz[2] && atom.xyz[2] < bc.oppo_xyz[2];
}

box_config describe_boundary(double center_x, double center_y, double center_z,
                             double size_x, double size_y, double size_z);
box_config decide_grid_dims(double center_x, double center_y, double center_z,
                            double size_x, double size_y, double size_z,
                            double sapcing);

struct val_and_grad {
    param_t first;
    atom_position second;
};

// First-order 3D cache for all ligand atom types
class receptor_cache {
public:
    explicit receptor_cache();  // For restoring only
    receptor_cache(const box_config& bc);

    // If `nthreads` < 2, no thread will be used
    void populate(const torsional_receptor& receptor, size_t nthreads = 0);

    /**
     * If `false` is returned, the atom is not within the box.
     * Both `atom_xyz` and `xyz_gradient` is updated with delta instead of replacement!
     */
    bool get(const index_t atom_type, const param_t partial_charge,
             const atom_position& atom_xyz,
             param_t& out_energy, atom_position& xyz_gradient) const;

    // Interactions with these atoms need be calcualted externally
    const std::unordered_set<index_t>& get_excluded_atoms() const { return excluded_; }
    const box_config& get_box_config() const { return bc_; }
    void set_slope(param_t value) { slope_ = value; }

    void save(const fs::path& map_dir, size_t nthreads = 0);
    void restore(const fs::path& map_dir, size_t nthreads = 0);

private:
    struct position_index {
        index_t p000;
        index_t p100;
        index_t p010;
        index_t p110;
        index_t p001;
        index_t p101;
        index_t p011;
        index_t p111;
        param_t x;
        param_t y;
        param_t z;
        param_t mx;
        param_t my;
        param_t mz;
    };

    class partial_worker {
    public:
        partial_worker(receptor_cache& cache) : cache_(cache) {}

        void populate(const torsional_receptor& receptor, progress_bar& pbar,
                      const index_t start_i, const index_t end_i);

        void save(const fs::path& map_dir, progress_bar& pbar,
                  const index_t start_i, const index_t end_i);

        void restore(const fs::path& map_dir, progress_bar& pbar,
                     const index_t start_i, const index_t end_i);

    private:
        receptor_cache& cache_;
    };

private:
    index_t get_offset(const index_t i, const index_t j, const index_t k) const {
        return (i * bc_.npoints[1] + j) * bc_.npoints[2] + k;
    }

    param_t get(const std::vector<val_and_grad>& points, const position_index& pi,
                param_t* xyz_gradient) const;

private:
    box_config bc_;
    std::vector<val_and_grad> coul_grid_;
    std::vector<std::vector<val_and_grad> > vdw_grids_;
    std::unordered_set<index_t> excluded_;  // Movable atoms
    param_t slope_;
};

}
