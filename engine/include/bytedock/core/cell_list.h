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

#include <algorithm>
#include <vector>

#include "bytedock/core/system.h"
#include "bytedock/lib/math.h"

namespace bytedock {

/**
 * Spatial cell list for O(N) neighbor queries within a cutoff distance.
 * Divides 3D space into cubic cells and assigns atoms to cells.
 * Queries return only atoms in neighboring cells (3x3x3 = 27 cells).
 */
class cell_list {
public:
    cell_list() : cutoff_(0), inv_cutoff_(0), nx_(0), ny_(0), nz_(0) {}

    /**
     * Build the cell list from atom positions.
     * @param positions  Array of atom positions
     * @param n_atoms    Number of atoms
     * @param cutoff     Cutoff distance — cells are this size
     */
    void build(const molecule_pose& positions, param_t cutoff) {
        cutoff_ = cutoff;
        cutoff_sq_ = cutoff * cutoff;
        inv_cutoff_ = 1.0 / cutoff;
        const size_t n_atoms = positions.size();

        if (n_atoms == 0) {
            nx_ = ny_ = nz_ = 0;
            cells_.clear();
            return;
        }

        // Find bounding box
        param_t min_x = positions[0].xyz[0], max_x = min_x;
        param_t min_y = positions[0].xyz[1], max_y = min_y;
        param_t min_z = positions[0].xyz[2], max_z = min_z;
        for (size_t i = 1; i < n_atoms; i++) {
            const param_t* p = positions[i].xyz;
            if (p[0] < min_x) min_x = p[0]; else if (p[0] > max_x) max_x = p[0];
            if (p[1] < min_y) min_y = p[1]; else if (p[1] > max_y) max_y = p[1];
            if (p[2] < min_z) min_z = p[2]; else if (p[2] > max_z) max_z = p[2];
        }

        origin_[0] = min_x - cutoff;
        origin_[1] = min_y - cutoff;
        origin_[2] = min_z - cutoff;

        nx_ = static_cast<index_t>((max_x - min_x) * inv_cutoff_) + 3;
        ny_ = static_cast<index_t>((max_y - min_y) * inv_cutoff_) + 3;
        nz_ = static_cast<index_t>((max_z - min_z) * inv_cutoff_) + 3;

        size_t n_cells = static_cast<size_t>(nx_) * ny_ * nz_;
        cells_.clear();
        cells_.resize(n_cells);

        // Assign atoms to cells
        for (size_t i = 0; i < n_atoms; i++) {
            index_t ci = cell_index(positions[i].xyz);
            if (ci >= 0 && static_cast<size_t>(ci) < n_cells) {
                cells_[ci].push_back(static_cast<index_t>(i));
            }
        }
    }

    /**
     * Build from a subset of atom positions (e.g., excluded atoms).
     */
    void build_subset(const molecule_pose& positions,
                      const std::vector<index_t>& atom_indices,
                      param_t cutoff) {
        cutoff_ = cutoff;
        cutoff_sq_ = cutoff * cutoff;
        inv_cutoff_ = 1.0 / cutoff;

        if (atom_indices.empty()) {
            nx_ = ny_ = nz_ = 0;
            cells_.clear();
            return;
        }

        // Find bounding box over subset
        const param_t* p0 = positions[atom_indices[0]].xyz;
        param_t min_x = p0[0], max_x = min_x;
        param_t min_y = p0[1], max_y = min_y;
        param_t min_z = p0[2], max_z = min_z;
        for (size_t i = 1; i < atom_indices.size(); i++) {
            const param_t* p = positions[atom_indices[i]].xyz;
            if (p[0] < min_x) min_x = p[0]; else if (p[0] > max_x) max_x = p[0];
            if (p[1] < min_y) min_y = p[1]; else if (p[1] > max_y) max_y = p[1];
            if (p[2] < min_z) min_z = p[2]; else if (p[2] > max_z) max_z = p[2];
        }

        origin_[0] = min_x - cutoff;
        origin_[1] = min_y - cutoff;
        origin_[2] = min_z - cutoff;

        nx_ = static_cast<index_t>((max_x - min_x) * inv_cutoff_) + 3;
        ny_ = static_cast<index_t>((max_y - min_y) * inv_cutoff_) + 3;
        nz_ = static_cast<index_t>((max_z - min_z) * inv_cutoff_) + 3;

        size_t n_cells = static_cast<size_t>(nx_) * ny_ * nz_;
        cells_.clear();
        cells_.resize(n_cells);

        for (size_t i = 0; i < atom_indices.size(); i++) {
            index_t idx = atom_indices[i];
            index_t ci = cell_index(positions[idx].xyz);
            if (ci >= 0 && static_cast<size_t>(ci) < n_cells) {
                cells_[ci].push_back(idx);
            }
        }
    }

    /**
     * Call visitor(atom_index) for each atom within cutoff of query_xyz.
     * Only considers atoms in the 27 neighboring cells.
     */
    template<typename Visitor>
    void for_each_neighbor(const param_t* query_xyz, Visitor&& visitor) const {
        if (nx_ == 0) return;

        index_t cx = static_cast<index_t>((query_xyz[0] - origin_[0]) * inv_cutoff_);
        index_t cy = static_cast<index_t>((query_xyz[1] - origin_[1]) * inv_cutoff_);
        index_t cz = static_cast<index_t>((query_xyz[2] - origin_[2]) * inv_cutoff_);

        index_t x0 = (cx > 0) ? cx - 1 : 0;
        index_t y0 = (cy > 0) ? cy - 1 : 0;
        index_t z0 = (cz > 0) ? cz - 1 : 0;
        index_t x1 = (cx + 1 < nx_) ? cx + 1 : nx_ - 1;
        index_t y1 = (cy + 1 < ny_) ? cy + 1 : ny_ - 1;
        index_t z1 = (cz + 1 < nz_) ? cz + 1 : nz_ - 1;

        for (index_t ix = x0; ix <= x1; ix++) {
            for (index_t iy = y0; iy <= y1; iy++) {
                for (index_t iz = z0; iz <= z1; iz++) {
                    size_t cell_idx = static_cast<size_t>(ix) * ny_ * nz_
                                    + static_cast<size_t>(iy) * nz_ + iz;
                    const auto& cell = cells_[cell_idx];
                    for (index_t atom_idx : cell) {
                        visitor(atom_idx);
                    }
                }
            }
        }
    }

    param_t cutoff_sq() const { return cutoff_sq_; }
    bool empty() const { return nx_ == 0; }

private:
    index_t cell_index(const param_t* xyz) const {
        index_t ix = static_cast<index_t>((xyz[0] - origin_[0]) * inv_cutoff_);
        index_t iy = static_cast<index_t>((xyz[1] - origin_[1]) * inv_cutoff_);
        index_t iz = static_cast<index_t>((xyz[2] - origin_[2]) * inv_cutoff_);
        if (ix < 0 || ix >= nx_ || iy < 0 || iy >= ny_ || iz < 0 || iz >= nz_) {
            return -1;
        }
        return ix * ny_ * nz_ + iy * nz_ + iz;
    }

    param_t cutoff_;
    param_t cutoff_sq_;
    param_t inv_cutoff_;
    param_t origin_[3];
    index_t nx_, ny_, nz_;
    std::vector<std::vector<index_t>> cells_;
};

}
