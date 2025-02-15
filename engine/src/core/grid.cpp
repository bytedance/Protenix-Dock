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

#include "bytedock/core/grid.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <thread>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/unordered_set.hpp>
#include <boost/serialization/vector.hpp>

#include "bytedock/core/constant.h"
#include "bytedock/core/pair.h"
#include "bytedock/ext/counter.h"
#include "bytedock/ext/pfile.h"
#include "bytedock/lib/error.h"

namespace bytedock {

namespace fs = boost::filesystem;

box_config describe_boundary(double center_x, double center_y, double center_z,
                             double size_x, double size_y, double size_z) {
    box_config bc;
    param_t length = size_x * 0.5_r;
    bc.init_xyz[0] = center_x - length;
    bc.oppo_xyz[0] = center_x + length;
    length = size_y * 0.5_r;
    bc.init_xyz[1] = center_y - length;
    bc.oppo_xyz[1] = center_y + length;
    length = size_z * 0.5_r;
    bc.init_xyz[2] = center_z - length;
    bc.oppo_xyz[2] = center_z + length;
    return bc;
}

inline void fill_derived_fields(box_config& bc) {
    bc.spacing_inv = 1_r / bc.spacing;
    for (size_t i = 0; i < 3; ++i) {
        bc.npoints_minus_one[i] = bc.npoints[i] - 1;
        bc.oppo_xyz[i] = bc.init_xyz[i] + bc.npoints_minus_one[i] * bc.spacing;
    }
}

box_config decide_grid_dims(double center_x, double center_y, double center_z,
                            double size_x, double size_y, double size_z,
                            double spacing) {
    index_t nvoxels[3];
    nvoxels[0] = static_cast<index_t>(std::ceil(size_x / spacing));
    nvoxels[1] = static_cast<index_t>(std::ceil(size_y / spacing));
    nvoxels[2] = static_cast<index_t>(std::ceil(size_z / spacing));
    box_config bc;
    bc.init_xyz[0] = center_x;
    bc.init_xyz[1] = center_y;
    bc.init_xyz[2] = center_z;
    param_t length;
    for (size_t i = 0; i < 3; ++i) {
        if (nvoxels[i] % 2 == 1) nvoxels[i] += 1;
        length = nvoxels[i] * spacing;
        bc.npoints[i] = nvoxels[i] + 1;
        bc.init_xyz[i] -= length * 0.5_r;
    }
    bc.spacing = spacing;
    fill_derived_fields(bc);
    return bc;
}

template<class Archive>
void serialize(Archive & ar, val_and_grad& vg, const unsigned int version) {
    ar & vg.first;
    ar & vg.second;
}

template<class Archive>
void serialize(Archive & ar, atom_position& ap, const unsigned int version) {
    ar & ap.xyz[0];
    ar & ap.xyz[1];
    ar & ap.xyz[2];
}

receptor_cache::receptor_cache() : slope_(0_r) {
    vdw_grids_.resize(kNumVdwTypes);
}

receptor_cache::receptor_cache(const box_config& bc) : bc_(bc), slope_(0_r) {
    vdw_grids_.resize(kNumVdwTypes);
    index_t total = bc_.npoints[0] * bc_.npoints[1] * bc_.npoints[2];
    if (total > 0) {
        // Zero initialization is ensured in C++11
        coul_grid_.resize(total);
        for (index_t i = 0; i < kNumVdwTypes; ++i) vdw_grids_[i].resize(total);
    }
}

void receptor_cache::populate(const torsional_receptor& receptor, size_t nthreads) {
    step_timer st(CONSTRUCT_RECEPTOR_CACHE_STEP);

    // Pick up movable atoms
    excluded_.clear();
    for (size_t i = 0; i < receptor.num_torsions(); ++i) {
        const auto& movable = receptor.get_leaf(i).movable_atoms;
        excluded_.insert(movable.begin(), movable.end());
    }

    // Precalculate with frozen atoms
    const size_t ni = bc_.npoints[0];
    if (nthreads > ni) nthreads = ni;
    progress_bar pbar(ni * bc_.npoints[1] * bc_.npoints[2]);
    if (nthreads < 2) {
        partial_worker builder(*this);
        builder.populate(receptor, pbar, 0, ni);
    } else {
        std::vector<std::thread> builders;
        index_t length = ni / nthreads;
        if (ni % nthreads > 0) length += 1;
        index_t start, end;
        for (size_t i = 0; i < nthreads; ++i) {
            start = i * length;
            end = start + length;
            if (end > ni) end = ni;
            builders.emplace_back([&, start, end]() {
                partial_worker builder(*this);
                builder.populate(receptor, pbar, start, end);
            });
        }
        for (size_t i = 0; i < builders.size(); ++i) builders[i].join();
    }
}

void receptor_cache::partial_worker::populate(
    const torsional_receptor& receptor, progress_bar& pbar,
    const index_t start_i, const index_t end_i
) {
    const auto& bc = cache_.bc_;
    const auto& excluded = cache_.excluded_;
    const auto& coords = receptor.get_positions();
    const auto& qlist = receptor.get_ffdata().partial_charges;
    const lj_vdw* params = receptor.get_ffdata().vdw_params.data();
    param_t point_xyz[3], vec[3], distance;
    atom_position tmp_grad = {0_r};
    val_and_grad coul_acc;
    std::array<val_and_grad, kNumVdwTypes> vdw_acc;
    index_t offset;
    for (size_t i = start_i; i < end_i; ++i) {
        point_xyz[0] = bc.init_xyz[0] + bc.spacing * i;
        for (size_t j = 0; j < bc.npoints[1]; ++j) {
            point_xyz[1] = bc.init_xyz[1] + bc.spacing * j;
            for (size_t k = 0; k < bc.npoints[2]; ++k) {
                point_xyz[2] = bc.init_xyz[2] + bc.spacing * k;
                coul_acc = {0_r};
                vdw_acc.fill({0_r});
                for (size_t atom_idx = 0; atom_idx < coords.size(); ++atom_idx) {
                    if (excluded.find(atom_idx) != excluded.end()) continue;
                    minus_3d(coords[atom_idx].xyz, point_xyz, vec);
                    distance = get_norm_3d(vec);
                    coul_acc.first += calculate_coul_pair(1_r, qlist[atom_idx],
                                                          vec, distance, 1_r,
                                                          coul_acc.second, tmp_grad);
                    for (size_t type_idx = 0; type_idx < kNumVdwTypes; ++type_idx) {
                        vdw_acc[type_idx].first += calculate_vdw_pair(
                            kVdwTypeDatabase[type_idx], params[atom_idx],
                            vec, distance, 1_r, vdw_acc[type_idx].second, tmp_grad
                        );
                    }
                }
                offset = cache_.get_offset(i, j, k);
                cache_.coul_grid_[offset] = coul_acc;
                for (size_t type_idx = 0; type_idx < kNumVdwTypes; ++type_idx) {
                    cache_.vdw_grids_[type_idx][offset] = vdw_acc[type_idx];
                }
                ++pbar;
            }
        }
    }
}

bool receptor_cache::get(const index_t atom_type,
                         const param_t partial_charge,
                         const atom_position& atom_xyz,
                         param_t& out_energy,
                         atom_position& xyz_gradient) const {
    param_t flb_xyz[3];
    minus_3d(atom_xyz.xyz, bc_.init_xyz, flb_xyz);
    scale_inplace_3d(flb_xyz, bc_.spacing_inv);

    // Ensure ligand atom is within receptor box
    index_t flb_ijk[3], brt_ijk[3];
    param_t grad_tmp[3], penalty = 0_r;
    for (size_t dim = 0; dim < 3; ++dim) {
        if (flb_xyz[dim] < 0_r) {
            flb_ijk[dim] = 0;
            grad_tmp[dim] = -1_r;
            penalty -= flb_xyz[dim];
            flb_xyz[dim] = 0_r;
        } else if (flb_xyz[dim] < bc_.npoints_minus_one[dim]) {
            flb_ijk[dim] = static_cast<index_t>(flb_xyz[dim]);
            grad_tmp[dim] = 0_r;
        } else {
            flb_ijk[dim] = bc_.npoints_minus_one[dim] - 1;
            grad_tmp[dim] = 1_r;
            penalty += flb_xyz[dim] - bc_.npoints_minus_one[dim];
            flb_xyz[dim] = bc_.npoints_minus_one[dim];
        }
        brt_ijk[dim] = flb_ijk[dim] + 1;
    }

    // Calculate position index
    position_index pi;
    pi.p000 = get_offset(flb_ijk[0], flb_ijk[1], flb_ijk[2]);
    pi.p100 = get_offset(brt_ijk[0], flb_ijk[1], flb_ijk[2]);
    pi.p010 = get_offset(flb_ijk[0], brt_ijk[1], flb_ijk[2]);
    pi.p110 = get_offset(brt_ijk[0], brt_ijk[1], flb_ijk[2]);
    pi.p001 = get_offset(flb_ijk[0], flb_ijk[1], brt_ijk[2]);
    pi.p101 = get_offset(brt_ijk[0], flb_ijk[1], brt_ijk[2]);
    pi.p011 = get_offset(flb_ijk[0], brt_ijk[1], brt_ijk[2]);
    pi.p111 = get_offset(brt_ijk[0], brt_ijk[1], brt_ijk[2]);
    pi.x = flb_xyz[0] - flb_ijk[0];
    pi.y = flb_xyz[1] - flb_ijk[1];
    pi.z = flb_xyz[2] - flb_ijk[2];
    pi.mx = 1_r - pi.x;
    pi.my = 1_r - pi.y;
    pi.mz = 1_r - pi.z;

    // Accumulate energy & gradient from penalty
    bool hit = true;
    if (penalty > 0_r) {
        out_energy += penalty * bc_.spacing * slope_;
        scale_inplace_3d(grad_tmp, slope_);  // bc_.spacing * bc_.spacing_inv = 1
        add_inplace_3d(xyz_gradient.xyz, grad_tmp);
        hit = false;
    }

    // Search energy & gradient
    out_energy += get(coul_grid_, pi, grad_tmp) * partial_charge;
    if (hit) {
        scale_inplace_3d(grad_tmp, partial_charge);
        add_inplace_3d(xyz_gradient.xyz, grad_tmp);
    }
    out_energy += get(vdw_grids_[atom_type], pi, grad_tmp);
    if (hit) add_inplace_3d(xyz_gradient.xyz, grad_tmp);
    return hit;
}

inline void scale_inplace(val_and_grad& vg, const param_t factor) {
    vg.first *= factor;
    scale_inplace_3d(vg.second.xyz, factor);
}
inline void add_inplace(val_and_grad& vg, const val_and_grad& other) {
    vg.first += other.first;
    add_inplace_3d(vg.second.xyz, other.second.xyz);
}

param_t receptor_cache::get(const std::vector<val_and_grad>& points,
                            const position_index& pi,
                            param_t* xyz_gradient) const {
    val_and_grad f000 = points[pi.p000];
    val_and_grad f100 = points[pi.p100];
    val_and_grad f010 = points[pi.p010];
    val_and_grad f110 = points[pi.p110];
    val_and_grad f001 = points[pi.p001];
    val_and_grad f101 = points[pi.p101];
    val_and_grad f011 = points[pi.p011];
    val_and_grad f111 = points[pi.p111];
    scale_inplace(f000, pi.mx * pi.my * pi.mz);
    scale_inplace(f100, pi.x  * pi.my * pi.mz);
    scale_inplace(f010, pi.mx * pi.y  * pi.mz);
    scale_inplace(f110, pi.x  * pi.y  * pi.mz);
    scale_inplace(f001, pi.mx * pi.my * pi.z);
    scale_inplace(f101, pi.x  * pi.my * pi.z);
    scale_inplace(f011, pi.mx * pi.y  * pi.z);
    scale_inplace(f111, pi.x  * pi.y  * pi.z);
    add_inplace(f000, f100);
    add_inplace(f000, f010);
    add_inplace(f000, f110);
    add_inplace(f000, f001);
    add_inplace(f000, f101);
    add_inplace(f000, f011);
    add_inplace(f000, f111);
    xyz_gradient[0] = f000.second.xyz[0];
    xyz_gradient[1] = f000.second.xyz[1];
    xyz_gradient[2] = f000.second.xyz[2];
    return f000.first;
}

const std::string kMapFileExtension = ".bdg";
const std::string kMetaFileName = "meta.dat";
const index_t kFormatVersion = 20241030;

void receptor_cache::save(const fs::path& map_dir, size_t nthreads) {
    step_timer st(SAVE_RECEPTOR_CACHE_STEP);

    // Save grids in parallel
    const size_t ni = kNumVdwTypes + 1;
    progress_bar pbar(ni);
    if (nthreads < 2) {
        partial_worker pw(*this);
        pw.save(map_dir, pbar, 0, ni);
    } else {
        std::vector<std::thread> savers;
        index_t length = ni / nthreads;
        if (ni % nthreads > 0) length += 1;
        index_t start, end;
        for (size_t i = 0; i < nthreads; ++i) {
            start = i * length;
            end = start + length;
            if (end > ni) end = ni;
            savers.emplace_back([&, start, end]() {
                partial_worker pw(*this);
                pw.save(map_dir, pbar, start, end);
            });
        }
        for (size_t i = 0; i < savers.size(); ++i) savers[i].join();
    }

    // Save metadata at last and it can be considered as a flag of completion
    auto meta_file = map_dir / kMetaFileName;
    auto ofs = open_for_write(meta_file.string());
    boost::archive::text_oarchive ar(*ofs);
    ar << kFormatVersion;
    ar << bc_.init_xyz;
    ar << bc_.npoints;
    ar << bc_.spacing;
    ar << excluded_;
}

void receptor_cache::restore(const fs::path& map_dir, size_t nthreads) {
    step_timer st(RESTORE_RECEPTOR_CACHE_STEP);

    // Check version compatibility first
    auto meta_file = map_dir / kMetaFileName;
    auto path = meta_file.string();
    if (!is_file(path)) {
        throw file_system_error("The metadata file of maps is missing!");
    }
    auto ifs = open_for_read(meta_file.string());
    boost::archive::text_iarchive ar(*ifs);
    index_t version;
    ar >> version;
    if (version != kFormatVersion) {
        std::ostringstream oss;
        oss << "Format version of map files are not compatible! "
            << "Expected [" << kFormatVersion << "] vs. Got [" << version << "]";
        throw failed_conf_error(oss.str());
    }
    ar >> bc_.init_xyz;
    ar >> bc_.npoints;
    ar >> bc_.spacing;
    ar >> excluded_;
    fill_derived_fields(bc_);

    // Restore grids in parallel
    const size_t ni = kNumVdwTypes + 1;
    progress_bar pbar(ni);
    if (nthreads < 2) {
        partial_worker pw(*this);
        pw.restore(map_dir, pbar, 0, ni);
    } else {
        std::vector<std::thread> loaders;
        index_t length = ni / nthreads;
        if (ni % nthreads > 0) length += 1;
        index_t start, end;
        for (size_t i = 0; i < nthreads; ++i) {
            start = i * length;
            end = start + length;
            if (end > ni) end = ni;
            loaders.emplace_back([&, start, end]() {
                partial_worker pw(*this);
                pw.restore(map_dir, pbar, start, end);
            });
        }
        for (size_t i = 0; i < loaders.size(); ++i) loaders[i].join();
    }
}

void receptor_cache::partial_worker::save(const fs::path& map_dir, progress_bar& pbar,
                                          const index_t start_i, const index_t end_i) {
    for (size_t i = start_i; i < end_i; ++i) {
        auto data_file = map_dir / (kVdwTypeTags[i] + kMapFileExtension);
        auto ofs = open_for_write(data_file.string());
        boost::archive::binary_oarchive oa(*ofs);
        if (i == kNumVdwTypes) {
            oa << cache_.coul_grid_;
        } else if (i < cache_.vdw_grids_.size()) {
            oa << cache_.vdw_grids_[i];
        }
        ++pbar;
    }
}

void receptor_cache::partial_worker::restore(
    const fs::path& map_dir, progress_bar& pbar,
    const index_t start_i, const index_t end_i
) {
    for (size_t i = start_i; i < end_i; ++i) {
        auto data_file = map_dir / (kVdwTypeTags[i] + kMapFileExtension);
        auto path = data_file.string();
        if (!is_file(path)) {
            std::ostringstream oss;
            oss << "Map file for atom type [" << kVdwTypeTags[i] << "] is missing!";
            throw file_system_error(oss.str());
        }
        auto ifs = open_for_read(path);
        boost::archive::binary_iarchive ia(*ifs);
        if (i == kNumVdwTypes) {
            ia >> cache_.coul_grid_;
        } else if (i < cache_.vdw_grids_.size()) {
            ia >> cache_.vdw_grids_[i];
        }
        ++pbar;
    }
}

}
