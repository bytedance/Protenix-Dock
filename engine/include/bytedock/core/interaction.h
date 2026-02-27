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

#include "bytedock/core/cell_list.h"
#include "bytedock/core/data.h"
#include "bytedock/core/grid.h"
#include "bytedock/ext/counter.h"

namespace bytedock {

// bytedock.score_function.molecular.bond.MMBondEnergy
class harmonic_bond_interaction {
public:
    param_t put_gradients(const molecule_pose& mol_xyz,
                          const force_field_params& mol_ffdata,
                          molecule_pose& xyz_gradient);
};

// bytedock.score_function.molecular.angle.MMAngleEnergy
class harmonic_angle_interaction {
public:
    param_t put_gradients(const molecule_pose& mol_xyz,
                          const force_field_params& mol_ffdata,
                          molecule_pose& xyz_gradient);
};

// bytedock.score_function.molecular.torsion.MMTorsionEnergy
class fourier_dihedral_interaction {
public:
    param_t put_gradients(const molecule_pose& mol_xyz,
                          const force_field_params& mol_ffdata,
                          molecule_pose& xyz_gradient);
};

class self_nonbonded_interactions {
public:
    self_nonbonded_interactions() : coul14_scale_(0.833333333),
                                    vdw14_scale_(0.5) {}

    param_t put_gradients(const molecule_pose& mol_xyz,
                          const force_field_params& mol_ffdata,
                          molecule_pose& xyz_gradient);
#if BDOCK_HAS_SIMD  // For testing
    param_t put_gradients_nosimd(const molecule_pose& mol_xyz,
                                 const force_field_params& mol_ffdata,
                                 molecule_pose& xyz_gradient);
#endif

private:
    param_t coul14_scale_;
    param_t vdw14_scale_;
};

// bytedock.score_function.byte_mmenergy.MMMoleculeEnergy
class intra_molecular_interactions {
public:
    param_t put_gradients(const molecule_pose& mol_xyz,
                          const force_field_params& mol_ffdata,
                          molecule_pose& xyz_gradient);

private:
    harmonic_bond_interaction bond_;
    harmonic_angle_interaction angle_;
    fourier_dihedral_interaction dihedral_;
    self_nonbonded_interactions nonbonded_;
};

// bytedock.score_function.byte_score.MMInteractionEnergy
class inter_molecular_interactions {
public:
    explicit inter_molecular_interactions();
    DISABLE_COPY_AND_ASSIGN(inter_molecular_interactions);

    param_t put_gradients(const molecule_pose& receptor_xyz,
                          const force_field_params& receptor_ffdata,
                          const receptor_cache* nonbonded_cache,
                          const molecule_pose& ligand_xyz,
                          const force_field_params& ligand_ffdata,
                          molecule_pose& receptor_gradient,
                          molecule_pose& ligand_gradient);

    void build_cell_list(const molecule_pose& receptor_xyz, param_t cutoff = 10.0);

    // Pre-gather excluded atom data for cache-friendly access in hot loop
    void cache_excluded_data(const receptor_cache* cache,
                             const force_field_params& receptor_ffdata);

    void set_box_penalty(const box_config& bc, param_t slope) {
        bc_ = bc;
        penalty_slope_ = slope;
        has_box_penalty_ = true;
    }

private:
    cache_reporter* reporter_;
    cell_list bruteforce_cells_;

    // Pre-gathered excluded atom data (contiguous for cache locality)
    bool excl_data_cached_ = false;
    std::vector<index_t> excl_indices_;
    std::vector<param_t> excl_charges_;
    std::vector<lj_vdw> excl_vdw_params_;

    // Box penalty for brute-force path (no cache)
    bool has_box_penalty_ = false;
    box_config bc_;
    param_t penalty_slope_ = 0_r;
};

// bytedock.score_function.byte_mmenergy.MMTotalEnergy
class binding_system_interactions {
public:
    explicit binding_system_interactions(
        std::shared_ptr<torsional_receptor> receptor,
        std::shared_ptr<free_ligand> ligand,
        std::shared_ptr<receptor_cache> cache = nullptr,
        const box_config* box_penalty = nullptr,
        param_t penalty_slope = 0_r
    );

    param_t put_gradients(const std::vector<param_t>& receptor_torsions,
                          const molecule_pose& ligand_xyz,
                          std::vector<param_t>& torsion_gradient,
                          molecule_pose& ligand_gradient);

    void init_gradients(std::vector<param_t>& torsion_gradient,
                        molecule_pose& ligand_gradient);  // Resize & fill zeros
    void reset_gradients(std::vector<param_t>& torsion_gradient,
                         molecule_pose& ligand_gradient);  // Fill zeros

private:
    std::shared_ptr<torsional_receptor> receptor_;
    std::shared_ptr<free_ligand> ligand_;
    std::shared_ptr<receptor_cache> cache_;

    // It will update receptor_xyz_ & rotation_jacobian_
    void receptor_torsion2xyz(const std::vector<param_t>& torsions);
    // It will update receptor_gradient_
    void receptor_xyz2torsion(std::vector<param_t>& torsion_gradient);

    molecule_pose receptor_xyz_;  // Intermediate variables for gradient descent
    molecule_pose receptor_gradient_;  // Shape is [num_receptor_atoms, 3]
    molecule_pose rotation_jacobian_;  // Shape is [num_movable_atoms, 3]

    inter_molecular_interactions inter_;
    intra_molecular_interactions intra_;
};

}
