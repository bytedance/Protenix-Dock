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

#include <iostream>
#include <memory>

#include "bytedock/core/system.h"
#include "bytedock/lib/math.h"

namespace bytedock {

struct rotatable_leaf {
    vector_3d axis_vector;
    atom_position axis_origin;
    std::vector<index_t> movable_atoms;
};

class torsional_receptor {
public:
    // Input is generated by Python module `bdock_opt.app.prepare_receptor`
    void parse(std::shared_ptr<std::istream> fin);

    // The i-th value of `torsions` is applied to i-th terminal of `geo_`
    molecule_pose apply_parameters(const std::vector<param_t>& torsions) const;

    const force_field_params& get_ffdata() const { return ffp_; }
    const molecule_pose& get_positions() const { return coords_; }
    const rotatable_leaf& get_leaf(size_t idx) const { return leaves_[idx]; }
    const size_t num_torsions() const { return leaves_.size(); }
    const size_t num_movable_atoms() const { return nmovable_; }

    /**
     * Follow definitions of `rot_bond_order` in parser. Then, coordinates of receptor
     * atoms can be calculated by feeding rearranged torsions to method
     * `ReceptorGeometry.torsion2xyz`.
     */
    std::vector<param_t> restore_bond_orders(
        const std::vector<param_t>& torsions
    ) const {
        std::vector<param_t> aligned(geo_.num_rotatable_bonds, 0_r);
        for (size_t i = 0; i < torsions.size(); ++i) {
            aligned[geo_.bond_mappings[i]] = torsions[i];
        }
        return aligned;
    }

private:
    force_field_params ffp_;
    molecule_pose coords_;
    rotable_geometry geo_;

    // Post-processed fields from `geo_`
    std::vector<rotatable_leaf> leaves_;
    size_t nmovable_;
};

class free_ligand {
public:
    // Input is generated by Python module `bdock_opt.app.prepare_ligand`
    void parse(std::shared_ptr<std::istream> fin, const bool geometry = true);

    /**
     * The i-th value of `torsions` is applied to i-th rotatable bond during
     * a pre-order traverse of the fragment tree.
     */
    molecule_pose apply_parameters(const molecule_pose& ligand_xyz,
                                   const atom_position& transition,
                                   const matrix_3x3& orientation,
                                   const std::vector<param_t>& torsions) const;

    const std::string& get_smiles() const { return smiles_; }
    const force_field_params& get_ffdata() const { return ffp_; }
    const size_t num_torsions() const { return geo_.fragments.size() - 1; }

    size_t get_nposes() const { return poses_.size(); }
    const molecule_pose& get_pose(index_t idx) const { return poses_[idx]; }

private:
    std::string smiles_;
    force_field_params ffp_;
    fragment_tree geo_;
    std::vector<molecule_pose> poses_;
};

inline atom_position calculate_geometric_center(const molecule_pose& coords) {
    atom_position center = {0_r};
    for (const auto& one : coords) add_inplace_3d(center.xyz, one.xyz);
    scale_inplace_3d(center.xyz, 1_r / coords.size());
    return center;
}

inline param_t calculate_gyration_radius(const molecule_pose& coords,
                                         const atom_position& center) {
    param_t acc = 0_r;
    for (auto& item : coords) {
        acc += get_distance_square_3d(item.xyz, center.xyz);
    }
    return coords.size() > 0 ? std::sqrt(acc / coords.size()) : 0_r;
}

}
