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

#include "bytedock/core/interaction.h"

#include "bytedock/core/constant.h"
#include "bytedock/core/pair.h"

#include <functional>

namespace bytedock {

param_t harmonic_bond_interaction::put_gradients(const molecule_pose& mol_xyz,
                                                 const force_field_params& mol_ffdata,
                                                 molecule_pose& xyz_gradient) {
    param_t energy = 0.;
    index_t ii, jj;
    param_t delta[3], r12, de_dr;
    for (auto& item : mol_ffdata.bonds) {
        ii = item.ids[0];
        jj = item.ids[1];
        minus_3d(mol_xyz[jj].xyz, mol_xyz[ii].xyz, delta);  // ii->jj
        r12 = get_norm_3d(delta);
        energy += 0.5_r * item.k * square(r12 - item.length);

        // Backward
        de_dr = item.k * (1 - safe_divide(item.length, r12));
        scale_inplace_3d(delta, de_dr);
        minus_inplace_3d(xyz_gradient[ii].xyz, delta);
        add_inplace_3d(xyz_gradient[jj].xyz, delta);
    }
    return energy;
}

param_t harmonic_angle_interaction::put_gradients(const molecule_pose& mol_xyz,
                                                  const force_field_params& mol_ffdata,
                                                  molecule_pose& xyz_gradient) {
    param_t energy = 0.;
    param_t theta, delta, de_dt, fi[3], fj[3], fk[3];
    index_t ii, jj, kk;
    for (auto& item : mol_ffdata.angles) {
        ii = item.ids[0];
        jj = item.ids[1];
        kk = item.ids[2];
        theta = get_angle_and_gradient_3d(mol_xyz[ii].xyz,
                                          mol_xyz[jj].xyz,
                                          mol_xyz[kk].xyz,
                                          fi, fk);
        delta = theta - item.radian;
        energy += 0.5_r * item.k * MATH_SQUARE_SCALAR(delta);
        de_dt = item.k * delta;
        scale_inplace_3d(fi, de_dt);
        scale_inplace_3d(fk, de_dt);
        add_inplace_3d(xyz_gradient[ii].xyz, fi);
        add_inplace_3d(xyz_gradient[kk].xyz, fk);
        add_3d(fi, fk, fj);
        minus_inplace_3d(xyz_gradient[jj].xyz, fj);
    }
    return energy;
}

inline param_t calculate_one_dihedral(const molecule_pose& mol_xyz,
                                      const fourier_dihedral& fd_params,
                                      molecule_pose& xyz_gradient) {
    const index_t ii = fd_params.ids[0], jj = fd_params.ids[1];
    const index_t kk = fd_params.ids[2], ll = fd_params.ids[3];
    param_t theta, grad_i[3], grad_j[3], grad_k[3], grad_l[3];
    theta = get_dihedral_angle_and_gradient_3d(
        mol_xyz[ii].xyz, mol_xyz[jj].xyz, mol_xyz[kk].xyz, mol_xyz[ll].xyz,
        grad_i, grad_j, grad_k, grad_l
    );
    param_t energy = 0., delta, de_dt, gg[3];
    for (size_t i = 0; i < FF_DIHEDRAL_NDIMS; i++) {
        if (fd_params.periodicity[i] == 0) break;
        delta = theta * fd_params.periodicity[i] - fd_params.phase[i];
        energy += fd_params.k[i] * (1 + std::cos(delta));
        de_dt = fd_params.k[i] * fd_params.periodicity[i] * std::sin(delta);
        scale_3d(grad_i, -de_dt, gg);
        add_inplace_3d(xyz_gradient[ii].xyz, gg);
        scale_3d(grad_j, -de_dt, gg);
        add_inplace_3d(xyz_gradient[jj].xyz, gg);
        scale_3d(grad_k, -de_dt, gg);
        add_inplace_3d(xyz_gradient[kk].xyz, gg);
        scale_3d(grad_l, -de_dt, gg);
        add_inplace_3d(xyz_gradient[ll].xyz, gg);
    }
    return energy;
}

param_t fourier_dihedral_interaction::put_gradients(
    const molecule_pose& mol_xyz,
    const force_field_params& mol_ffdata,
    molecule_pose& xyz_gradient
) {
    param_t energy = 0.;
    for (auto& item : mol_ffdata.impropers) {
        energy += calculate_one_dihedral(mol_xyz, item, xyz_gradient);
    }
    for (auto& item : mol_ffdata.propers) {
        energy += calculate_one_dihedral(mol_xyz, item,xyz_gradient);
    }
    return energy;
}

const lj_vdw* get_ligand_vdw_param(const index_t idx, const force_field_params& ffp) {
    return kVdwTypeDatabase + ffp.vdw_types[idx];
}

const lj_vdw* get_receptor_vdw_param(const index_t idx, const force_field_params& ffp) {
    return ffp.vdw_params.data() + idx;
}

param_t self_nonbonded_interactions::put_gradients(
    const molecule_pose& mol_xyz,
    const force_field_params& mol_ffdata,
    molecule_pose& xyz_gradient
) {
    param_t e1 = 0., e2 = 0.;
    index_t ii, jj;
    param_t r_ij[3], distance;
    std::function<
        const lj_vdw*(const index_t, const force_field_params&)
    > get_vdw_param = mol_ffdata.vdw_types.size() > 0
                    ? get_ligand_vdw_param : get_receptor_vdw_param;
    for (auto& item : mol_ffdata.pairs) {
        ii = item.ids[0];
        jj = item.ids[1];
        minus_3d(mol_xyz[jj].xyz, mol_xyz[ii].xyz, r_ij);  // i->j
        distance = get_norm_3d(r_ij);
        e1 += calculate_coul_pair(
            mol_ffdata.partial_charges[ii], mol_ffdata.partial_charges[jj],
            r_ij, distance, coul14_scale_, xyz_gradient[ii], xyz_gradient[jj]
        );
        e2 += calculate_vdw_pair(
            *get_vdw_param(ii, mol_ffdata), *get_vdw_param(jj, mol_ffdata),
            r_ij, distance, vdw14_scale_, xyz_gradient[ii], xyz_gradient[jj]
        );
    }
    for (auto& item : mol_ffdata.others) {
        ii = item.ids[0];
        jj = item.ids[1];
        minus_3d(mol_xyz[jj].xyz, mol_xyz[ii].xyz, r_ij);  // i->j
        distance = get_norm_3d(r_ij);
        e1 += calculate_coul_pair(
            mol_ffdata.partial_charges[ii], mol_ffdata.partial_charges[jj],
            r_ij, distance, 1., xyz_gradient[ii], xyz_gradient[jj]
        );
        e2 += calculate_vdw_pair(
            *get_vdw_param(ii, mol_ffdata), *get_vdw_param(jj, mol_ffdata),
            r_ij, distance, 1., xyz_gradient[ii], xyz_gradient[jj]
        );
    }
    return e1 + e2;
}

param_t intra_molecular_interactions::put_gradients(
    const molecule_pose& mol_xyz,
    const force_field_params& mol_ffdata,
    molecule_pose& xyz_gradient
) {
    /**
     * Activated items:
     * - "Bond_Energy" => MMBondEnergy(),
     * - "Angle_Energy"=> MMAngleEnergy(),
     * - "ProperTorsion_Energy" => MMTorsionEnergy(field="ProperTorsions"),
     * - "ImproperTorsion_Energy" => MMTorsionEnergy(field="ImproperTorsions"),
     * - "VDW_Energy" => MMVDWEnergy(index_key="FF_NonbondedAll_atomidx", scale=1.0),
     * - "VDW14_Energy" => MMVDWEnergy(index_key="FF_Nonbonded14_atomidx", scale=0.5),
     * - "Coulomb_Energy" => MMCoulombEnergy(
     *       index_key="FF_NonbondedAll_atomidx", scale=1.0
     *   ),
     * - "Coulomb14_Energy": MMCoulombEnergy(
     *       index_key="FF_Nonbonded14_atomidx", scale=0.833333333
     *   )
     */
    param_t energy = 0.;
    energy += bond_.put_gradients(mol_xyz, mol_ffdata, xyz_gradient);
    energy += angle_.put_gradients(mol_xyz, mol_ffdata, xyz_gradient);
    energy += dihedral_.put_gradients(mol_xyz, mol_ffdata, xyz_gradient);
    energy += nonbonded_.put_gradients(mol_xyz, mol_ffdata, xyz_gradient);
    return energy;
}

inter_molecular_interactions::inter_molecular_interactions() {
    reporter_ = cache_counter::singleton().get_or_create();
}

param_t inter_molecular_interactions::put_gradients(
    const molecule_pose& receptor_xyz,
    const force_field_params& receptor_ffdata,
    const receptor_cache* nonbonded_cache,
    const molecule_pose& ligand_xyz,
    const force_field_params& ligand_ffdata,
    molecule_pose& receptor_gradient,
    molecule_pose& ligand_gradient
) {
    /**
     * Activated items:
     * - "InteractionVDW_Energy" => MMInteractionVDWEnergy(method="nocutoff")
     * - "InteractionCoulomb_Energy" => MMInteractionCoulombEnergy()
     */
    param_t energy = 0., q_i, r_ij[3], distance;
    index_t type_idx;
    for (size_t i = 0; i < ligand_xyz.size(); i++) {
        q_i = ligand_ffdata.partial_charges[i];
        type_idx = ligand_ffdata.vdw_types[i];
        const lj_vdw& type_i = kVdwTypeDatabase[type_idx];
        if (nonbonded_cache) {  // Better speed but less precision
            if (nonbonded_cache->get(type_idx, q_i, ligand_xyz[i],
                                     energy, ligand_gradient[i])) {
                reporter_->hit();
            } else {
                reporter_->miss();
            }
            for (const auto& j : nonbonded_cache->get_excluded_atoms()) {
                minus_3d(receptor_xyz[j].xyz, ligand_xyz[i].xyz, r_ij);  // i->j
                distance = get_norm_3d(r_ij);
                energy += calculate_coul_pair(
                    q_i, receptor_ffdata.partial_charges[j],
                    r_ij, distance, 1., ligand_gradient[i], receptor_gradient[j]
                );
                energy += calculate_vdw_pair(
                    type_i, receptor_ffdata.vdw_params[j], r_ij, distance, 1.,
                    ligand_gradient[i], receptor_gradient[j]
                );
            }
        } else {  // Apply the brute-force calculation
            for (size_t j = 0; j < receptor_xyz.size(); j++) {
                minus_3d(receptor_xyz[j].xyz, ligand_xyz[i].xyz, r_ij);  // i->j
                distance = get_norm_3d(r_ij);
                energy += calculate_coul_pair(
                    q_i, receptor_ffdata.partial_charges[j],
                    r_ij, distance, 1., ligand_gradient[i], receptor_gradient[j]
                );
                energy += calculate_vdw_pair(
                    type_i, receptor_ffdata.vdw_params[j],
                    r_ij, distance, 1., ligand_gradient[i], receptor_gradient[j]
                );
            }
        }
    }
    return energy;
}

binding_system_interactions::binding_system_interactions(
    std::shared_ptr<torsional_receptor> receptor,
    std::shared_ptr<free_ligand> ligand,
    std::shared_ptr<receptor_cache> cache
) : receptor_(receptor), ligand_(ligand), cache_(cache) {
    receptor_xyz_ = receptor_->get_positions();  // Copy init coordinates
    receptor_gradient_.resize(receptor_xyz_.size());
    rotation_jacobian_.resize(receptor_->num_movable_atoms());
}

const atom_position kZeroPoint = {0., 0., 0.};

param_t binding_system_interactions::put_gradients(
    const std::vector<param_t>& receptor_torsions,
    const molecule_pose& ligand_xyz,
    std::vector<param_t>& torsion_gradient,
    molecule_pose& ligand_gradient
) {
    param_t energy = 0.;
    receptor_torsion2xyz(receptor_torsions);
    // TODO: Enable sparse update on position gradients of receptor atoms
    std::fill(receptor_gradient_.begin(), receptor_gradient_.end(), kZeroPoint);

    /**
     * Activated items:
     * - MMMoleculeEnergy()
     * - MMInteractionVDWEnergy()
     */
    const force_field_params& receptor_ffdata = receptor_->get_ffdata();
    const force_field_params& ligand_ffdata = ligand_->get_ffdata();
    energy += intra_.put_gradients(receptor_xyz_, receptor_ffdata, receptor_gradient_);
    energy += intra_.put_gradients(ligand_xyz, ligand_ffdata, ligand_gradient);
    energy += inter_.put_gradients(receptor_xyz_, receptor_ffdata, cache_.get(),
                                   ligand_xyz, ligand_ffdata,
                                   receptor_gradient_, ligand_gradient);

    // Apply chain rule
    receptor_xyz2torsion(torsion_gradient);
    return energy;
}

void binding_system_interactions::init_gradients(
    std::vector<param_t>& torsion_gradient, molecule_pose& ligand_gradient
) {
    torsion_gradient.clear();
    ligand_gradient.clear();
    torsion_gradient.resize(receptor_->num_torsions(), 0.);
    ligand_gradient.resize(ligand_->get_ffdata().vdw_types.size(), kZeroPoint);
}

void binding_system_interactions::reset_gradients(
    std::vector<param_t>& torsion_gradient, molecule_pose& ligand_gradient
) {
    std::fill(torsion_gradient.begin(), torsion_gradient.end(), 0.);
    std::fill(ligand_gradient.begin(), ligand_gradient.end(), kZeroPoint);
}

void binding_system_interactions::receptor_torsion2xyz(
    const std::vector<param_t>& torsions
) {
    matrix_3x3 rot, dv_dt;
    const molecule_pose& init_xyz = receptor_->get_positions();
    vector_3d tmp;
    size_t offset = 0;
    for(size_t i = 0; i < torsions.size(); i++) {
        const auto& leaf = receptor_->get_leaf(i);
        get_rotation_matrix_and_gradient(leaf.axis_vector, torsions[i], rot, dv_dt);

        // Apply rotation matrix to rotatable terminals
        for (const auto& j : leaf.movable_atoms) {
            minus_3d(init_xyz[j].xyz, leaf.axis_origin.xyz, tmp.xyz);
            multiply_3x3_3d(rot.data, tmp.xyz, receptor_xyz_[j].xyz);
            add_inplace_3d(receptor_xyz_[j].xyz, leaf.axis_origin.xyz);
            multiply_3x3_3d(dv_dt.data, tmp.xyz, rotation_jacobian_[offset].xyz);
            offset++;
        }
    }
}

void binding_system_interactions::receptor_xyz2torsion(
    std::vector<param_t>& torsion_gradient
) {
    size_t offset = 0;
    for(size_t i = 0; i < torsion_gradient.size(); i++) {
        const auto& leaf = receptor_->get_leaf(i);
        for (const auto& j : leaf.movable_atoms) {
            torsion_gradient[i] += dot_product_3d(rotation_jacobian_[offset].xyz,
                                                  receptor_gradient_[j].xyz);
            offset++;
        }
    }
}

}
