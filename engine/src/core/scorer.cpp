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

#include "bytedock/core/scorer.h"

#include "bytedock/lib/utility.h"

namespace bytedock {

std::unordered_map<std::string, param_t> abstract_scorer::summary(
    const molecule_pose& receptor_xyz,
    const force_field_params& receptor_ffdata,
    const molecule_pose& ligand_xyz,
    const force_field_params& ligand_ffdata,
    const instant_cache* pose_memo
) const {
    std::unordered_map<std::string, param_t> merged;
    for (auto& scorer: children_) {
        auto subset = scorer->summary(receptor_xyz, receptor_ffdata,
                                      ligand_xyz, ligand_ffdata, pose_memo);
        merged.insert(subset.begin(), subset.end());
    }
    return merged;
}

std::unordered_map<std::string, param_t> leaf_scorer::summary(
    const molecule_pose& receptor_xyz,
    const force_field_params& receptor_ffdata,
    const molecule_pose& ligand_xyz,
    const force_field_params& ligand_ffdata,
    const instant_cache* pose_memo
) const {
    std::unordered_map<std::string, param_t> scores;
    scores[name_] = report(receptor_xyz, receptor_ffdata,
                           ligand_xyz, ligand_ffdata, pose_memo);
    return scores;
}

root_scorer::root_scorer(
    const std::string& name, const ScoringFunctionParams& sfp
) : sf_name_(name) {
    if CHECK_PARAMETER_NON_ZERO(sfp.interaction_coulomb_energy_coef) {
        add_coefficient("InteractionCoulomb_oriEnergy",
                        sfp.interaction_coulomb_energy_coef);
        add_child(new mm_interaction_coulomb_energy("InteractionCoulomb_oriEnergy"));
    }
    if CHECK_PARAMETER_NON_ZERO(sfp.interaction_vdw_energy_coef) {
        add_coefficient("InteractionVDW_Energy",
                        sfp.interaction_vdw_energy_coef);
        add_child(new mm_interaction_vdw_energy(
            "InteractionVDW_Energy", sfp.interaction_vdw_energy
        ));
    }
    if CHECK_PARAMETER_NON_ZERO(sfp.ionic_energy_coef) {
        add_coefficient("Ionic_Energy", sfp.ionic_energy_coef);
        add_child(new ionic_interaction_total_energy(
            "Ionic_Energy", sfp.ionic_energy
        ));
    }
    if CHECK_PARAMETER_NON_ZERO(sfp.pi_cationic_energy_coef) {
        add_coefficient("Pi_Cationic_Energy", sfp.pi_cationic_energy_coef);
        add_child(new cation_pi_interaction_total_energy(
            "Pi_Cationic_Energy", sfp.pi_cationic_energy
        ));
    }
    if CHECK_PARAMETER_NON_ZERO(sfp.f2f_pi_stacking_energy_coef) {
        add_coefficient("F2FPiStacking_Energy", sfp.f2f_pi_stacking_energy_coef);
        add_child(new pi_stacking_total_energy(
            "F2FPiStacking_Energy", 0_r, sfp.f2f_pi_stacking_energy
        ));
    }
    if CHECK_PARAMETER_NON_ZERO(sfp.e2f_pi_stacking_energy_coef) {
        add_coefficient("E2FPiStacking_Energy", sfp.e2f_pi_stacking_energy_coef);
        add_child(new pi_stacking_total_energy(
            "E2FPiStacking_Energy", 90_r, sfp.e2f_pi_stacking_energy
        ));
    }
    if CHECK_PARAMETER_NON_ZERO(sfp.hbond_energy_coef) {
        add_coefficient("Hbond_Energy", sfp.hbond_energy_coef);
        add_child(new hbond_energy("Hbond_Energy", sfp.hbond_energy));
    }
    if CHECK_PARAMETER_NON_ZERO(sfp.hydrophobic_energy_coef) {
        add_coefficient("Hydrophobic_Energy", sfp.hydrophobic_energy_coef);
        add_child(new hydrophobic_energy(
            "Hydrophobic_Energy", sfp.hydrophobic_energy
        ));
    }
    if CHECK_PARAMETER_NON_ZERO(sfp.torsion_strain_penalty_coef) {
        add_coefficient("TorsionStrain_Penalty", sfp.torsion_strain_penalty_coef);
        add_child(new torsion_strain_penalty_scorer("TorsionStrain_Penalty"));
    }
    if CHECK_PARAMETER_NON_ZERO(sfp.rotatable_energy_coef) {
        add_coefficient("TorsionStrain_Penalty", sfp.rotatable_energy_coef);
        add_child(new rotatable_energy("rotatable_energy"));
    }
    if CHECK_PARAMETER_NON_ZERO(sfp.gbsa_energy_delta_coef) {
        add_coefficient("Gbsa_Com_Energy-Gbsa_Lig_Energy-Gbsa_Pro_Energy",
                        sfp.gbsa_energy_delta_coef);
        add_child(new gbsa_total_energy);
    }
}

param_t root_scorer::combine(
    const std::unordered_map<std::string, param_t>& scores
) const {
    param_t merged = 0.;
    for (auto& iter: coefficients_) {
        auto value = scores.find(iter.first);
        if (value != scores.end()) merged += iter.second * value->second;
    }
    return merged;
}

inline void to_ring_centroid(const molecule_pose& coords,
                             const pi_ring& pr,
                             atom_position& out) {
    param_t* xyz = out.xyz;
    xyz[0] = xyz[1] = xyz[2] = 0_r;
    for (size_t i = 0; i < pr.size; i++) add_inplace_3d(xyz, coords[pr.ids[i]].xyz);
    xyz[0] /= pr.size;
    xyz[1] /= pr.size;
    xyz[2] /= pr.size;
}

inline std::vector<atom_position> get_ring_centroids(
    const molecule_pose& pose, const std::vector<pi_ring>& rings
) {
    std::vector<atom_position> centers(rings.size());
    size_t i = 0;
    while (i < centers.size()) {
        to_ring_centroid(pose, rings[i], centers[i]);
        ++i;
    }
    return centers;
}

inline void to_ring_normal(const molecule_pose& coords, const pi_ring& pr,
                           vector_3d& out) {
    const atom_position& pos1 = coords[pr.ids[0]];
    const atom_position& pos2 = coords[pr.ids[1]];
    const atom_position& pos3 = coords[pr.ids[2]];
    param_t edge1[3], edge2[3];
    minus_3d(pos2.xyz, pos1.xyz, edge1);
    minus_3d(pos2.xyz, pos3.xyz, edge2);
    param_t* xyz = out.xyz;
    cross_product_3d(edge1, edge2, xyz);
    normalize_3d(xyz);
}

inline std::vector<vector_3d> get_ring_normals(
    const molecule_pose& pose, const std::vector<pi_ring>& rings
) {
    std::vector<vector_3d> normals(rings.size());
    size_t i = 0;
    while (i < normals.size()) {
        to_ring_normal(pose, rings[i], normals[i]);
        i++;
    }
    return normals;
}

instant_cache root_scorer::precalculate(
    const molecule_pose& receptor_xyz, const force_field_params& receptor_ffdata,
    const molecule_pose& ligand_xyz, const force_field_params& ligand_ffdata
) const {
    instant_cache memo;
    memo.distance_map.resize(ligand_xyz.size() * receptor_xyz.size());
    size_t offset = 0;
    for (size_t i = 0; i < ligand_xyz.size(); i++) {
        const atom_position& xyz1 = ligand_xyz[i];
        for (size_t j = 0; j < receptor_xyz.size(); j++) {
            memo.distance_map[offset] = get_distance_3d(xyz1.xyz, receptor_xyz[j].xyz);
            offset++;
        }
    }
    memo.receptor_ring5_centroids = get_ring_centroids(receptor_xyz,
                                                       receptor_ffdata.pi5_rings);
    memo.receptor_ring6_centroids = get_ring_centroids(receptor_xyz,
                                                       receptor_ffdata.pi6_rings);
    memo.ligand_ring5_centroids = get_ring_centroids(ligand_xyz,
                                                     ligand_ffdata.pi5_rings);
    memo.ligand_ring6_centroids = get_ring_centroids(ligand_xyz,
                                                     ligand_ffdata.pi6_rings);
    memo.receptor_ring5_normals = get_ring_normals(receptor_xyz,
                                                   receptor_ffdata.pi5_rings);
    memo.receptor_ring6_normals = get_ring_normals(receptor_xyz,
                                                   receptor_ffdata.pi6_rings);
    memo.ligand_ring5_normals = get_ring_normals(ligand_xyz, ligand_ffdata.pi5_rings);
    memo.ligand_ring6_normals = get_ring_normals(ligand_xyz, ligand_ffdata.pi6_rings);
    return memo;
}

mm_interaction_vdw_energy::mm_interaction_vdw_energy(
    const std::string& name, const InteractionVdwEnergyParams& params
) : leaf_scorer(name) {
    scale_ = params.scale;
    energy_cutoff_ = params.cutoff;
    const lj_vdw* origin;
    for (size_t i = 0; i < kNumVdwTypes; ++i) {
        origin = kVdwTypeDatabase + i;
        scaled_vdw_types_[i].sigma = origin->sigma * params.types[i].sigma_scale;
        scaled_vdw_types_[i].epsilon = origin->epsilon * params.types[i].epsilon_scale;
    }
}

#define FROM_2D_INDEX(idx, idy, dimy) ((idx)*(dimy)+(idy))

param_t mm_interaction_vdw_energy::report(
    const molecule_pose& receptor_xyz,
    const force_field_params& receptor_ffdata,
    const molecule_pose& ligand_xyz,
    const force_field_params& ligand_ffdata,
    const instant_cache* pose_memo
) const {
    param_t energy = 0.;
    param_t distance, sigma1, epsilon1, fused, u6;
    for (size_t i = 0; i < ligand_xyz.size(); i++) {
        const lj_vdw& type1 = scaled_vdw_types_[ligand_ffdata.vdw_types[i]];
        sigma1 = type1.sigma;
        epsilon1 = type1.epsilon;
        for (size_t j = 0; j < receptor_xyz.size(); j++) {
            const lj_vdw& type2 = receptor_ffdata.vdw_params[j];
            if (pose_memo == nullptr) {
                distance = get_distance_3d(ligand_xyz[i].xyz, receptor_xyz[j].xyz);
            } else {
                distance = pose_memo->distance_map[
                    FROM_2D_INDEX(i, j, receptor_xyz.size())
                ];
            }
            fused = (sigma1 + type2.sigma) * 0.5_r;
            u6 = std::pow(safe_divide(fused, distance), 6_r);
            fused = std::sqrt(epsilon1 * type2.epsilon) * 4_r;
            distance = (MATH_SQUARE_SCALAR(u6) - u6) * fused; // Used as `tmp`
            energy += clamp(distance, -energy_cutoff_, energy_cutoff_);
        }
    }
    return energy * scale_;
}

param_t mm_interaction_coulomb_energy::report(
    const molecule_pose& receptor_xyz, const force_field_params& receptor_ffdata,
    const molecule_pose& ligand_xyz, const force_field_params& ligand_ffdata,
    const instant_cache* pose_memo
) const {
    param_t energy = 0., distance, q1;
    for (size_t i = 0; i < ligand_xyz.size(); i++) {
        q1 = ligand_ffdata.partial_charges[i];
        for (size_t j = 0; j < receptor_xyz.size(); j++) {
            if (pose_memo == nullptr) {
                distance = get_distance_3d(ligand_xyz[i].xyz, receptor_xyz[j].xyz);
            } else {
                distance = pose_memo->distance_map[
                    FROM_2D_INDEX(i, j, receptor_xyz.size())
                ];
            }
            energy += safe_divide(
                q1 * receptor_ffdata.partial_charges[j], distance
            ) * kCoulombFactor;
        }
    }
    return energy;
}

param_t hydrophobic_energy::report(const molecule_pose& receptor_xyz,
                                           const force_field_params& receptor_ffdata,
                                           const molecule_pose& ligand_xyz,
                                           const force_field_params& ligand_ffdata,
                                           const instant_cache* pose_memo
) const {
    param_t energy = 0_r;
    const auto& ligand_hydro_atoms = ligand_ffdata.hydrophobic_atoms;
    const auto& receptor_hydro_atoms = receptor_ffdata.hydrophobic_atoms;
    param_t sigma1, sigma2, threshold_r1, threshold_r2, distance;
    index_t ligand_atomidx, receptor_atomidx;
    for (size_t i = 0; i < ligand_hydro_atoms.size(); i++) {
        ligand_atomidx = ligand_hydro_atoms[i];
        sigma1 = kVdwTypeDatabase[ligand_ffdata.vdw_types[ligand_atomidx]].sigma;
        for (size_t j = 0; j < receptor_hydro_atoms.size(); j++) {
            receptor_atomidx = receptor_hydro_atoms[j];
            sigma2 = receptor_ffdata.vdw_params[receptor_atomidx].sigma;
            threshold_r1 = (sigma1 + sigma2) * scale_factor_ + dist_min_;
            threshold_r2 = threshold_r1 + dist_max_;
            if (pose_memo == nullptr) {
                distance = get_distance_3d(ligand_xyz[ligand_atomidx].xyz,
                                           receptor_xyz[receptor_atomidx].xyz);
            } else {
                distance = pose_memo->distance_map[FROM_2D_INDEX(
                    ligand_atomidx, receptor_atomidx, receptor_xyz.size()
                )];
            }
            energy -= get_reversed_linear_percent(distance,
                                                  threshold_r1,
                                                  threshold_r2);
        }
    }
    return energy * scale_;
}

ionic_interaction_total_energy::ionic_interaction_total_energy(
    const std::string& name, const IonicEnergyParams& params
) : leaf_scorer(name) {
    lc_coef_ = params.lc_coef;
    pc_coef_ = params.pc_coef;
    lc_inter_.reset(new ionic_interaction_batch(params.lc_inter, true));
    pc_inter_.reset(new ionic_interaction_batch(params.pc_inter, false));
}

param_t ionic_interaction_total_energy::report(
    const molecule_pose& receptor_xyz,
    const force_field_params& receptor_ffdata,
    const molecule_pose& ligand_xyz,
    const force_field_params& ligand_ffdata,
    const instant_cache* pose_memo
) const {
    param_t energy = 0.;
    energy += lc_coef_ * lc_inter_->calculate(
        ligand_xyz, ligand_ffdata.cation_atoms,
        receptor_xyz, receptor_ffdata.anion_atoms, pose_memo
    );
    energy += pc_coef_ * pc_inter_->calculate(
        receptor_xyz, receptor_ffdata.cation_atoms,
        ligand_xyz, ligand_ffdata.anion_atoms, pose_memo
    );
    return energy;
}

param_t ionic_interaction_total_energy::ionic_interaction_batch::calculate(
    const molecule_pose& cation_mol_xyz,
    const std::vector<index_t>& cation_atoms,
    const molecule_pose& anion_mol_xyz,
    const std::vector<index_t>& anion_atoms,
    const instant_cache* pose_memo
) const {
    param_t energy = 0.;
    param_t distance;
    for (size_t i = 0; i < cation_atoms.size(); i++) {
        for (size_t j = 0; j < anion_atoms.size(); j++) {
            if (pose_memo == nullptr) {
                distance = get_distance_3d(cation_mol_xyz[cation_atoms[i]].xyz,
                                           anion_mol_xyz[anion_atoms[j]].xyz);
            } else if (from_ligand_) {
                distance = pose_memo->distance_map[FROM_2D_INDEX(
                    cation_atoms[i], anion_atoms[j], anion_mol_xyz.size()
                )];
            } else {  // Cations belong to receptor
                distance = pose_memo->distance_map[FROM_2D_INDEX(
                    anion_atoms[j], cation_atoms[i], cation_mol_xyz.size()
                )];
            }
            energy -= get_reversed_linear_percent(distance, dist_min_, dist_max_);
        }
    }
    return energy;
}

cation_pi_interaction_total_energy::cation_pi_interaction_total_energy(
    const std::string& name, const PiCationicEnergyParams& params
) : leaf_scorer(name) {
    lc5_coef_ = params.lc5_coef;
    lc6_coef_ = params.lc6_coef;
    pc5_coef_ = params.pc5_coef;
    pc6_coef_ = params.pc6_coef;
    lc5_inter_.reset(new cation_pi_interaction_batch(params.lc5_inter));
    lc6_inter_.reset(new cation_pi_interaction_batch(params.lc6_inter));
    pc5_inter_.reset(new cation_pi_interaction_batch(params.pc5_inter));
    pc6_inter_.reset(new cation_pi_interaction_batch(params.pc6_inter));
}

param_t cation_pi_interaction_total_energy::report(
    const molecule_pose& receptor_xyz,
    const force_field_params& receptor_ffdata,
    const molecule_pose& ligand_xyz,
    const force_field_params& ligand_ffdata,
    const instant_cache* pose_memo
) const {
    param_t energy = 0.;
    if (pose_memo == nullptr) {
        energy += lc5_coef_ * lc5_inter_->calculate(
            ligand_xyz, ligand_ffdata.cation_atoms,
            receptor_xyz, receptor_ffdata.pi5_rings
        );
        energy += lc6_coef_ * lc6_inter_->calculate(
            ligand_xyz, ligand_ffdata.cation_atoms,
            receptor_xyz, receptor_ffdata.pi6_rings
        );
        energy += pc5_coef_ * pc5_inter_->calculate(
            receptor_xyz, receptor_ffdata.cation_atoms,
            ligand_xyz, ligand_ffdata.pi5_rings
        );
        energy += pc6_coef_ * pc6_inter_->calculate(
            receptor_xyz, receptor_ffdata.cation_atoms,
            ligand_xyz, ligand_ffdata.pi6_rings
        );
    } else {
        energy += lc5_coef_ * lc5_inter_->calculate(
            ligand_xyz, ligand_ffdata.cation_atoms,
            pose_memo->receptor_ring5_centroids
        );
        energy += lc6_coef_ * lc6_inter_->calculate(
            ligand_xyz, ligand_ffdata.cation_atoms,
            pose_memo->receptor_ring6_centroids
        );
        energy += pc5_coef_ * pc5_inter_->calculate(
            receptor_xyz, receptor_ffdata.cation_atoms,
            pose_memo->ligand_ring5_centroids
        );
        energy += pc6_coef_ * pc6_inter_->calculate(
            receptor_xyz, receptor_ffdata.cation_atoms,
            pose_memo->ligand_ring6_centroids
        );
    }
    return energy;
}

param_t cation_pi_interaction_total_energy::cation_pi_interaction_batch::calculate(
    const molecule_pose& cation_mol_xyz,
    const std::vector<index_t>& cation_atoms,
    const molecule_pose& ring_mol_xyz,
    const std::vector<pi_ring>& rings
) const {
    param_t energy = 0_r, distance;
    atom_position center;
    for (size_t i = 0; i < rings.size(); i++) {
        to_ring_centroid(ring_mol_xyz, rings[i], center);
        for (size_t j = 0; j < cation_atoms.size(); j++) {
            distance = get_distance_3d(center.xyz, cation_mol_xyz[cation_atoms[j]].xyz);
            energy -= get_reversed_linear_percent(distance, dist_min_, dist_max_);
        }
    }
    return energy;
}

param_t cation_pi_interaction_total_energy::cation_pi_interaction_batch::calculate(
    const molecule_pose& cation_mol_xyz,
    const std::vector<index_t>& cation_atoms,
    const std::vector<atom_position>& ring_centroids
) const {
    param_t energy = 0_r, distance;
    const param_t *center;
    for (size_t i = 0; i < ring_centroids.size(); i++) {
        center = ring_centroids[i].xyz;
        for (size_t j = 0; j < cation_atoms.size(); j++) {
            distance = get_distance_3d(center, cation_mol_xyz[cation_atoms[j]].xyz);
            energy -= get_reversed_linear_percent(distance, dist_min_, dist_max_);
        }
    }
    return energy;
}

hbond_energy::hbond_energy(
    const std::string& name, const HbondEnergyParams& params
) : leaf_scorer(name) {
    scale_ = params.scale;
    hldcc_coef_ = params.hldcc_coef;
    hldcn_coef_ = params.hldcn_coef;
    hldnc_coef_ = params.hldnc_coef;
    hldnn_coef_ = params.hldnn_coef;
    hpdcc_coef_ = params.hpdcc_coef;
    hpdcn_coef_ = params.hpdcn_coef;
    hpdnc_coef_ = params.hpdnc_coef;
    hpdnn_coef_ = params.hpdnn_coef;
    hldcc_inter_.reset(new hbond_interactions_ga_tuning_param(params.hldcc_inter));
    hldcn_inter_.reset(new hbond_interactions_ga_tuning_param(params.hldcn_inter));
    hldnc_inter_.reset(new hbond_interactions_ga_tuning_param(params.hldnc_inter));
    hldnn_inter_.reset(new hbond_interactions_ga_tuning_param(params.hldnn_inter));
    hpdcc_inter_.reset(new hbond_interactions_ga_tuning_param(params.hpdcc_inter));
    hpdcn_inter_.reset(new hbond_interactions_ga_tuning_param(params.hpdcn_inter));
    hpdnc_inter_.reset(new hbond_interactions_ga_tuning_param(params.hpdnc_inter));
    hpdnn_inter_.reset(new hbond_interactions_ga_tuning_param(params.hpdnn_inter));
}

param_t hbond_energy::report(const molecule_pose& receptor_xyz,
                                     const force_field_params& receptor_ffdata,
                                     const molecule_pose& ligand_xyz,
                                     const force_field_params& ligand_ffdata,
                                     const instant_cache* pose_memo) const {
    param_t energy = 0.;
    energy += hldcc_coef_ * hldcc_inter_->calculate(
        ligand_xyz, ligand_ffdata.hbonddon_charged,
        receptor_xyz, receptor_ffdata.hbondacc_charged,
        true, pose_memo
    );
    energy += hldcn_coef_ * hldcn_inter_->calculate(
        ligand_xyz, ligand_ffdata.hbonddon_charged,
        receptor_xyz, receptor_ffdata.hbondacc_neutral,
        true, pose_memo
    );
    energy += hldnc_coef_ * hldnc_inter_->calculate(
        ligand_xyz, ligand_ffdata.hbonddon_neutral,
        receptor_xyz, receptor_ffdata.hbondacc_charged,
        true, pose_memo
    );
    energy += hldnn_coef_ * hldnn_inter_->calculate(
        ligand_xyz, ligand_ffdata.hbonddon_neutral,
        receptor_xyz, receptor_ffdata.hbondacc_neutral,
        true, pose_memo
    );
    energy += hpdcc_coef_ * hpdcc_inter_->calculate(
        receptor_xyz, receptor_ffdata.hbonddon_charged,
        ligand_xyz, ligand_ffdata.hbondacc_charged,
        false, pose_memo
    );
    energy += hpdcn_coef_ * hpdcn_inter_->calculate(
        receptor_xyz, receptor_ffdata.hbonddon_charged,
        ligand_xyz, ligand_ffdata.hbondacc_neutral,
        false, pose_memo
    );
    energy += hpdnc_coef_ * hpdnc_inter_->calculate(
        receptor_xyz, receptor_ffdata.hbonddon_neutral,
        ligand_xyz, ligand_ffdata.hbondacc_charged,
        false, pose_memo
    );
    energy += hpdnn_coef_ * hpdnn_inter_->calculate(
        receptor_xyz, receptor_ffdata.hbonddon_neutral,
        ligand_xyz, ligand_ffdata.hbondacc_neutral,
        false, pose_memo
    );
    return energy * scale_;
}

param_t hbond_energy::hbond_interactions_ga_tuning_param::calculate(
    const molecule_pose& donor_coords,
    const std::vector<atom_pair>& donor_pairs,
    const molecule_pose& acceptor_coords,
    const std::vector<index_t>& acceptor_atoms,
    const bool donor_from_ligand,
    const instant_cache* pose_memo
) const {
    param_t energy = 0.;
    param_t distance, angle;
    vector_3d h2don_vector, h2acc_vector;
    size_t h_atomidx, acc_atomidx;
    for (size_t i = 0; i < donor_pairs.size(); i++) {
        const atom_position& donor_pos = donor_coords[donor_pairs[i].ids[0]];
        h_atomidx = donor_pairs[i].ids[1];
        const atom_position& h_pos = donor_coords[h_atomidx];
        for (size_t j = 0; j < acceptor_atoms.size(); j++) {
            acc_atomidx = acceptor_atoms[j];
            const atom_position& acc_pos = acceptor_coords[acc_atomidx];
            if (pose_memo == nullptr) {
                distance = get_distance_3d(h_pos.xyz, acc_pos.xyz);
            } else if (donor_from_ligand) {
                distance = pose_memo->distance_map[
                    FROM_2D_INDEX(h_atomidx, acc_atomidx, acceptor_coords.size())
                ];
            } else {  // Donor belongs to receptor
                distance = pose_memo->distance_map[
                    FROM_2D_INDEX(acc_atomidx, h_atomidx, donor_coords.size())
                ];
            }
            distance = get_reversed_linear_percent(distance, dist_center_,
                                                   dist_min_, dist_max_);
            if (distance < kParamEpsilon) continue;
            minus_3d(donor_pos.xyz, h_pos.xyz, h2don_vector.xyz);
            minus_3d(acc_pos.xyz, h_pos.xyz, h2acc_vector.xyz);
            angle = MATH_TO_DEGREE(get_vector_angle(h2don_vector, h2acc_vector));
            angle = get_reversed_linear_percent(angle, angle_center_,
                                                angle_min_, angle_max_);
            energy -= distance * angle;
        }
    }
    return energy;
}

const index_t kTspPeriods[] = {1, 2, 3, 4};
const param_t kTspPhases[] = {0_r, kMathPi, 0_r, kMathPi};

param_t torsion_strain_penalty_scorer::report(
    const molecule_pose& receptor_xyz, const force_field_params& receptor_ffdata,
    const molecule_pose& ligand_xyz, const force_field_params& ligand_ffdata,
    const instant_cache* pose_memo
) const {
    const index_t* atoms;
    param_t score = 0_r, theta, tmp;
    size_t j;
    for (size_t i = 0; i < ligand_ffdata.bad_torsions.size(); ++i) {
        const auto& bt = ligand_ffdata.bad_torsions[i];
        tmp = 0_r;
        for (j = 0; j < bt.matched.size(); ++j) {
            atoms = bt.matched[j].ids;
            if (atoms[0] == atoms[1]) break;  // Same as `match_mask`
            theta = get_dihedral_angle_3d(
                ligand_xyz[atoms[0]].xyz, ligand_xyz[atoms[1]].xyz,
                ligand_xyz[atoms[2]].xyz, ligand_xyz[atoms[3]].xyz
            );
            tmp += (
                1 + std::cos(theta * kTspPeriods[0] - kTspPhases[0])
            ) * bt.coefs[0] + (
                1 + std::cos(theta * kTspPeriods[1] - kTspPhases[1])
            ) * bt.coefs[1] + (
                1 + std::cos(theta * kTspPeriods[2] - kTspPhases[2])
            ) * bt.coefs[2] + (
                1 + std::cos(theta * kTspPeriods[3] - kTspPhases[3])
            ) * bt.coefs[3];
        }
        if (j > 0) score += tmp / j;
    }
    return score;
}

inline param_t get_projected_centroid_distance(const vector_3d& normal1,
                                               const atom_position& center1,
                                               const atom_position& center2) {
    param_t tmp[3];
    minus_3d(center2.xyz, center1.xyz, tmp);
    param_t dp = dot_product_3d(normal1.xyz, tmp);
    return std::sqrt(square_sum_3d(tmp) - MATH_SQUARE_SCALAR(dp));
}

inline param_t get_sigma_distance(const vector_3d& receptor_normal,
                                  const atom_position& receptor_center,
                                  const vector_3d& ligand_normal,
                                  const atom_position& ligand_center) {
    param_t sigma1, sigma2;
    sigma1 = get_projected_centroid_distance(
        receptor_normal, receptor_center, ligand_center
    );
    sigma2 = get_projected_centroid_distance(
        ligand_normal, ligand_center, receptor_center
    );
    return MATH_PAIR_MIN(sigma1, sigma2);
}

pi_stacking_total_energy::pi_stacking_total_energy(
    const std::string& name, const param_t center, const PiStackingEnergyParams& params
) : leaf_scorer(name) {
    l5r5_coef_ = params.l5r5_coef;
    l5r6_coef_ = params.l5r6_coef;
    l6r5_coef_ = params.l6r5_coef;
    l6r6_coef_ = params.l6r6_coef;
    l5r5_inter_.reset(new pi_stacking_batch(center, params.l5r5_inter));
    l5r6_inter_.reset(new pi_stacking_batch(center, params.l5r6_inter));
    l6r5_inter_.reset(new pi_stacking_batch(center, params.l6r5_inter));
    l6r6_inter_.reset(new pi_stacking_batch(center, params.l6r6_inter));
}

param_t pi_stacking_total_energy::report(
    const molecule_pose& receptor_xyz,
    const force_field_params& receptor_ffdata,
    const molecule_pose& ligand_xyz,
    const force_field_params& ligand_ffdata,
    const instant_cache* pose_memo
) const {
    param_t energy = 0_r;
    if (pose_memo == nullptr) {
        auto receptor_r5c = get_ring_centroids(receptor_xyz, receptor_ffdata.pi5_rings);
        auto receptor_r6c = get_ring_centroids(receptor_xyz, receptor_ffdata.pi6_rings);
        auto ligand_r5c = get_ring_centroids(ligand_xyz, ligand_ffdata.pi5_rings);
        auto ligand_r6c = get_ring_centroids(ligand_xyz, ligand_ffdata.pi6_rings);
        auto receptor_r5n = get_ring_normals(receptor_xyz, receptor_ffdata.pi5_rings);
        auto receptor_r6n = get_ring_normals(receptor_xyz, receptor_ffdata.pi6_rings);
        auto ligand_r5n = get_ring_normals(ligand_xyz, ligand_ffdata.pi5_rings);
        auto ligand_r6n = get_ring_normals(ligand_xyz, ligand_ffdata.pi6_rings);
        energy += l5r5_coef_ * l5r5_inter_->calculate(receptor_r5c, receptor_r5n,
                                                      ligand_r5c, ligand_r5n);
        energy += l5r6_coef_ * l5r6_inter_->calculate(receptor_r6c, receptor_r6n,
                                                      ligand_r5c, ligand_r5n);
        energy += l6r5_coef_ * l6r5_inter_->calculate(receptor_r5c, receptor_r5n,
                                                      ligand_r6c, ligand_r6n);
        energy += l6r6_coef_ * l6r6_inter_->calculate(receptor_r6c, receptor_r6n,
                                                      ligand_r6c, ligand_r6n);
    } else {
        energy += l5r5_coef_ * l5r5_inter_->calculate(
            pose_memo->receptor_ring5_centroids, pose_memo->receptor_ring5_normals,
            pose_memo->ligand_ring5_centroids, pose_memo->ligand_ring5_normals
        );
        energy += l5r6_coef_ * l5r6_inter_->calculate(
            pose_memo->receptor_ring6_centroids, pose_memo->receptor_ring6_normals,
            pose_memo->ligand_ring5_centroids, pose_memo->ligand_ring5_normals
        );
        energy += l6r5_coef_ * l6r5_inter_->calculate(
            pose_memo->receptor_ring5_centroids, pose_memo->receptor_ring5_normals,
            pose_memo->ligand_ring6_centroids, pose_memo->ligand_ring6_normals
        );
        energy += l6r6_coef_ * l6r6_inter_->calculate(
            pose_memo->receptor_ring6_centroids, pose_memo->receptor_ring6_normals,
            pose_memo->ligand_ring6_centroids, pose_memo->ligand_ring6_normals
        );
    }
    return energy;
}

param_t pi_stacking_total_energy::pi_stacking_batch::calculate(
    const std::vector<atom_position>& receptor_centers,
    const std::vector<vector_3d>& receptor_normals,
    const std::vector<atom_position>& ligand_centers,
    const std::vector<vector_3d>& ligand_normals
) const {
    param_t energy = 0_r, distance, sigma, theta;
    for (size_t i = 0; i < receptor_centers.size(); i++) {
        auto& receptor_center = receptor_centers[i];
        auto& receptor_normal = receptor_normals[i];
        for (size_t j = 0; j < ligand_centers.size(); j++) {
            distance = get_distance_3d(receptor_center.xyz, ligand_centers[j].xyz);
            if (distance > dist_max_) continue;
            distance = get_reversed_linear_percent(distance, dist_min_, dist_max_);
            sigma = get_sigma_distance(receptor_normal, receptor_center,
                                       ligand_normals[j], ligand_centers[j]);
            if (sigma > sigma_max_) continue;
            sigma = get_reversed_linear_percent(sigma, sigma_min_, sigma_max_);
            theta = MATH_TO_DEGREE(get_normal_angle(receptor_normal,
                                                    ligand_normals[j]));
            theta = std::min(theta, 180_r - theta);
            theta = get_reversed_linear_percent(theta, center_,
                                                angle_min_, angle_max_);
            energy -= distance * sigma * theta;
        }
    }
    return energy;
}

param_t rotatable_energy::report(const molecule_pose& receptor_xyz,
                                    const force_field_params& receptor_ffdata,
                                    const molecule_pose& ligand_xyz,
                                    const force_field_params& ligand_ffdata,
                                    const instant_cache* pose_memo) const {
    auto nbonds = static_cast<param_t>(ligand_ffdata.rotatable_bonds.size());
    return -1_r / (2.620643_r + nbonds);
}

inline void fill_atom_properties(const force_field_params& mol_ffdata,
                                 std::vector<param_t>& radii,
                                 std::vector<param_t>& screens) {
    const size_t natoms = mol_ffdata.atomic_numbers.size();
    radii.resize(natoms, -1.);
    index_t elem_a, elem_b;
    for (auto& pair : mol_ffdata.all_bonds) {
        elem_a = mol_ffdata.atomic_numbers[pair.ids[0]];
        elem_b = mol_ffdata.atomic_numbers[pair.ids[1]];
        if (elem_a == 1 && elem_b == 7) {
            radii[pair.ids[0]] = 1.3;
        } else if (elem_a == 7 && elem_b == 1) {
            radii[pair.ids[1]] = 1.3;
        }
    }
    screens.resize(natoms);
    for (size_t i = 0; i < natoms; ++i) {
        auto& elem = mol_ffdata.atomic_numbers[i];
        if (radii[i] < 0_r) {
            radii[i] = elem <= kNumRadii ? kRadiusDatabase[elem] : kRadiusDatabase[0];
        }
        screens[i] = elem <= kNumScreens ? kScreenDatabase[elem] : kScreenDatabase[0];
    }
}

std::unordered_map<std::string, param_t> gbsa_total_energy::summary(
    const molecule_pose& receptor_xyz,
    const force_field_params& receptor_ffdata,
    const molecule_pose& ligand_xyz,
    const force_field_params& ligand_ffdata,
    const instant_cache* pose_memo
) const {
    // Calculate partial terms
    std::vector<param_t> receptor_radii, receptor_screens;
    fill_atom_properties(receptor_ffdata, receptor_radii, receptor_screens);
    param_t gbsa_pro_energy = calculate(receptor_xyz, receptor_ffdata.partial_charges,
                                        receptor_radii, receptor_screens);
    std::vector<param_t> radii, screens;
    fill_atom_properties(ligand_ffdata, radii, screens);
    param_t gbsa_lig_energy = calculate(ligand_xyz, ligand_ffdata.partial_charges,
                                        radii, screens);

    // Calculate the complex term
    molecule_pose coords(ligand_xyz);  // Copy
    coords.insert(coords.end(), receptor_xyz.begin(), receptor_xyz.end());
    std::vector<param_t> charges(ligand_ffdata.partial_charges);  // Copy
    charges.insert(charges.end(), receptor_ffdata.partial_charges.begin(),
                                  receptor_ffdata.partial_charges.end());
    radii.insert(radii.end(), receptor_radii.begin(), receptor_radii.end());
    screens.insert(screens.end(), receptor_screens.begin(), receptor_screens.end());
    param_t gbsa_com_energy = calculate(coords, charges, radii, screens);

    // Calculate the energy delta
    std::unordered_map<std::string, param_t> scores;
    scores["Gbsa_Pro_Energy"] = gbsa_pro_energy;
    scores["Gbsa_Lig_Energy"] = gbsa_lig_energy;
    scores["Gbsa_Com_Energy"] = gbsa_com_energy;
    scores["Gbsa_Com_Energy-Gbsa_Lig_Energy-Gbsa_Pro_Energy"] = gbsa_com_energy
                                                              - gbsa_lig_energy
                                                              - gbsa_pro_energy;
    return scores;
}

param_t gbsa_total_energy::calculate(
    const molecule_pose& coords,
    const std::vector<param_t>& charges,
    const std::vector<param_t>& radii,
    const std::vector<param_t>& screens
) const {
    const size_t natoms = radii.size();
    std::vector<param_t> offset_radii(natoms);
    for (size_t i = 0; i < natoms; ++i) offset_radii[i] = radii[i] - offset_;

    // Fill values of `sum_I`
    param_t screened_radius, screened_radius_squared;
    param_t r12, r12_inv, U, D, L, U_inv, L_inv;
    std::vector<param_t> sum_I(natoms, 0.);
    for (size_t i = 0; i < natoms; ++i) {  // The second dim in shape [natoms, natoms]
        screened_radius = offset_radii[i] * screens[i];
        screened_radius_squared = MATH_SQUARE_SCALAR(screened_radius);
        for (size_t j = 0; j < natoms; ++j) {  // The first dim in shape
            if (i == j) continue;
            r12 = get_distance_3d(coords[i].xyz, coords[j].xyz);
            U = r12 + screened_radius;
            if (U < offset_radii[j]) continue;
            r12_inv = 1_r / r12;
            D = std::abs(r12 - screened_radius);
            L = MATH_PAIR_MAX(offset_radii[j], D);
            U_inv = 1_r / U;
            L_inv = 1_r / L;
            D = 0.25_r * (r12 - screened_radius_squared * r12_inv)
              * (MATH_SQUARE_SCALAR(U_inv) - MATH_SQUARE_SCALAR(L_inv));  // => `I`
            D += (L_inv - U_inv + 0.5_r * safe_ln(L * U_inv) * r12_inv);
            D *= 0.5_r;
            sum_I[j] += D;
        }
    }

    // Fill values of `born_radii`
    std::vector<param_t> born_radii(natoms);
    for (size_t i = 0; i < natoms; ++i) {
        U = sum_I[i] * offset_radii[i];  // Used as `psi`
        D = MATH_SQUARE_SCALAR(U);  // Used as `psi**2`
        L = std::tanh(alpha_ * U - beta_ * D + gamma_ * U * D);
        born_radii[i] = 1_r / (1_r / offset_radii[i] - L / radii[i]);
    }

    // Calculate egb & esa
    param_t egb = 0.;
    for (size_t i = 0; i < natoms; ++i) {
        for (size_t j = 0; j < natoms; ++j) {
            if (i == j) {
                U = born_radii[i];  // U => f
            } else {
                // r12 => r12_squared, r12_inv => born_radius_product, U => f
                r12 = get_distance_square_3d(coords[i].xyz, coords[j].xyz);
                r12_inv = born_radii[i] * born_radii[j];
                U = safe_sqrt(r12 + r12_inv * std::exp(-0.25_r * r12 / r12_inv));
            }
            L = 0.5_r * kCoulombFactor * (
                1_r / solute_dielectric_ - std::exp(-kappa_ * U) / solvent_dielectric_
            );  // => `prefactor`
            egb -= charges[i] * charges[j] / U * L;
        }
    }
    param_t esa = 0.;
    for (size_t i = 0; i < natoms; ++i) {
        esa += square(radii[i] + probe_radius_)
             * std::pow(radii[i] / born_radii[i], 6_r);
    }
    return egb + esa * surface_energy_;
}

}
