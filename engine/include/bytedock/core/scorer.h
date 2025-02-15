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

#include <memory>
#include <unordered_map>

#include "bytedock/core/conf.h"
#include "bytedock/core/constant.h"
#include "bytedock/core/data.h"

namespace bytedock {

// If coordinates of receptor or ligand change, cache should be reconstructed
struct instant_cache {
    // Shape is [ligand_natoms, receptor_natoms]
    std::vector<param_t> distance_map;

    // Shape is [pi_nrings]
    std::vector<atom_position> receptor_ring5_centroids;
    std::vector<atom_position> receptor_ring6_centroids;
    std::vector<atom_position> ligand_ring5_centroids;
    std::vector<atom_position> ligand_ring6_centroids;
    std::vector<vector_3d> receptor_ring5_normals;
    std::vector<vector_3d> receptor_ring6_normals;
    std::vector<vector_3d> ligand_ring5_normals;
    std::vector<vector_3d> ligand_ring6_normals;
};

class abstract_scorer {
public:
    // Only leaf scorers override this method
    virtual std::unordered_map<std::string, param_t> summary(
        const molecule_pose& receptor_xyz,
        const force_field_params& receptor_ffdata,
        const molecule_pose& ligand_xyz,
        const force_field_params& ligand_ffdata,
        const instant_cache* pose_memo = nullptr
    ) const;

protected:
    // Ownership of the pointer is transferred as well
    void add_child(abstract_scorer* scorer) {
        children_.push_back(std::unique_ptr<abstract_scorer>(scorer));
    }

private:
    std::vector<std::unique_ptr<abstract_scorer> > children_;
};

class leaf_scorer : public abstract_scorer {
public:
    leaf_scorer(const std::string& name) : name_(name) {}

    std::unordered_map<std::string, param_t> summary(
        const molecule_pose& receptor_xyz,
        const force_field_params& receptor_ffdata,
        const molecule_pose& ligand_xyz,
        const force_field_params& ligand_ffdata,
        const instant_cache* pose_memo = nullptr
    ) const override;

    virtual param_t report(const molecule_pose& receptor_xyz,
                           const force_field_params& receptor_ffdata,
                           const molecule_pose& ligand_xyz,
                           const force_field_params& ligand_ffdata,
                           const instant_cache* pose_memo = nullptr) const = 0;

private:
    const std::string name_;
};

/**
 * This is not desgined for thread safety! Please ensure each call of method `summary`
 * is from and only from one thread!
 */
class root_scorer : public abstract_scorer {
public:
    root_scorer(const std::string& name, const ScoringFunctionParams& sfp);

    param_t combine(const std::unordered_map<std::string, param_t>& scores) const;
    param_t combine(const molecule_pose& receptor_xyz,
                    const force_field_params& receptor_ffdata,
                    const molecule_pose& ligand_xyz,
                    const force_field_params& ligand_ffdata,
                    const instant_cache& pose_memo) const {
        auto scores = summary(receptor_xyz, receptor_ffdata,
                              ligand_xyz, ligand_ffdata, &pose_memo);
        return combine(scores);
    }

    // Alway enable precalculation for better speed
    param_t combine(const molecule_pose& receptor_xyz,
                    const force_field_params& receptor_ffdata,
                    const molecule_pose& ligand_xyz,
                    const force_field_params& ligand_ffdata) const {
        auto pose_memo = precalculate(receptor_xyz, receptor_ffdata,
                                      ligand_xyz, ligand_ffdata);
        auto scores = summary(receptor_xyz, receptor_ffdata,
                              ligand_xyz, ligand_ffdata, &pose_memo);
        return combine(scores);
    }

    /**
     * Each top-level scorer knows which filed can be shared across sub scorers. Thus,
     * it is in charge of filling instant cache.
     */
    instant_cache precalculate(const molecule_pose& receptor_xyz,
                               const force_field_params& receptor_ffdata,
                               const molecule_pose& ligand_xyz,
                               const force_field_params& ligand_ffdata) const;

protected:
    void add_coefficient(const std::string& name, param_t value) {
        coefficients_[name] = value;
    }

private:
    const std::string sf_name_;  // For debug
    std::unordered_map<std::string, param_t> coefficients_;
};

class mm_interaction_vdw_energy : public leaf_scorer {
public:
    mm_interaction_vdw_energy(const std::string& name);
    mm_interaction_vdw_energy(const std::string& name,
                                      const InteractionVdwEnergyParams& params);

    param_t report(
        const molecule_pose& receptor_xyz,
        const force_field_params& receptor_ffdata,
        const molecule_pose& ligand_xyz,
        const force_field_params& ligand_ffdata,
        const instant_cache* pose_memo = nullptr
    ) const override;

private:
    param_t scale_;
    param_t energy_cutoff_;

    // It is enabled by option `adopt_tuned_ligvdw_params=True` in pscore settings
    lj_vdw scaled_vdw_types_[kNumVdwTypes];
};

class mm_interaction_coulomb_energy : public leaf_scorer {
public:
    mm_interaction_coulomb_energy(const std::string& name) : leaf_scorer(name) {}

    param_t report(const molecule_pose& receptor_xyz,
                   const force_field_params& receptor_ffdata,
                   const molecule_pose& ligand_xyz,
                   const force_field_params& ligand_ffdata,
                   const instant_cache* pose_memo = nullptr) const override;
};

class hydrophobic_energy : public leaf_scorer {
public:
    hydrophobic_energy(
        const std::string& name, const HydrophobicEnergyParams& params
    ) : leaf_scorer(name), scale_(params.scale), dist_min_(params.dist_min),
        dist_max_(params.dist_max), scale_factor_(params.scale_factor) {}

    param_t report(const molecule_pose& receptor_xyz,
                   const force_field_params& receptor_ffdata,
                   const molecule_pose& ligand_xyz,
                   const force_field_params& ligand_ffdata,
                   const instant_cache* pose_memo = nullptr) const override;

private:
    param_t scale_;
    param_t dist_min_;
    param_t dist_max_;
    param_t scale_factor_;
};

class ionic_interaction_total_energy : public leaf_scorer {
public:
    ionic_interaction_total_energy(const std::string& name,
                                      const IonicEnergyParams& params);

    param_t report(const molecule_pose& receptor_xyz,
                   const force_field_params& receptor_ffdata,
                   const molecule_pose& ligand_xyz,
                   const force_field_params& ligand_ffdata,
                   const instant_cache* pose_memo = nullptr) const override;

private:
    class ionic_interaction_batch {
    public:
        ionic_interaction_batch(
            const IonicInteractionParams& params, bool from_ligand
        ) : dist_min_(params.dist_min), dist_max_(params.dist_max),
            from_ligand_(from_ligand) {}

        param_t calculate(const molecule_pose& cation_mol_xyz,
                          const std::vector<index_t>& cation_atoms,
                          const molecule_pose& anion_mol_xyz,
                          const std::vector<index_t>& anion_atoms,
                          const instant_cache* pose_memo) const;

    private:
        param_t dist_min_;
        param_t dist_max_;
        bool from_ligand_;
    };

    param_t lc_coef_;
    param_t pc_coef_;
    std::unique_ptr<ionic_interaction_batch> lc_inter_;
    std::unique_ptr<ionic_interaction_batch> pc_inter_;
};

class cation_pi_interaction_total_energy : public leaf_scorer {
public:
    cation_pi_interaction_total_energy(const std::string& name,
                                       const PiCationicEnergyParams& params);

    param_t report(const molecule_pose& receptor_xyz,
                   const force_field_params& receptor_ffdata,
                   const molecule_pose& ligand_xyz,
                   const force_field_params& ligand_ffdata,
                   const instant_cache* pose_memo = nullptr) const override;

private:
    class cation_pi_interaction_batch {
    public:
        cation_pi_interaction_batch(
            const CationPiInteractionParams& params
        ) : dist_min_(params.dist_min), dist_max_(params.dist_max) {}

        param_t calculate(const molecule_pose& cation_mol_xyz,
                          const std::vector<index_t>& cation_atoms,
                          const molecule_pose& ring_mol_xyz,
                          const std::vector<pi_ring>& rings) const;

        param_t calculate(const molecule_pose& cation_mol_xyz,
                          const std::vector<index_t>& cation_atoms,
                          const std::vector<atom_position>& ring_centroids) const;

    private:
        param_t dist_min_;
        param_t dist_max_;
    };

    param_t lc5_coef_;
    param_t lc6_coef_;
    param_t pc5_coef_;
    param_t pc6_coef_;
    std::unique_ptr<cation_pi_interaction_batch> lc5_inter_;
    std::unique_ptr<cation_pi_interaction_batch> lc6_inter_;
    std::unique_ptr<cation_pi_interaction_batch> pc5_inter_;
    std::unique_ptr<cation_pi_interaction_batch> pc6_inter_;
};

class hbond_energy : public leaf_scorer {
public:
    hbond_energy(const std::string& name, const HbondEnergyParams& params);

    param_t report(const molecule_pose& receptor_xyz,
                   const force_field_params& receptor_ffdata,
                   const molecule_pose& ligand_xyz,
                   const force_field_params& ligand_ffdata,
                   const instant_cache* pose_memo = nullptr) const override;

private:
    class hbond_interactions_ga_tuning_param {
    public:
        hbond_interactions_ga_tuning_param(
            const HbondInteractionParams& params
        ) : dist_min_(params.dist_min), dist_max_(params.dist_max),
            angle_min_(params.angle_min), angle_max_(params.angle_max),
            dist_center_(params.dist_center), angle_center_(params.angle_center) {}

        param_t calculate(const molecule_pose& donor_coords,
                          const std::vector<atom_pair>& donor_pairs,
                          const molecule_pose& acceptor_coords,
                          const std::vector<index_t>& acceptor_atoms,
                          const bool donor_from_ligand,
                          const instant_cache* pose_memo) const;

    private:
        param_t dist_min_;
        param_t dist_max_;
        param_t angle_min_;
        param_t angle_max_;
        param_t dist_center_;
        param_t angle_center_;
    };

private:
    param_t scale_;
    param_t hldcc_coef_;
    param_t hldcn_coef_;
    param_t hldnc_coef_;
    param_t hldnn_coef_;
    param_t hpdcc_coef_;
    param_t hpdcn_coef_;
    param_t hpdnc_coef_;
    param_t hpdnn_coef_;
    std::unique_ptr<hbond_interactions_ga_tuning_param> hldcc_inter_;
    std::unique_ptr<hbond_interactions_ga_tuning_param> hldcn_inter_;
    std::unique_ptr<hbond_interactions_ga_tuning_param> hldnc_inter_;
    std::unique_ptr<hbond_interactions_ga_tuning_param> hldnn_inter_;
    std::unique_ptr<hbond_interactions_ga_tuning_param> hpdcc_inter_;
    std::unique_ptr<hbond_interactions_ga_tuning_param> hpdcn_inter_;
    std::unique_ptr<hbond_interactions_ga_tuning_param> hpdnc_inter_;
    std::unique_ptr<hbond_interactions_ga_tuning_param> hpdnn_inter_;
};

class torsion_strain_penalty_scorer : public leaf_scorer {
public:
    torsion_strain_penalty_scorer(const std::string& name) : leaf_scorer(name) {}

    param_t report(const molecule_pose& receptor_xyz,
                   const force_field_params& receptor_ffdata,
                   const molecule_pose& ligand_xyz,
                   const force_field_params& ligand_ffdata,
                   const instant_cache* pose_memo = nullptr) const override;
};

class pi_stacking_total_energy : public leaf_scorer {
public:
    pi_stacking_total_energy(const std::string& name,
                             const param_t center,
                             const PiStackingEnergyParams& params);

    param_t report(
        const molecule_pose& receptor_xyz,
        const force_field_params& receptor_ffdata,
        const molecule_pose& ligand_xyz,
        const force_field_params& ligand_ffdata,
        const instant_cache* pose_memo = nullptr
    ) const override;

protected:
    class pi_stacking_batch {
    public:
        pi_stacking_batch(
            const param_t center, const PiStackingInteractionParams& params
        ) : center_(center), dist_min_(params.dist_min), dist_max_(params.dist_max),
            angle_min_(params.angle_min), angle_max_(params.angle_max),
            sigma_min_(params.sigma_min), sigma_max_(params.sigma_max) {}

        param_t calculate(
            const std::vector<atom_position>& receptor_centers,
            const std::vector<vector_3d>& receptor_normals,
            const std::vector<atom_position>& ligand_centers,
            const std::vector<vector_3d>& ligand_normals
        ) const;

    private:
        param_t center_;
        param_t dist_min_;
        param_t dist_max_;
        param_t angle_min_;
        param_t angle_max_;
        param_t sigma_min_;
        param_t sigma_max_;
    };

    param_t l5r5_coef_;
    param_t l5r6_coef_;
    param_t l6r5_coef_;
    param_t l6r6_coef_;
    std::unique_ptr<pi_stacking_batch> l5r5_inter_;
    std::unique_ptr<pi_stacking_batch> l5r6_inter_;
    std::unique_ptr<pi_stacking_batch> l6r5_inter_;
    std::unique_ptr<pi_stacking_batch> l6r6_inter_;
};

class rotatable_energy : public leaf_scorer {
public:
    rotatable_energy(const std::string& name) : leaf_scorer(name) {}

    param_t report(const molecule_pose& receptor_xyz,
                   const force_field_params& receptor_ffdata,
                   const molecule_pose& ligand_xyz,
                   const force_field_params& ligand_ffdata,
                   const instant_cache* pose_memo = nullptr) const override;
};


class gbsa_total_energy : public abstract_scorer {
public:
    gbsa_total_energy() : solvent_dielectric_(78.5),
                          solute_dielectric_(1.0),
                          offset_(0.09),
                          alpha_(1.0),
                          beta_(0.8),
                          gamma_(4.85),
                          probe_radius_(1.4),
                          surface_energy_(0.0678584012906),
                          kappa_(0.09293963686178248) {}

    virtual std::unordered_map<std::string, param_t> summary(
        const molecule_pose& receptor_xyz,
        const force_field_params& receptor_ffdata,
        const molecule_pose& ligand_xyz,
        const force_field_params& ligand_ffdata,
        const instant_cache* pose_memo = nullptr
    ) const override;

private:
    param_t calculate(const molecule_pose& coords,
                      const std::vector<param_t>& charges,
                      const std::vector<param_t>& radii,
                      const std::vector<param_t>& screens) const;

    param_t solvent_dielectric_;
    param_t solute_dielectric_;
    param_t offset_;
    param_t alpha_;
    param_t beta_;
    param_t gamma_;
    param_t probe_radius_;
    param_t surface_energy_;
    param_t kappa_;
};

}
