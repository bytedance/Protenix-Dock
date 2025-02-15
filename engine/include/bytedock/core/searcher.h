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

#pragma

#include "bytedock/core/optimizer.h"
#include "bytedock/core/factory.h"
#include "bytedock/lib/random.h"

namespace bytedock {

// Each Monte-Carlo searcher only uses 1 CPU thread
class monte_carlo_searcher {
public:
    monte_carlo_searcher(
        const scoring_function_factory& sf_manager,
        const box_config& box, size_t seed, size_t max_nsteps, size_t relax_nsteps,
        double mc_mmenergy_threshold, blocking_queue<name_and_task>& in_queue
    ): full_step_(240), fast_step_(relax_nsteps), amplitude_(2_r), temperature_(1.2_r),
       sf_mgr_(sf_manager), bc_(box), seed_(seed), max_nsteps_(max_nsteps),
       mmenergy_threshold_(mc_mmenergy_threshold), task_queue_(in_queue) {
        auto& sf_config = sf_mgr_.get_config();
        vdw_calc_ = std::make_unique<mm_interaction_vdw_energy>(
            "InteractionVDW_Energy",
            sf_config.pose_selection.interaction_vdw_energy
        );
    }

    void fill(blocking_queue<name_and_task>& out_queue);

private:
    bool mutate_and_optimize(binding_input& in, optimized_result& out);
    molecule_pose randomize(const free_ligand& ligand, const size_t pose_id,
                            const torsional_receptor& receptor,
                            random_generator& rdgen);
    bool mutate_inplace(const free_ligand& ligand, random_generator& rdgen,
                        molecule_pose& ligand_xyz);

    bool metropolis_accept(param_t old_e, param_t new_e, random_generator& rdgen) {
        if (new_e < old_e) return true;
        return rdgen.uniform_real(0, 1) < std::exp((old_e - new_e) / temperature_);
    }

    bool check_all_atoms_in_box(const molecule_pose& ligand_xyz) {
        for (const auto& atom : ligand_xyz) {
            if (!check_atom_in_box(atom, bc_)) return false;
        }
        return true;
    }

    bool check_geometric_center_in_box(const molecule_pose& ligand_xyz) {
        atom_position center = calculate_geometric_center(ligand_xyz);
        return check_atom_in_box(center, bc_);
    }

    std::unique_ptr<mm_interaction_vdw_energy> vdw_calc_;
    lbfgs_step full_step_;  // For pose opitimization to local minima
    lbfgs_step fast_step_;  // For quick relaxation after pose mutation
    const param_t amplitude_;
    const param_t temperature_;

    const scoring_function_factory& sf_mgr_;
    const box_config& bc_;
    const size_t seed_;
    const size_t max_nsteps_;
    const param_t mmenergy_threshold_;
    blocking_queue<name_and_task>& task_queue_;
};

}
