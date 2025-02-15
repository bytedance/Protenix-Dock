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

#include "bytedock/core/searcher.h"

#include "bytedock/ext/logging.h"
#include "bytedock/lib/math.h"
#include "bytedock/lib/printer.h"

namespace bytedock {

// init_max_energy 1e8

std::mutex energy_info_mutex;

void monte_carlo_searcher::fill(blocking_queue<name_and_task>& out_queue) {
    name_and_task pair = std::move(task_queue_.pop());
    while (!task_queue_.is_eoq(pair)) {
        LOG_INFO << "Task [" << pair.first << "] has started...";
        {
            step_timer t(SEARCH_BINDING_POSE_STEP);
            pair.second.fetch.name = remove_pose_suffix_value(pair.first);
            mutate_and_optimize(pair.second.feed, pair.second.fetch);
        }
        out_queue.push(std::move(pair));
        pair = std::move(task_queue_.pop());
    }
    task_queue_.close();  // Put back the EOF flag
    out_queue.close();
}

molecule_pose monte_carlo_searcher::randomize(
    const free_ligand& ligand, const size_t pose_id,
    const torsional_receptor& receptor, random_generator& rdgen
) {
    atom_position new_center;
    vector_4d quaternion;
    matrix_3x3 orientation;
    std::vector<param_t> torsions(ligand.num_torsions(), 0.0_r);
    molecule_pose candidate;
    auto& receptor_xyz = receptor.get_positions();
    auto& receptor_ffdata = receptor.get_ffdata();
    auto& ligand_xyz = ligand.get_pose(pose_id);
    auto& ligand_ffdata = ligand.get_ffdata();
    size_t num = 0;
    while (true) {
        rdgen.random_in_box(bc_.init_xyz, bc_.oppo_xyz, new_center.xyz);
        rdgen.random_orientation(quaternion.abcd);
        for (size_t i = 0; i < torsions.size(); ++i) {
            if (pose_id == 0) {
                torsions[i] = rdgen.uniform_torsion();
            }

        }
        get_rotation_matrix(quaternion, orientation);
        candidate = ligand.apply_parameters(ligand_xyz, new_center,
                                            orientation, torsions);  // Move
        ++num;
        if (num < 2000 && !check_all_atoms_in_box(candidate)) continue;
        if (!check_geometric_center_in_box(candidate)) continue;
        if (num > 1999 ||
            vdw_calc_->report(receptor_xyz, receptor_ffdata,
                              candidate, ligand_ffdata) < 100_r) break;
    }
    LOG_DEBUG << "New center for intial conformer: [" << new_center.xyz[0]
            << ", " << new_center.xyz[1] << ", " << new_center.xyz[2] << "]";
    LOG_DEBUG << "Orientation for intial conformer: [" << quaternion.abcd[0]
            << ", " << quaternion.abcd[1] << ", " << quaternion.abcd[2]
            << ", " << quaternion.abcd[3] << "]";
    LOG_DEBUG << "Torsions for intial conformer: [" << std::endl
              << format<param_t, 10>(torsions) << "]";
    return candidate;
}

inline bool check_ligand_xyz(const molecule_pose& coords) {
    for (auto& atom : coords) {
        if (atom.xyz[0] > 1e8_r || atom.xyz[0] < -1e8_r ||
            atom.xyz[1] > 1e8_r || atom.xyz[1] < -1e8_r ||
            atom.xyz[2] > 1e8_r || atom.xyz[2] < -1e8_r) return false;
    }
    return true;
}

bool monte_carlo_searcher::mutate_and_optimize(binding_input& in, optimized_result& out) {
    random_generator rg(seed_ + in.pose_id);
    out.offset = in.pose_id;
    auto& receptor_ffdata = in.receptor->get_ffdata();
    auto& ligand_ffdata = in.ligand->get_ffdata();
    std::string current_ligand_name = out.name;
    // Prepare initial conformer
    if (in.pose_id == 0) {
        out.ligand_xyz = in.ligand->get_pose(0);  // Copy
    } else {
        size_t current_candidate_idx = 0;
        if (in.ligand->get_nposes() > in.pose_id) {
            current_candidate_idx = in.pose_id;
            LOG_WARNING << "Monte-Carlo searcher randomize init pose from pose_id:" << current_candidate_idx;
        }
        out.ligand_xyz = randomize(*in.ligand, current_candidate_idx, *in.receptor, rg);  // Move

    }
    out.torsions.resize(in.receptor->num_torsions(), 0_r);
    binding_system_interactions model(in.receptor, in.ligand, in.cache);
    bool converged = full_step_.apply(model, out);
    out.receptor_xyz = in.receptor->apply_parameters(out.torsions);  // Move
    auto& pp = sf_mgr_.get_pose_selection();
    out.pscore = pp.combine(out.receptor_xyz, receptor_ffdata,
                            out.ligand_xyz, ligand_ffdata);  // `out` is the best
    if (in.pose_id == 0) return converged;  // Walker 0 disables Monte-Carlo searching

    std::vector<param_t> tmp_torsion_gradient;
    molecule_pose tmp_ligand_gradient;
    model.init_gradients(tmp_torsion_gradient, tmp_ligand_gradient);

    // Walk in conformer space
    if (check_all_atoms_in_box(out.ligand_xyz)) {
        optimized_result tmp(out);  // Copy
        optimized_result cdd;
        index_t total_nevals = out.nevals;
        for (size_t step = 0; step < max_nsteps_; ++step) {
            cdd = tmp;  // Copy
            if (!mutate_inplace(*in.ligand, rg, cdd.ligand_xyz)) {
                LOG_WARNING << "Monte-Carlo searcher#" << in.pose_id << " skipped "
                            << "orientation mutation at step#" << step << " due to too "
                            << "small gyration radius!";
            }
            if(mmenergy_threshold_ > 0.0_r) {
                param_t tmp_energy = 0.0_r;
                {
                    // Evaluate initial f(x) and df/dx, but only use f(x)
                    tmp_energy = model.put_gradients(cdd.torsions, cdd.ligand_xyz,
                                                    tmp_torsion_gradient, tmp_ligand_gradient);
                    model.reset_gradients(tmp_torsion_gradient, tmp_ligand_gradient);
                }
                if(tmp_energy < global_energy_info[current_ligand_name]["minimum_pre_opt_energy"]) {
                    LOG_WARNING << "Monte-Carlo searcher update minimum_pre_opt_energy to :" << tmp_energy;
                    std::lock_guard<std::mutex> lock(energy_info_mutex);
                    global_energy_info[current_ligand_name]["minimum_pre_opt_energy"] = tmp_energy;
                }
                if(tmp_energy > global_energy_info[current_ligand_name]["minimum_pre_opt_energy"] + mmenergy_threshold_) {
                    continue;
                }
            }

            converged = fast_step_.apply(model, cdd);
            total_nevals += cdd.nevals;
            cdd.receptor_xyz = in.receptor->apply_parameters(cdd.torsions);  // Move
            cdd.pscore = pp.combine(cdd.receptor_xyz, receptor_ffdata,
                                    cdd.ligand_xyz, ligand_ffdata);
            if ((step == 0 || metropolis_accept(tmp.pscore, cdd.pscore, rg)) &&
                check_all_atoms_in_box(cdd.ligand_xyz)) {
                tmp = std::move(cdd);
                if (tmp.pscore < out.pscore) {
                    converged = full_step_.apply(model, tmp);
                    total_nevals += tmp.nevals;
                    tmp.receptor_xyz = in.receptor->apply_parameters(tmp.torsions);
                    tmp.pscore = pp.combine(tmp.receptor_xyz, receptor_ffdata,
                                            tmp.ligand_xyz, ligand_ffdata);
                    if (tmp.pscore < out.pscore) out = tmp;  // Copy
                }
            }
        }
        out.nevals = total_nevals;
    } else if (!check_ligand_xyz(out.ligand_xyz)) {
        LOG_WARNING << "Monte-Carlo searcher#" << in.pose_id << " crashed due to "
                    << "too large coordinates in the relaxed initial conformer!";
        out.nevals = 0;  // Tell any writer to skip this conformer
    } else {
        LOG_WARNING << "Monte-Carlo searcher#" << in.pose_id << " stopped due to "
                    << "out-of-box atoms in the relaxed initial conformer!";
    }
    return converged;
}

bool monte_carlo_searcher::mutate_inplace(const free_ligand& ligand,
                                          random_generator& rdgen,
                                          molecule_pose& ligand_xyz) {
    int which = rdgen.uniform_int(0, ligand.num_torsions() + 1);
    vector_3d delta;

    // Position
    atom_position center = calculate_geometric_center(ligand_xyz);
    if (which == 0) {
        rdgen.random_in_sphere(delta.xyz);
        scale_inplace_3d(delta.xyz, amplitude_);
        add_inplace_3d(center.xyz, delta.xyz);
    }

    // Orientation
    matrix_3x3 orientation = {1_r, 0_r, 0_r, 0_r, 1_r, 0_r, 0_r, 0_r, 1_r};
    if (which == 1) {
        param_t gr = calculate_gyration_radius(ligand_xyz, center);
        if (gr < 1e-6) return false;  // or kParamEpsilon
        rdgen.random_in_sphere(delta.xyz);
        scale_inplace_3d(delta.xyz, amplitude_ / gr);
        gr = get_norm_3d(delta.xyz);  // Used as `angle`
        scale_inplace_3d(delta.xyz, 1_r / gr);  // Used as `axis`
        get_rotation_matrix(delta, gr, orientation);
    }

    // Torsions
    std::vector<param_t> torsions(ligand.num_torsions(), 0_r);
    if (which > 1) torsions[which - 2] = rdgen.uniform_torsion();

    // Update ligand conformer
    auto new_xyz = ligand.apply_parameters(ligand_xyz, center, orientation, torsions);
    ligand_xyz.swap(new_xyz);
    return true;
}

}
