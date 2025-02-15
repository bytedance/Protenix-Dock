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

#include "bytedock/core/factory.h"

#include <boost/filesystem/operations.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "bytedock/ext/logging.h"
#include "bytedock/ext/pfile.h"
#include "bytedock/lib/error.h"
#include "bytedock/lib/utility.h"

namespace bytedock {

namespace fs = boost::filesystem;
namespace pt = boost::property_tree;

inline void fill(const pt::ptree& node, InteractionCoulombEnergyParams& dest) {}

inline void fill(const pt::ptree& node, VdwTypeScaleParams& dest) {
    dest.sigma_scale = node.get<param_t>("sigma_scale");
    dest.epsilon_scale = node.get<param_t>("epsilon_scale");
}

static void fill(const pt::ptree& node, InteractionVdwEnergyParams& dest) {
    dest.scale = node.get<param_t>("scale");
    dest.cutoff = node.get<param_t>("cutoff");
    auto& types_node = node.get_child("types");
    dest.types.resize(types_node.size());
    size_t idx = 0;
    for (auto& ii : types_node) {
        fill(ii.second, dest.types[idx]);
        idx++;
    }
}

inline void fill(const pt::ptree& node, IonicInteractionParams& dest) {
    dest.dist_min = node.get<param_t>("dist_min");
    dest.dist_max = node.get<param_t>("dist_max");
}

inline void fill(const pt::ptree& node, IonicEnergyParams& dest) {
    dest.lc_coef = node.get<param_t>("lc_coef");
    dest.pc_coef = node.get<param_t>("pc_coef");
    fill(node.get_child("lc_inter"), dest.lc_inter);
    fill(node.get_child("pc_inter"), dest.pc_inter);
}

inline void fill(const pt::ptree& node, CationPiInteractionParams& dest) {
    dest.dist_min = node.get<param_t>("dist_min");
    dest.dist_max = node.get<param_t>("dist_max");
}

inline void fill(const pt::ptree& node, PiCationicEnergyParams& dest) {
    dest.lc5_coef = node.get<param_t>("lc5_coef");
    dest.lc6_coef = node.get<param_t>("lc6_coef");
    dest.pc5_coef = node.get<param_t>("pc5_coef");
    dest.pc6_coef = node.get<param_t>("pc6_coef");
    fill(node.get_child("lc5_inter"), dest.lc5_inter);
    fill(node.get_child("lc6_inter"), dest.lc6_inter);
    fill(node.get_child("pc5_inter"), dest.pc5_inter);
    fill(node.get_child("pc6_inter"), dest.pc6_inter);
}

inline void fill(const pt::ptree& node, PiStackingInteractionParams& dest) {
    dest.dist_min = node.get<param_t>("dist_min");
    dest.dist_max = node.get<param_t>("dist_max");
    dest.angle_min = node.get<param_t>("angle_min");
    dest.angle_max = node.get<param_t>("angle_max");
    dest.sigma_min = node.get<param_t>("sigma_min");
    dest.sigma_max = node.get<param_t>("sigma_max");
}

inline void fill(const pt::ptree& node, PiStackingEnergyParams& dest) {
    dest.l5r5_coef = node.get<param_t>("l5r5_coef");
    dest.l5r6_coef = node.get<param_t>("l5r6_coef");
    dest.l6r5_coef = node.get<param_t>("l6r5_coef");
    dest.l6r6_coef = node.get<param_t>("l6r6_coef");
    fill(node.get_child("l5r5_inter"), dest.l5r5_inter);
    fill(node.get_child("l5r6_inter"), dest.l5r6_inter);
    fill(node.get_child("l6r5_inter"), dest.l6r5_inter);
    fill(node.get_child("l6r6_inter"), dest.l6r6_inter);
}

inline void fill(const pt::ptree& node, HbondInteractionParams& dest) {
    dest.dist_min = node.get<param_t>("dist_min");
    dest.dist_max = node.get<param_t>("dist_max");
    dest.angle_min = node.get<param_t>("angle_min");
    dest.angle_max = node.get<param_t>("angle_max");
    dest.dist_center = node.get<param_t>("dist_center");
    dest.angle_center = node.get<param_t>("angle_center");
}

inline void fill(const pt::ptree& node, HbondEnergyParams& dest) {
    dest.scale = node.get<param_t>("scale");
    dest.hldcc_coef = node.get<param_t>("hldcc_coef");
    dest.hldcn_coef = node.get<param_t>("hldcn_coef");
    dest.hldnc_coef = node.get<param_t>("hldnc_coef");
    dest.hldnn_coef = node.get<param_t>("hldnn_coef");
    dest.hpdcc_coef = node.get<param_t>("hpdcc_coef");
    dest.hpdcn_coef = node.get<param_t>("hpdcn_coef");
    dest.hpdnc_coef = node.get<param_t>("hpdnc_coef");
    dest.hpdnn_coef = node.get<param_t>("hpdnn_coef");
    fill(node.get_child("hldcc_inter"), dest.hldcc_inter);
    fill(node.get_child("hldcn_inter"), dest.hldcn_inter);
    fill(node.get_child("hldnc_inter"), dest.hldnc_inter);
    fill(node.get_child("hldnn_inter"), dest.hldnn_inter);
    fill(node.get_child("hpdcc_inter"), dest.hpdcc_inter);
    fill(node.get_child("hpdcn_inter"), dest.hpdcn_inter);
    fill(node.get_child("hpdnc_inter"), dest.hpdnc_inter);
    fill(node.get_child("hpdnn_inter"), dest.hpdnn_inter);
}

inline void fill(const pt::ptree& node, HydrophobicEnergyParams& dest) {
    dest.scale = node.get<param_t>("scale");
    dest.dist_min = node.get<param_t>("dist_min");
    dest.dist_max = node.get<param_t>("dist_max");
    dest.scale_factor = node.get<param_t>("scale_factor");
}

inline void fill(const pt::ptree& node, TorsionStrainPenaltyParams& dest) {}

inline void fill(const pt::ptree& node, RotatableEnergyParams& dest) {}

inline void fill(const pt::ptree& node, MolecularWeightRewardParams& dest) {}

inline void fill(const pt::ptree& node, GbsaEnergyDeltaParams& dest) {}

inline void fill(const pt::ptree& node, ScoringFunctionParams& sfp) {
    sfp.ligand_ff = static_cast<LigandForceFieldType>(node.get<index_t>("ligand_ff"));
    sfp.interaction_coulomb_energy_coef = 
        node.get<param_t>("interaction_coulomb_energy_coef", 0_r);
    sfp.interaction_vdw_energy_coef = 
        node.get<param_t>("interaction_vdw_energy_coef", 0_r);
    sfp.ionic_energy_coef = node.get<param_t>("ionic_energy_coef", 0_r);
    sfp.pi_cationic_energy_coef = node.get<param_t>("pi_cationic_energy_coef", 0_r);
    sfp.f2f_pi_stacking_energy_coef = 
        node.get<param_t>("f2f_pi_stacking_energy_coef", 0_r);
    sfp.e2f_pi_stacking_energy_coef = 
        node.get<param_t>("e2f_pi_stacking_energy_coef", 0_r);
    sfp.hbond_energy_coef = node.get<param_t>("hbond_energy_coef", 0_r);
    sfp.hydrophobic_energy_coef = node.get<param_t>("hydrophobic_energy_coef", 0_r);
    sfp.torsion_strain_penalty_coef = 
        node.get<param_t>("torsion_strain_penalty_coef", 0_r);
    sfp.rotatable_energy_coef = node.get<param_t>("rotatable_energy_coef", 0_r);
    sfp.molecular_weight_reward_coef = 
        node.get<param_t>("molecular_weight_reward_coef", 0_r);
    sfp.gbsa_energy_delta_coef = node.get<param_t>("gbsa_energy_delta_coef", 0_r);
    if CHECK_PARAMETER_NON_ZERO(sfp.interaction_coulomb_energy_coef) {
        fill(node.get_child("interaction_coulomb_energy"),
             sfp.interaction_coulomb_energy);
    }
    if CHECK_PARAMETER_NON_ZERO(sfp.interaction_vdw_energy_coef) {
        fill(node.get_child("interaction_vdw_energy"),
             sfp.interaction_vdw_energy);
    }
    if CHECK_PARAMETER_NON_ZERO(sfp.ionic_energy_coef) {
        fill(node.get_child("ionic_energy"), sfp.ionic_energy);
    }
    if CHECK_PARAMETER_NON_ZERO(sfp.pi_cationic_energy_coef) {
        fill(node.get_child("pi_cationic_energy"), sfp.pi_cationic_energy);
    }
    if CHECK_PARAMETER_NON_ZERO(sfp.f2f_pi_stacking_energy_coef) {
        fill(node.get_child("f2f_pi_stacking_energy"), sfp.f2f_pi_stacking_energy);
    }
    if CHECK_PARAMETER_NON_ZERO(sfp.e2f_pi_stacking_energy_coef) {
        fill(node.get_child("e2f_pi_stacking_energy"), sfp.e2f_pi_stacking_energy);
    }
    if CHECK_PARAMETER_NON_ZERO(sfp.hbond_energy_coef) {
        fill(node.get_child("hbond_energy"), sfp.hbond_energy);
    }
    if CHECK_PARAMETER_NON_ZERO(sfp.hydrophobic_energy_coef) {
        fill(node.get_child("hydrophobic_energy"), sfp.hydrophobic_energy);
    }
    if CHECK_PARAMETER_NON_ZERO(sfp.torsion_strain_penalty_coef) {
        fill(node.get_child("torsion_strain_penalty"), sfp.torsion_strain_penalty);
    }
    if CHECK_PARAMETER_NON_ZERO(sfp.rotatable_energy_coef) {
        fill(node.get_child("rotatable_energy"), sfp.rotatable_energy);
    }
    if CHECK_PARAMETER_NON_ZERO(sfp.molecular_weight_reward_coef) {
        fill(node.get_child("molecular_weight_reward"), sfp.molecular_weight_reward);
    }
    if CHECK_PARAMETER_NON_ZERO(sfp.gbsa_energy_delta_coef) {
        fill(node.get_child("gbsa_energy_delta"), sfp.gbsa_energy_delta);
    }
}

static void parse_scoring_function_config(
    std::shared_ptr<std::istream> fin, ScoringFunctionConfig& sfc
) {
    pt::ptree node;
    pt::read_json(*fin, node);
    fill(node.get_child("pose_selection"), sfc.pose_selection);
    fill(node.get_child("affinity_ranking"), sfc.affinity_ranking);
}

scoring_function_factory::scoring_function_factory(const std::string& config_path) {
    std::shared_ptr<std::istream> fin;
    if (config_path.empty()) {
        // exec path -> "bin" dir -> root dir -> "conf" dir -> param file
        fs::path exe_path(get_executable_path());
        auto param_file = exe_path.parent_path().parent_path()
                          / "conf" / "pscore-v7_and_bscore-fake.json";
        if (!fs::exists(param_file))  throw failed_conf_error(
            "Failed to find the default parameter file of scoring function. Please "
            "explicitly provide its location!."
        );
        auto fall_back = param_file.string();
        LOG_INFO << "Since no scoring function parameter file provided, the default "
                 << "one is enabled at: " << fall_back;
        fin = open_for_read(fall_back);
    } else {
        LOG_INFO << "Intializing scoring function with parameters in: " << config_path;
        fin = open_for_read(config_path);
    }
    parse_scoring_function_config(fin, sf_config_);
    LOG_DEBUG << "pscore = " << to_string(sf_config_.pose_selection);
    LOG_DEBUG << "bscore = " << to_string(sf_config_.affinity_ranking);
    pose_scorer_ = std::make_shared<root_scorer>("pscore", sf_config_.pose_selection);
    affinity_scorer_ = std::make_shared<root_scorer>("bscore",
                                                     sf_config_.affinity_ranking);
}

}
