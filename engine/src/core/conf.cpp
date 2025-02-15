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

#include "bytedock/core/conf.h"

#include "bytedock/lib/utility.h"

#include <sstream>

namespace bytedock {

std::string to_string(const InteractionCoulombEnergyParams& params,
                      const std::string& prefix) {
    return "{}";
}

static std::string to_string(const VdwTypeScaleParams& params,
                             const std::string& prefix) {
    std::ostringstream oss;
    oss << "{" << std::endl;
    oss << prefix << "  \"sigma_scale\": " << params.sigma_scale << "," << std::endl;
    oss << prefix << "  \"epsilon_scale\": " << params.epsilon_scale << std::endl;
    oss << prefix << "}";
    return oss.str();
}

std::string to_string(const InteractionVdwEnergyParams& params,
                      const std::string& prefix) {
    std::ostringstream oss;
    std::string indent(prefix + "    ");
    oss << "{" << std::endl;
    oss << prefix << "  \"scale\": "
        << params.scale << "," << std::endl;
    oss << prefix << "  \"cutoff\": "
        << params.cutoff << "," << std::endl;
    oss << prefix << "  \"types\": [";
    size_t i;
    for (i = 0; i < params.types.size() - 1; ++i) {
        oss << indent << to_string(params.types[i], indent)
            << "," << std::endl;
    }
    oss << indent << to_string(params.types[i], indent) << std::endl;
    oss << prefix << "  ]" << std::endl;
    oss << prefix << "}";
    return oss.str();
}

static std::string to_string(const IonicInteractionParams& params,
                             const std::string& prefix) {
    std::ostringstream oss;
    oss << "{" << std::endl;
    oss << prefix << "  \"dist_min\": " << params.dist_min << "," << std::endl;
    oss << prefix << "  \"dist_max\": " << params.dist_max << std::endl;
    oss << prefix << "}";
    return oss.str();
}

std::string to_string(const IonicEnergyParams& params,
                      const std::string& prefix) {
    std::ostringstream oss;
    std::string indent(prefix + "  ");
    oss << "{" << std::endl;
    oss << prefix << "  \"lc_coef\": "
        << params.lc_coef << "," << std::endl;
    oss << prefix << "  \"pc_coef\": "
        << params.pc_coef << "," << std::endl;
    oss << prefix << "  \"lc_inter\": " << to_string(params.lc_inter, indent)
        << "," << std::endl;
    oss << prefix << "  \"pc_inter\": " << to_string(params.pc_inter, indent)
        << std::endl;
    oss << prefix << "}";
    return oss.str();
}

static std::string to_string(const CationPiInteractionParams& params,
                             const std::string& prefix) {
    std::ostringstream oss;
    oss << "{" << std::endl;
    oss << prefix << "  \"dist_min\": " << params.dist_min << "," << std::endl;
    oss << prefix << "  \"dist_max\": " << params.dist_max << std::endl;
    oss << prefix << "}";
    return oss.str();
}

std::string to_string(const PiCationicEnergyParams& params,
                      const std::string& prefix) {
    std::ostringstream oss;
    std::string indent(prefix + "  ");
    oss << "{" << std::endl;
    oss << prefix << "  \"lc5_coef\": "
        << params.lc5_coef << "," << std::endl;
    oss << prefix << "  \"lc6_coef\": "
        << params.lc6_coef << "," << std::endl;
    oss << prefix << "  \"pc5_coef\": "
        << params.pc5_coef << "," << std::endl;
    oss << prefix << "  \"pc6_coef\": "
        << params.pc6_coef << "," << std::endl;
    oss << prefix << "  \"lc5_inter\": " << to_string(params.lc5_inter, indent)
        << "," << std::endl;
    oss << prefix << "  \"lc6_inter\": " << to_string(params.lc6_inter, indent)
        << "," << std::endl;
    oss << prefix << "  \"pc5_inter\": " << to_string(params.pc5_inter, indent)
        << "," << std::endl;
    oss << prefix << "  \"pc6_inter\": " << to_string(params.pc6_inter, indent)
        << std::endl;
    oss << prefix << "}";
    return oss.str();
}

static std::string to_string(const PiStackingInteractionParams& params,
                             const std::string& prefix) {
    std::ostringstream oss;
    oss << "{" << std::endl;
    oss << prefix << "  \"dist_min\": " << params.dist_min << "," << std::endl;
    oss << prefix << "  \"dist_max\": " << params.dist_max << "," << std::endl;
    oss << prefix << "  \"angle_min\": " << params.angle_min << "," << std::endl;
    oss << prefix << "  \"angle_max\": " << params.angle_max << "," << std::endl;
    oss << prefix << "  \"sigma_min\": " << params.sigma_min << "," << std::endl;
    oss << prefix << "  \"sigma_max\": " << params.sigma_max << std::endl;
    oss << prefix << "}";
    return oss.str();
}

std::string to_string(const PiStackingEnergyParams& params,
                      const std::string& prefix) {
    std::ostringstream oss;
    std::string indent(prefix + "  ");
    oss << "{" << std::endl;
    oss << prefix << "  \"l5r5_coef\": "
        << params.l5r5_coef << "," << std::endl;
    oss << prefix << "  \"l5r6_coef\": "
        << params.l5r6_coef << "," << std::endl;
    oss << prefix << "  \"l6r5_coef\": "
        << params.l6r5_coef << "," << std::endl;
    oss << prefix << "  \"l6r6_coef\": "
        << params.l6r6_coef << "," << std::endl;
    oss << prefix << "  \"l5r5_inter\": " << to_string(params.l5r5_inter, indent)
        << "," << std::endl;
    oss << prefix << "  \"l5r6_inter\": " << to_string(params.l5r6_inter, indent)
        << "," << std::endl;
    oss << prefix << "  \"l6r5_inter\": " << to_string(params.l6r5_inter, indent)
        << "," << std::endl;
    oss << prefix << "  \"l6r6_inter\": " << to_string(params.l6r6_inter, indent)
        << std::endl;
    oss << prefix << "}";
    return oss.str();
}

static std::string to_string(const HbondInteractionParams& params,
                             const std::string& prefix) {
    std::ostringstream oss;
    oss << "{" << std::endl;
    oss << prefix << "  \"dist_min\": " << params.dist_min << "," << std::endl;
    oss << prefix << "  \"dist_max\": " << params.dist_max << "," << std::endl;
    oss << prefix << "  \"angle_min\": " << params.angle_min << "," << std::endl;
    oss << prefix << "  \"angle_max\": " << params.angle_max << "," << std::endl;
    oss << prefix << "  \"dist_center\": " << params.dist_center << "," << std::endl;
    oss << prefix << "  \"angle_center\": " << params.angle_center << std::endl;
    oss << prefix << "}";
    return oss.str();
}

std::string to_string(const HbondEnergyParams& params,
                      const std::string& prefix) {
    std::ostringstream oss;
    std::string indent(prefix + "  ");
    oss << "{" << std::endl;
    oss << prefix << "  \"scale\": "
        << params.scale << "," << std::endl;
    oss << prefix << "  \"hldcc_coef\": "
        << params.hldcc_coef << "," << std::endl;
    oss << prefix << "  \"hldcn_coef\": "
        << params.hldcn_coef << "," << std::endl;
    oss << prefix << "  \"hldnc_coef\": "
        << params.hldnc_coef << "," << std::endl;
    oss << prefix << "  \"hldnn_coef\": "
        << params.hldnn_coef << "," << std::endl;
    oss << prefix << "  \"hpdcc_coef\": "
        << params.hpdcc_coef << "," << std::endl;
    oss << prefix << "  \"hpdcn_coef\": "
        << params.hpdcn_coef << "," << std::endl;
    oss << prefix << "  \"hpdnc_coef\": "
        << params.hpdnc_coef << "," << std::endl;
    oss << prefix << "  \"hpdnn_coef\": "
        << params.hpdnn_coef << "," << std::endl;
    oss << prefix << "  \"hldcc_inter\": " << to_string(params.hldcc_inter, indent)
        << "," << std::endl;
    oss << prefix << "  \"hldcn_inter\": " << to_string(params.hldcn_inter, indent)
        << "," << std::endl;
    oss << prefix << "  \"hldnc_inter\": " << to_string(params.hldnc_inter, indent)
        << "," << std::endl;
    oss << prefix << "  \"hldnn_inter\": " << to_string(params.hldnn_inter, indent)
        << "," << std::endl;
    oss << prefix << "  \"hpdcc_inter\": " << to_string(params.hpdcc_inter, indent)
        << "," << std::endl;
    oss << prefix << "  \"hpdcn_inter\": " << to_string(params.hpdcn_inter, indent)
        << "," << std::endl;
    oss << prefix << "  \"hpdnc_inter\": " << to_string(params.hpdnc_inter, indent)
        << "," << std::endl;
    oss << prefix << "  \"hpdnn_inter\": " << to_string(params.hpdnn_inter, indent)
        << std::endl;
    oss << prefix << "}";
    return oss.str();
}

std::string to_string(const HydrophobicEnergyParams& params,
                      const std::string& prefix) {
    std::ostringstream oss;
    oss << "{" << std::endl;
    oss << prefix << "  \"scale\": " << params.scale << "," << std::endl;
    oss << prefix << "  \"dist_min\": " << params.dist_min << "," << std::endl;
    oss << prefix << "  \"dist_max\": " << params.dist_max << "," << std::endl;
    oss << prefix << "  \"scale_factor\": " << params.scale_factor << std::endl;
    oss << prefix << "}";
    return oss.str();
}

std::string to_string(const TorsionStrainPenaltyParams& params,
                      const std::string& prefix) {
    return "{}";
}

std::string to_string(const RotatableEnergyParams& params,
                      const std::string& prefix) {
    return "{}";
}

std::string to_string(const MolecularWeightRewardParams& params,
                      const std::string& prefix) {
    return "{}";
}

std::string to_string(const GbsaEnergyDeltaParams& params,
                      const std::string& prefix) {
    return "{}";
}

std::string to_string(const ScoringFunctionParams& params,
                      const std::string& prefix) {
    std::ostringstream oss;
    std::string indent(prefix + "  ");
    oss << "{" << std::endl;
    oss << prefix << "  \"interaction_coulomb_energy_coef\": "
        << params.interaction_coulomb_energy_coef << "," << std::endl;
    if CHECK_PARAMETER_NON_ZERO(params.interaction_coulomb_energy_coef) {
        oss << prefix << "  \"interaction_coulomb_energy\": "
            << to_string(params.interaction_coulomb_energy, indent)
            << "," << std::endl;
    }
    oss << prefix << "  \"interaction_vdw_energy_coef\": "
        << params.interaction_vdw_energy_coef << "," << std::endl;
    if CHECK_PARAMETER_NON_ZERO(params.interaction_vdw_energy_coef) {
        oss << prefix << "  \"interaction_vdw_energy\": "
            << to_string(params.interaction_vdw_energy, indent)
            << "," << std::endl;
    }
    oss << prefix << "  \"ionic_energy_coef\": "
        << params.ionic_energy_coef << "," << std::endl;
    if CHECK_PARAMETER_NON_ZERO(params.ionic_energy_coef) {
        oss << prefix << "  \"ionic_energy\": "
            << to_string(params.ionic_energy, indent) << "," << std::endl;
    }
    oss << prefix << "  \"pi_cationic_energy_coef\": "
        << params.pi_cationic_energy_coef << "," << std::endl;
    if CHECK_PARAMETER_NON_ZERO(params.pi_cationic_energy_coef) {
        oss << prefix << "  \"pi_cationic_energy\": "
            << to_string(params.pi_cationic_energy, indent) << "," << std::endl;
    }
    oss << prefix << "  \"f2f_pi_stacking_energy_coef\": "
        << params.f2f_pi_stacking_energy_coef << "," << std::endl;
    if CHECK_PARAMETER_NON_ZERO(params.f2f_pi_stacking_energy_coef) {
        oss << prefix << "  \"f2f_pi_stacking_energy\": "
            << to_string(params.f2f_pi_stacking_energy, indent) << "," << std::endl;
    }
    oss << prefix << "  \"e2f_pi_stacking_energy_coef\": "
        << params.e2f_pi_stacking_energy_coef << "," << std::endl;
    if CHECK_PARAMETER_NON_ZERO(params.e2f_pi_stacking_energy_coef) {
        oss << prefix << "  \"e2f_pi_stacking_energy\": "
            << to_string(params.e2f_pi_stacking_energy, indent) << "," << std::endl;
    }
    oss << prefix << "  \"hbond_energy_coef\": " << params.hbond_energy_coef
        << "," << std::endl;
    if CHECK_PARAMETER_NON_ZERO(params.hbond_energy_coef) {
        oss << prefix << "  \"hbond_energy\": "
            << to_string(params.hbond_energy, indent) << "," << std::endl;
    }
    oss << prefix << "  \"hydrophobic_energy_coef\": "
        << params.hydrophobic_energy_coef << "," << std::endl;
    if CHECK_PARAMETER_NON_ZERO(params.hydrophobic_energy_coef) {
        oss << prefix << "  \"hydrophobic_energy\": "
            << to_string(params.hydrophobic_energy, indent) << "," << std::endl;
    }
    oss << prefix << "  \"torsion_strain_penalty_coef\": "
        << params.torsion_strain_penalty_coef << "," << std::endl;
    if CHECK_PARAMETER_NON_ZERO(params.torsion_strain_penalty_coef) {
        oss << prefix << "  \"torsion_strain_penalty\": "
            << to_string(params.torsion_strain_penalty, indent) << "," << std::endl;
    }
    oss << prefix << "  \"rotatable_energy_coef\": "
        << params.rotatable_energy_coef << "," << std::endl;
    if CHECK_PARAMETER_NON_ZERO(params.rotatable_energy_coef) {
        oss << prefix << "  \"rotatable_energy\": "
            << to_string(params.rotatable_energy, indent) << "," << std::endl;
    }
    oss << prefix << "  \"molecular_weight_reward_coef\": "
        << params.molecular_weight_reward_coef << "," << std::endl;
    if CHECK_PARAMETER_NON_ZERO(params.molecular_weight_reward_coef) {
        oss << prefix << "  \"molecular_weight_reward\": "
            << to_string(params.molecular_weight_reward, indent) << "," << std::endl;
    }
    oss << prefix << "  \"gbsa_energy_delta_coef\": "
        << params.gbsa_energy_delta_coef << "," << std::endl;
    if CHECK_PARAMETER_NON_ZERO(params.gbsa_energy_delta_coef) {
        oss << prefix << "  \"gbsa_energy_delta\": "
            << to_string(params.gbsa_energy_delta, indent) << "," << std::endl;
    }
    oss << prefix << "  \"ligand_ff\": " << params.ligand_ff << std::endl;
    oss << prefix << "}";
    return oss.str();
}

}
