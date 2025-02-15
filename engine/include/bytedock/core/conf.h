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

#include "bytedock/lib/dtype.h"

#include <string>
#include <vector>

namespace bytedock {

enum LigandForceFieldType {
    ByteFF1 = 0,  // 0.0.3-trained
    GAFF2
};

struct InteractionCoulombEnergyParams {};

struct VdwTypeScaleParams {
    param_t sigma_scale;
    param_t epsilon_scale;
};

struct InteractionVdwEnergyParams {
    param_t scale;
    param_t cutoff;
    std::vector<VdwTypeScaleParams> types;
};

struct IonicInteractionParams {
    param_t dist_min;
    param_t dist_max;
};

struct IonicEnergyParams {
    param_t lc_coef;
    param_t pc_coef;
    IonicInteractionParams lc_inter;
    IonicInteractionParams pc_inter;
};

struct CationPiInteractionParams {
    param_t dist_min;
    param_t dist_max;
};

struct PiCationicEnergyParams {
    param_t lc5_coef;
    param_t lc6_coef;
    param_t pc5_coef;
    param_t pc6_coef;
    CationPiInteractionParams lc5_inter;
    CationPiInteractionParams lc6_inter;
    CationPiInteractionParams pc5_inter;
    CationPiInteractionParams pc6_inter;
};

struct PiStackingInteractionParams {
    param_t dist_min;
    param_t dist_max;
    param_t angle_min;
    param_t angle_max;
    param_t sigma_min;
    param_t sigma_max;
};

struct PiStackingEnergyParams {  // F2F & E2F
    param_t l5r5_coef;
    param_t l5r6_coef;
    param_t l6r5_coef;
    param_t l6r6_coef;
    PiStackingInteractionParams l5r5_inter;
    PiStackingInteractionParams l5r6_inter;
    PiStackingInteractionParams l6r5_inter;
    PiStackingInteractionParams l6r6_inter;
};

struct HbondInteractionParams {
    param_t dist_min;
    param_t dist_max;
    param_t angle_min;
    param_t angle_max;
    param_t dist_center;
    param_t angle_center;
};

struct HbondEnergyParams {
    param_t scale;
    param_t hldcc_coef;
    param_t hldcn_coef;
    param_t hldnc_coef;
    param_t hldnn_coef;
    param_t hpdcc_coef;
    param_t hpdcn_coef;
    param_t hpdnc_coef;
    param_t hpdnn_coef;
    HbondInteractionParams hldcc_inter;
    HbondInteractionParams hldcn_inter;
    HbondInteractionParams hldnc_inter;
    HbondInteractionParams hldnn_inter;
    HbondInteractionParams hpdcc_inter;
    HbondInteractionParams hpdcn_inter;
    HbondInteractionParams hpdnc_inter;
    HbondInteractionParams hpdnn_inter;
};

struct HydrophobicEnergyParams {
    param_t scale;
    param_t dist_min;
    param_t dist_max;
    param_t scale_factor;
};

struct TorsionStrainPenaltyParams {};

struct RotatableEnergyParams {};

struct MolecularWeightRewardParams {};

struct GbsaEnergyDeltaParams {};

struct ScoringFunctionParams {
    LigandForceFieldType ligand_ff;

    param_t interaction_coulomb_energy_coef;
    param_t interaction_vdw_energy_coef;
    param_t ionic_energy_coef;
    param_t pi_cationic_energy_coef;
    param_t f2f_pi_stacking_energy_coef;
    param_t e2f_pi_stacking_energy_coef;
    param_t hbond_energy_coef;
    param_t hydrophobic_energy_coef;
    param_t torsion_strain_penalty_coef;
    param_t rotatable_energy_coef;
    param_t molecular_weight_reward_coef;
    param_t gbsa_energy_delta_coef;

    InteractionCoulombEnergyParams interaction_coulomb_energy;
    InteractionVdwEnergyParams interaction_vdw_energy;
    IonicEnergyParams ionic_energy;
    PiCationicEnergyParams pi_cationic_energy;
    PiStackingEnergyParams f2f_pi_stacking_energy;
    PiStackingEnergyParams e2f_pi_stacking_energy;
    HbondEnergyParams hbond_energy;
    HydrophobicEnergyParams hydrophobic_energy;
    TorsionStrainPenaltyParams torsion_strain_penalty;
    RotatableEnergyParams rotatable_energy;
    MolecularWeightRewardParams molecular_weight_reward;
    GbsaEnergyDeltaParams gbsa_energy_delta;
};

struct ScoringFunctionConfig {
    ScoringFunctionParams pose_selection;
    ScoringFunctionParams affinity_ranking;
};

std::string to_string(const InteractionCoulombEnergyParams& params,
                      const std::string& prefix = "");

std::string to_string(const InteractionVdwEnergyParams& params,
                      const std::string& prefix = "");

std::string to_string(const IonicEnergyParams& params,
                      const std::string& prefix = "");

std::string to_string(const PiCationicEnergyParams& params,
                      const std::string& prefix = "");

std::string to_string(const PiStackingEnergyParams& params,
                      const std::string& prefix = "");

std::string to_string(const HbondEnergyParams& params,
                      const std::string& prefix = "");

std::string to_string(const HydrophobicEnergyParams& params,
                      const std::string& prefix = "");

std::string to_string(const TorsionStrainPenaltyParams& params,
                      const std::string& prefix = "");

std::string to_string(const RotatableEnergyParams& params,
                      const std::string& prefix = "");

std::string to_string(const MolecularWeightRewardParams& params,
                      const std::string& prefix = "");

std::string to_string(const GbsaEnergyDeltaParams& params,
                      const std::string& prefix = "");

std::string to_string(const ScoringFunctionParams& sfp,
                      const std::string& prefix = "");

}
