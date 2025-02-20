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

#include "bytedock/core/optimizer.h"
#include "bytedock/ext/logging.h"

#include "test_lib.h"

namespace bytedock {

class OptimizerTest : public testing::Test {
protected:
    void SetUp() override {
        sink_ = enable_global_logging("-", 1);
        auto ls = std::make_shared<std::stringstream>();
        (*ls) << R"ligand_json({
    "xyz": [[
        [ 9.781, 25.678, 5.445],
        [10.482, 24.527, 5.659],
        [10.255, 23.716, 6.568]
    ]],
    "ffdata": {
        "FF_Bonds_atomidx": [],
        "FF_Bonds_k": [],
        "FF_Bonds_length": [],
        "FF_Angles_atomidx": [],
        "FF_Angles_k": [],
        "FF_Angles_angle": [],
        "FF_ProperTorsions_atomidx": [],
        "FF_ProperTorsions_periodicity": [],
        "FF_ProperTorsions_k": [],
        "FF_ProperTorsions_phase": [],
        "FF_ImproperTorsions_atomidx": [],
        "FF_ImproperTorsions_periodicity": [],
        "FF_ImproperTorsions_k": [],
        "FF_ImproperTorsions_phase": [],
        "FF_Nonbonded14_atomidx": [],
        "FF_NonbondedAll_atomidx": [],
        "hydrophobic_atomidx": [],
        "cation_atomidx": [],
        "anion_atomidx": [],
        "piring5_atomidx": [],
        "piring6_atomidx": [],
        "tstrain_atomidx": [],
        "tstrain_params_pose_selection": [],
        "hbonddon_charged_atomidx": [],
        "hbonddon_neut_atomidx": [],
        "hbondacc_charged_atomidx": [],
        "hbondacc_neut_atomidx": [],
        "partial_charges": [-0.5409, 0.6301, -0.6651],
        "FF_vdW_paraidx": [1019, 1013, 1016],
        "rotatable_bond_index": [],
        "bond_index": [],
        "atomic_numbers": [],
        "hydrophobic_group": []
    },
    "geometry": {
        "num_atoms": 3,
        "permute_index": [0, 1, 2],
        "frag_split_index": [],
        "frag_traverse_levels": [
            {
                "0": {
                    "children": [],
                    "parent": []
                }
            }
        ],
        "edge_to_bond": {}
    }
})ligand_json";
        auto rs = std::make_shared<std::stringstream>();
        (*rs) << R"receptor_json({
    "xyz": [
        [28.800, 40.380,  5.300],
        [15.410, 27.400,  6.380],
        [15.100, 27.430,  5.060],
        [14.600, 26.630,  4.840],
        [16.630, 31.420,  1.350],
        [17.020, 30.230,  0.680],
        [17.570, 30.450, -0.070],
        [ 2.680, 30.610, 11.500],
        [ 1.600, 31.410, 10.900],
        [ 1.520, 32.290, 11.390],
        [ 1.820, 31.600,  9.930],
        [ 0.730, 30.910, 10.960]
    ],
    "ffdata": {
        "FF_Bonds_atomidx": [],
        "FF_Bonds_k": [],
        "FF_Bonds_length": [],
        "FF_Angles_atomidx": [],
        "FF_Angles_k": [],
        "FF_Angles_angle": [],
        "FF_ProperTorsions_atomidx": [],
        "FF_ProperTorsions_periodicity": [],
        "FF_ProperTorsions_k": [],
        "FF_ProperTorsions_phase": [],
        "FF_ImproperTorsions_atomidx": [],
        "FF_ImproperTorsions_periodicity": [],
        "FF_ImproperTorsions_k": [],
        "FF_ImproperTorsions_phase": [],
        "FF_Nonbonded14_atomidx": [],
        "FF_NonbondedAll_atomidx": [],
        "hydrophobic_atomidx": [],
        "cation_atomidx": [],
        "anion_atomidx": [],
        "piring5_atomidx": [],
        "piring6_atomidx": [],
        "hbonddon_charged_atomidx": [],
        "hbonddon_neut_atomidx": [],
        "hbondacc_charged_atomidx": [],
        "hbondacc_neut_atomidx": [],
        "partial_charges": [
            -0.202, 0.6462, -0.6376, 0.4747, 0.3654, -0.6761,
            0.4102, -0.0143, -0.3854, 0.34, 0.34, 0.34
        ],
        "FF_vdW_sigma": [
            3.250, 3.400, 3.066, 0.000, 3.400, 3.066,
            0.000, 3.400, 3.250, 1.069, 1.069, 1.069
        ],
        "FF_vdW_epsilon": [
            0.1700, 0.0860, 0.2104, 0.0000, 0.1094, 0.2104,
            0.0000, 0.1094, 0.1700, 0.0157, 0.0157, 0.0157
        ],
        "rotatable_bond_index": [],
        "bond_index": [],
        "atomic_numbers": []
    },
    "geometry": {
        "num_atoms": 12,
        "data": {
            "2": {
                "fragments": [[2, 3], [5, 6]],
                "rot_bond_1": [1, 4],
                "rot_bond_2": [2, 5],
                "rot_bond_order": [0, 1]
            },
            "4": {
                "fragments": [[8, 9, 10, 11]],
                "rot_bond_1": [7],
                "rot_bond_2": [8],
                "rot_bond_order": [2]
            }
        }
    }
})receptor_json";
        ligand_ = std::make_shared<free_ligand>();
        ligand_->parse(ls);
        receptor_ = std::make_shared<torsional_receptor>();
        receptor_->parse(rs);
    }

    void TearDown() override {
        disable_sink(sink_);
    }

    boost::shared_ptr<sink_t> sink_;
    std::shared_ptr<free_ligand> ligand_;
    std::shared_ptr<torsional_receptor> receptor_;
};

TEST_F(OptimizerTest, ProbeZeroThetas) {
    const size_t ntorsions = receptor_->num_torsions();
    binding_system_interactions model(receptor_, ligand_);
    optimized_result out;
    out.ligand_xyz = ligand_->get_pose(0);  // Copy
    out.torsions.resize(ntorsions, 0.);
    lbfgs_step o1(1, 8);
    o1.apply(model, out);
    EXPECT_NEAR(-65.1168, out.energy, 1e-4);
    std::vector<param_t> ref_tosions = {1.7078, -0.1169, 0.0151};
    for (size_t i = 0; i < ntorsions; i++) {
        EXPECT_NEAR(ref_tosions[i], out.torsions[i], 1e-4);
    }
    molecule_pose ref_xyz = {
        {12.3282, 26.6937, 6.1659},
        { 7.2854, 22.2396, 5.1075},
        {12.6085, 26.1862, 6.4895}
    };
    for (size_t i = 0; i < out.ligand_xyz.size(); i++) {
        check_3d(ref_xyz[i].xyz, out.ligand_xyz[i].xyz, 1e-4);
    }
}

TEST_F(OptimizerTest, RunMoreIterations) {
    const size_t ntorsions = receptor_->num_torsions();
    binding_system_interactions model(receptor_, ligand_);
    optimized_result out;
    out.ligand_xyz = ligand_->get_pose(0);  // Copy
    out.torsions = {1., 2., 3.};
    lbfgs_step o2(2, 12);
    o2.apply(model, out);
    EXPECT_NEAR(-75.4551, out.energy, 1e-4);
    std::vector<param_t> ref_tosions = {0.3757, 2.6070, 2.9846};
    for (size_t i = 0; i < ntorsions; i++) {
        EXPECT_NEAR(ref_tosions[i], out.torsions[i], 2e-4);
    }
    molecule_pose ref_xyz = {
        {12.6585, 27.4573,  6.1878},
        { 7.5778, 21.4340,  5.0107},
        {12.3453, 26.6342,  6.5759}
    };
    for (size_t i = 0; i < out.ligand_xyz.size(); i++) {
        check_3d(ref_xyz[i].xyz, out.ligand_xyz[i].xyz, 2e-4);
    }
}

}
