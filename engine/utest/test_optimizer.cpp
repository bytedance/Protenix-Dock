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
        "tstrain_params": [],
        "hbonddon_charged_atomidx": [],
        "hbonddon_neut_atomidx": [],
        "hbondacc_charged_atomidx": [],
        "hbondacc_neut_atomidx": [],
        "partial_charges": [-0.5409, 0.6301, -0.6651],
        "FF_vdW_paraidx": [19, 13, 16],
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
    EXPECT_NEAR(-66.2703, out.energy, 1e-4);
    std::vector<param_t> ref_tosions = {1.7399, -0.1192, 0.0154};
    for (size_t i = 0; i < ntorsions; i++) {
        EXPECT_NEAR(ref_tosions[i], out.torsions[i], 1e-4);
    }
    molecule_pose ref_xyz = {
        {12.3818, 26.7147, 6.1795},
        { 7.2128, 22.1888, 5.0975},
        {12.6529, 26.2327, 6.4881}
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
    EXPECT_NEAR(-75.9029, out.energy, 1e-4);
    std::vector<param_t> ref_tosions = {0.3791, 2.6009, 2.9846};
    for (size_t i = 0; i < ntorsions; i++) {
        EXPECT_NEAR(ref_tosions[i], out.torsions[i], 2e-4);
    }
    molecule_pose ref_xyz = {
        {12.6351, 27.4424,  6.1808},
        { 7.5941, 21.4569,  5.0175},
        {12.3234, 26.6055,  6.5765}
    };
    for (size_t i = 0; i < out.ligand_xyz.size(); i++) {
        check_3d(ref_xyz[i].xyz, out.ligand_xyz[i].xyz, 2e-4);
    }
}

}
