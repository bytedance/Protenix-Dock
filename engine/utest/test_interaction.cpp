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
#include "bytedock/ext/logging.h"

#include "test_lib.h"

namespace bytedock {

TEST(IntraMolecularTest, HarmonicBond) {
    molecule_pose mol_xyz = {
        {7.908, 27.763, 9.369},
        {8.890, 26.662, 8.959},
        {8.190, 25.303, 9.005},
        {9.472, 26.977, 7.565},
        {8.806, 26.234, 6.393},
        {9.781, 25.678, 5.445},
    };
    molecule_pose xyz_gradient(mol_xyz.size(), {0., 0., 0.});

    force_field_params mol_ffdata;
    // ids[2], k, length
    mol_ffdata.bonds = {
        {0, 1, 520.7307, 1.5364},
        {1, 2, 520.7307, 1.5364},
        {1, 3, 520.7307, 1.5364},
        {3, 4, 520.7307, 1.5364},
        {4, 5, 697.9920, 1.4652},
    };

    harmonic_bond_interaction hbi;
    double energy = hbi.put_gradients(mol_xyz, mol_ffdata, xyz_gradient);

    double ref_energy = 0.039156;
    molecule_pose ref_gradient = {
        { 1.7266, -1.9359, -0.7209},
        {-4.7169, -2.0213,  3.9922},
        { 1.6703,  3.2427, -0.1098},
        { 1.9575,  1.4257, -2.0395},
        {-2.4871,  0.3434,  0.6763},
        { 1.8495, -1.0547, -1.7983},
    };
    EXPECT_NEAR(ref_energy, energy, 1e-4);
    for (size_t i = 0; i < mol_xyz.size(); i++) {
        check_3d(ref_gradient[i].xyz, xyz_gradient[i].xyz, 2e-2);
    }
}

TEST(IntraMolecularTest, HarmonicAngle) {
    molecule_pose mol_xyz = {
        { 7.908, 27.763, 9.369},
        { 8.890, 26.662, 8.959},
        { 8.190, 25.303, 9.005},
        { 9.472, 26.977, 7.565},
        { 8.806, 26.234, 6.393},
        { 9.781, 25.678, 5.445},
        {10.482, 24.527, 5.659},
    };
    molecule_pose xyz_gradient(mol_xyz.size(), {0., 0., 0.});

    force_field_params mol_ffdata;
    // ids[3], k, radian, degree
    mol_ffdata.angles = {
        {0, 1, 2,  81.2141, MATH_TO_RADIAN(109.8749)},
        {0, 1, 3,  81.2141, MATH_TO_RADIAN(109.8749)},
        {1, 3, 4,  81.2141, MATH_TO_RADIAN(109.8749)},
        {2, 1, 3,  81.2141, MATH_TO_RADIAN(109.8749)},
        {3, 4, 5,  81.2141, MATH_TO_RADIAN(109.8749)},
        {4, 5, 6, 109.2473, MATH_TO_RADIAN(118.8684)},
    };

    harmonic_angle_interaction hai;
    double energy = hai.put_gradients(mol_xyz, mol_ffdata, xyz_gradient);

    double ref_energy = 0.971852;
    molecule_pose ref_gradient = {
        {-0.0771, -0.0219, -0.1260},
        { 3.2393,  0.5275, -0.9226},
        {-0.6397,  0.4076,  2.3082},
        {-7.9720, -2.5917,  2.9012},
        { 6.3865,  7.7524, -6.7594},
        {-4.1822, -6.8758,  8.9158},
        { 3.2451,  0.8019, -6.3172},
    };
    EXPECT_NEAR(ref_energy, energy, 2e-5);
    for (size_t i = 0; i < mol_xyz.size(); i++) {
        check_3d(ref_gradient[i].xyz, xyz_gradient[i].xyz, 2e-4);
    }
}

TEST(IntraMolecularTest, FourierDihedral) {
    molecule_pose mol_xyz = {
        { 7.908, 27.763, 9.369},
        { 8.890, 26.662, 8.959},
        { 8.190, 25.303, 9.005},
        { 9.472, 26.977, 7.565},
        { 8.806, 26.234, 6.393},
        { 9.781, 25.678, 5.445},
        {10.482, 24.527, 5.659},
        {10.255, 23.716, 6.568},
        {11.714, 24.349, 4.739},
        {11.305, 23.973, 3.332},
    };
    molecule_pose xyz_gradient(mol_xyz.size(), {0., 0., 0.});

    force_field_params mol_ffdata;
    // ids[4], periodicity[6], phase[6], k[6], phase_raw[6]
    mol_ffdata.propers = {
        {0, 1, 3, 4, 1, 2, 3, 4, 0, 0, 0., M_PI, 0., M_PI, 0., 0., -0.0237973,  0.0580221,  0.0825244,  0.0666989},
        {1, 3, 4, 5, 1, 2, 3, 4, 0, 0, 0., M_PI, 0., M_PI, 0., 0.,  0.1194276, -0.0281046,  0.1573266,  0.1047968},
        {2, 1, 3, 4, 1, 2, 3, 4, 0, 0, 0., M_PI, 0., M_PI, 0., 0., -0.0237973,  0.0580221,  0.0825244,  0.0666989},
        {3, 4, 5, 6, 1, 2, 3, 4, 0, 0, 0., M_PI, 0., M_PI, 0., 0., -0.2879912, -0.0000743, -0.0120200, -0.0986383},
        {4, 5, 6, 7, 1, 2, 3, 4, 0, 0, 0., M_PI, 0., M_PI, 0., 0., -0.5587948,  2.5157614, -0.2466097,  0.1436590},
        {4, 5, 6, 8, 1, 2, 3, 4, 0, 0, 0., M_PI, 0., M_PI, 0., 0., -0.0613361,  1.6208833,  0.1640346,  0.1171070},
        {5, 6, 8, 9, 1, 2, 3, 4, 0, 0, 0., M_PI, 0., M_PI, 0., 0.,  0.0861966, -0.0945933, -0.2458158, -0.0589235},
        {7, 6, 8, 9, 1, 2, 3, 4, 0, 0, 0., M_PI, 0., M_PI, 0., 0.,  0.0312411,  0.3406516, -0.2255443, -0.0649818},
    };
    mol_ffdata.impropers = {
        {6, 5, 7, 8, 2, 0, 0, 0, 0, 0, M_PI, 0., 0., 0., 0., 0., 3.5},
        {6, 7, 8, 5, 2, 0, 0, 0, 0, 0, M_PI, 0., 0., 0., 0., 0., 3.5},
        {6, 8, 5, 7, 2, 0, 0, 0, 0, 0, M_PI, 0., 0., 0., 0., 0., 3.5},
    };

    fourier_dihedral_interaction fdi;
    double energy = fdi.put_gradients(mol_xyz, mol_ffdata, xyz_gradient);

    double ref_energy = -0.445926;
    molecule_pose ref_gradient = {
        { 0.1888,  0.1282,  0.1078},
        {-0.1391, -0.3687, -0.1414},
        { 0.0773, -0.0390,  0.0235},
        {-0.3187,  0.7261, -0.1640},
        { 3.1714,  1.8067,  1.9171},
        {-4.4598, -3.2939, -2.4912},
        { 1.6537,  1.0094,  0.7008},
        {-1.4417, -1.1764, -1.4134},
        { 1.2504,  1.6017,  1.3607},
        { 0.0178, -0.3940,  0.1001},
    };
    EXPECT_NEAR(ref_energy, energy, 2e-6);
    for (size_t i = 0; i < mol_xyz.size(); i++) {
        check_3d(ref_gradient[i].xyz, xyz_gradient[i].xyz, 1e-4);
    }
}

TEST(IntraMolecularTest, SelfNonbonded) {
    molecule_pose mol_xyz = {
        { 7.908, 27.763, 9.369},
        { 8.890, 26.662, 8.959},
        { 8.190, 25.303, 9.005},
        { 9.472, 26.977, 7.565},
        { 8.806, 26.234, 6.393},
        { 9.781, 25.678, 5.445},
        {10.482, 24.527, 5.659},
        {10.255, 23.716, 6.568},
    };
    molecule_pose xyz_gradient(mol_xyz.size(), {0., 0., 0.});

    force_field_params mol_ffdata;
    mol_ffdata.pairs = {  // ids[2]
        {0, 4},
        {1, 5},
        {2, 4},
        {3, 6},
        {4, 7},
    };
    mol_ffdata.others = {  // ids[2]
        {0, 5},
        {0, 6},
        {0, 7},
        {1, 6},
        {1, 7},
        {2, 5},
        {2, 6},
        {2, 7},
        {3, 7},
    };
    mol_ffdata.partial_charges = {
        -0.0861, -0.0617, -0.0861, -0.0604, 0.0017, -0.5409, 0.6301, -0.6651
    };
    mol_ffdata.vdw_params = {  // sigma, epsilon
        {3.3795, 0.1088},
        {3.3795, 0.1088},
        {3.3795, 0.1088},
        {3.3795, 0.1088},
        {3.3795, 0.1088},
        {3.2069, 0.1677},
        {3.4806, 0.0869},
        {3.0398, 0.2102},
    };

    self_nonbonded_interactions sni;
    double energy = sni.put_gradients(mol_xyz, mol_ffdata, xyz_gradient);

    double ref_energy = 12.38822;
    molecule_pose ref_gradient = {
        { 0.2816, -0.4702, -0.6080},
        { 0.1016, -0.3013, -0.4143},
        { 1.9576,  1.6140, -6.5721},
        { 0.2167, -0.9959, -0.2429},
        { 0.0525, -4.3827,  5.9701},
        {-0.7626,  0.3425,  1.9301},
        { 1.2506, -1.0887, -2.0243},
        {-3.0981,  5.2822,  1.9613},
    };
    EXPECT_NEAR(ref_energy, energy, 1e-3);
    for (size_t i = 0; i < mol_xyz.size(); i++) {
        check_3d(ref_gradient[i].xyz, xyz_gradient[i].xyz, 3e-3);
    }
}

TEST(InterMolecularTest, InterMolecular) {
    molecule_pose receptor_xyz = {
        {28.80, 40.38, 5.30},
        {28.08, 40.41, 4.60}
    };
    molecule_pose receptor_gradient(receptor_xyz.size(), {0., 0., 0.});
    molecule_pose ligand_xyz = {
        { 9.781, 25.678, 5.445},
        {10.482, 24.527, 5.659},
        {10.255, 23.716, 6.568},
    };
    molecule_pose ligand_gradient(ligand_xyz.size(), {0., 0., 0.});

    force_field_params receptor_ffdata;
    receptor_ffdata.partial_charges = {
        -0.2020, 0.3120
    };
    receptor_ffdata.vdw_params = {  // sigma, epsilon
        {3.2500, 0.1700},
        {1.0691, 0.0157},
    };
    force_field_params ligand_ffdata;
    ligand_ffdata.partial_charges = {
        -0.5409, 0.6301, -0.6651
    };
    ligand_ffdata.vdw_types = {  // sigma, epsilon
        19, 13, 16,
    };

    inter_molecular_interactions imi;
    double energy = imi.put_gradients(receptor_xyz, receptor_ffdata, nullptr,
                                      ligand_xyz, ligand_ffdata,
                                      receptor_gradient, ligand_gradient);

    double ref_energy = -0.8934;
    EXPECT_NEAR(ref_energy, energy, 2e-5);
    molecule_pose ref_gradient = {
        {-0.0484, -0.0391,  0.0029},
        { 0.0765,  0.0642, -0.0077},
    };
    for (size_t i = 0; i < receptor_xyz.size(); i++) {
        check_3d(ref_gradient[i].xyz, receptor_gradient[i].xyz, 1e-4);
    }
    ref_gradient = {
        {-0.0293, -0.0252,  0.0033},
        { 0.0315,  0.0305, -0.0041},
        {-0.0303, -0.0304,  0.0056},
    };
    for (size_t i = 0; i < ligand_xyz.size(); i++) {
        check_3d(ref_gradient[i].xyz, ligand_gradient[i].xyz, 1e-4);
    }
}

TEST(InterMolecularTest, ReceptorXyzToTorison) {
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
        "num_rotatable_bonds": 5,
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

    // Prepare inputs
    boost::shared_ptr<sink_t> sink = enable_global_logging("-", 1);
    auto receptor = std::make_shared<torsional_receptor>();
    receptor->parse(rs);
    auto ligand = std::make_shared<free_ligand>();
    ligand->parse(ls);
    size_t ntorsions = receptor->num_torsions();

    // Test with theta=0
    binding_system_interactions bsi(receptor, ligand);
    std::vector<param_t> torsion_gradient;
    molecule_pose ligand_gradient;
    bsi.init_gradients(torsion_gradient, ligand_gradient);
    std::vector<param_t> receptor_torsions(ntorsions, 0.);
    molecule_pose ligand_xyz = std::move(ligand->get_pose(0));
    double energy = bsi.put_gradients(receptor_torsions, ligand_xyz,
                                      torsion_gradient, ligand_gradient);
    double ref_energy = -24.8695;
    EXPECT_NEAR(ref_energy, energy, 4e-5);
    std::vector<param_t> ref_gradient = {-1.6258,  0.1114, -0.0144};
    for (size_t i = 0; i < ntorsions; i++) {
        EXPECT_NEAR(ref_gradient[i], torsion_gradient[i], 3e-5);
    }

    // Test with different thetas
    bsi.reset_gradients(torsion_gradient, ligand_gradient);
    for (size_t i = 0; i < ntorsions; i++) {
        receptor_torsions[i] = MATH_TO_RADIAN((i*360.) / ntorsions);
    }
    energy = bsi.put_gradients(receptor_torsions, ligand_xyz,
                               torsion_gradient, ligand_gradient);
    ref_energy = -25.51645;
    EXPECT_NEAR(ref_energy, energy, 5e-5);
    ref_gradient = {-1.6258, -0.5845, -0.0027};
    for (size_t i = 0; i < ntorsions; i++) {
        EXPECT_NEAR(ref_gradient[i], torsion_gradient[i], 4e-5);
    }
    disable_sink(sink);
}

}
