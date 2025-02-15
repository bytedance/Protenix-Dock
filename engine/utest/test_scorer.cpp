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
#include "bytedock/lib/math.h"

#include "test_lib.h"

namespace bytedock {

TEST(FastScorerTest, TorsionPenalty) {
    molecule_pose mol_xyz = {
        {-19.492, 19.720, -7.201},
        {-19.063, 18.558, -7.938},
        {-18.960, 18.607, -9.386},
        {-20.077, 17.159, -7.527},
        {-21.147, 17.375, -6.637},
        {-21.982, 16.325, -6.257},
        {-21.750, 15.047, -6.749},
        {-20.675, 14.819, -7.609},
        {-19.827, 15.867, -7.998},
        {-18.539, 15.433, -9.053},
        {-17.500, 18.112, -7.433},
        {-16.455, 17.929, -8.463},
        {-15.587, 19.154, -8.657},
        {-14.229, 19.120, -8.308},
        {-13.422, 20.244, -8.487},
        {-13.959, 21.408, -9.034},
        {-15.298, 21.443, -9.417},
        {-16.107, 20.321, -9.239},
        {-17.198, 18.289, -5.976},
        {-16.256, 17.156, -5.513},
        {-15.287, 17.773, -4.563},
        {-15.163, 19.208, -4.949},
        {-16.534, 19.601, -5.524},
        {-17.337, 20.458, -4.579},
        {-17.723, 20.058, -3.207},
        {-20.147, 19.632, -3.822},
        {-17.359, 22.127, -4.945},
        {-16.891, 22.222, -6.310},
        {-18.666, 22.636, -4.594},
        {-16.156, 22.842, -3.870},
        {-16.567, 23.231, -2.585}
    };
    force_field_params mol_ffdata;
    mol_ffdata.bad_torsions.resize(4);
    mol_ffdata.bad_torsions[0].coefs[0] = 0.1;
    mol_ffdata.bad_torsions[0].coefs[1] = 0.2;
    mol_ffdata.bad_torsions[0].coefs[2] = 0.3;
    mol_ffdata.bad_torsions[0].coefs[3] = 0.4;
    mol_ffdata.bad_torsions[0].matched = {
        {4, 3, 1, 10}
    };
    mol_ffdata.bad_torsions[1].coefs[0] = 0.4;
    mol_ffdata.bad_torsions[1].coefs[1] = 0.3;
    mol_ffdata.bad_torsions[1].coefs[2] = 0.2;
    mol_ffdata.bad_torsions[1].coefs[3] = 0.1;
    mol_ffdata.bad_torsions[1].matched = {
        {30, 29, 26, 23}
    };
    mol_ffdata.bad_torsions[2].coefs[0] = 0.2;
    mol_ffdata.bad_torsions[2].coefs[1] = 0.1;
    mol_ffdata.bad_torsions[2].coefs[2] = 0.4;
    mol_ffdata.bad_torsions[2].coefs[3] = 0.3;
    mol_ffdata.bad_torsions[2].matched = {
        {3, 1, 10, 11},
        {3, 1, 10, 18}
    };
    mol_ffdata.bad_torsions[3].coefs[0] = 0.3;
    mol_ffdata.bad_torsions[3].coefs[1] = 0.4;
    mol_ffdata.bad_torsions[3].coefs[2] = 0.1;
    mol_ffdata.bad_torsions[3].coefs[3] = 0.2;
    mol_ffdata.bad_torsions[3].matched = {
        {29, 26, 23, 22},
        {29, 26, 23, 24}
    };
    torsion_strain_penalty_scorer term("TorsionStrain_Penalty");
    auto score = term.report(mol_xyz, mol_ffdata, mol_xyz, mol_ffdata);
    EXPECT_NEAR(5.07775, score, 1e-5);
}

TEST(AffinityScorerTest, RotatableEnergy) {
    molecule_pose mol_xyz;
    force_field_params mol_ffdata;
    mol_ffdata.rotatable_bonds = {
        {0, 1}, {1, 2}, {2, 3}
    };
    rotatable_energy term("rotatable_energy");
    auto score = term.report(mol_xyz, mol_ffdata, mol_xyz, mol_ffdata);
    EXPECT_NEAR(-0.1779, score, 1e-4);
}

TEST(AffinityScorerTest, GbsaTotalEnergy) {
    molecule_pose ligand_xyz = {
        { 7.9080, 27.7630,  9.3690},
        { 8.8900, 26.6620,  8.9590},
        { 8.1900, 25.3030,  9.0050},
        { 9.4720, 26.9770,  7.5650},
        { 8.8060, 26.2340,  6.3930},
        { 9.7810, 25.6780,  5.4450},
        {10.4820, 24.5270,  5.6590},
        {10.2550, 23.7160,  6.5680},
        { 7.8577, 27.8194, 10.4664},
        { 8.2507, 28.7281,  8.9676},
        { 6.9102, 27.5324,  8.9676},
        { 9.7322, 26.6217,  9.6655},
        { 8.0365, 25.0049, 10.0527},
        { 7.2166, 25.3746,  8.4977},
        { 8.8135, 24.5520,  8.4977},
        { 9.3648, 28.0573,  7.3878},
        {10.5178, 26.6366,  7.5871},
        { 8.2962, 25.3296,  6.7565},
        { 9.8020, 26.0140,  4.4710}
    };
    force_field_params ligand_ffdata;
    ligand_ffdata.partial_charges = {
        -0.0861, -0.0617, -0.0861, -0.0604, 0.0017, -0.5409, 0.6301, -0.6651, 0.0304,
         0.0304,  0.0304,  0.0117,  0.0304, 0.0304,  0.0304, 0.0442,  0.0442, 0.0477,
         0.3655
    };
    ligand_ffdata.atomic_numbers = {
        6, 6, 6, 6, 6, 7, 6, 8, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
    };
    ligand_ffdata.all_bonds = {
        {0,  1},
        {1,  2},
        {1,  3},
        {3,  4},
        {4,  5},
        {5,  6},
        {6,  7},
        {0,  8},
        {0,  9},
        {0, 10},
        {1, 11},
        {2, 12},
        {2, 13},
        {2, 14},
        {3, 15},
        {3, 16},
        {4, 17},
        {5, 18}
    };

    molecule_pose receptor_xyz;
    force_field_params receptor_ffdata;
    gbsa_total_energy term;
    auto scores = term.summary(receptor_xyz, receptor_ffdata,
                               ligand_xyz, ligand_ffdata);
    param_t refv = -14.4667, zero = 0.;
    EXPECT_NEAR(refv, scores["Gbsa_Lig_Energy"], 1e-4);
    EXPECT_NEAR(zero, scores["Gbsa_Pro_Energy"], 1e-4);
    EXPECT_NEAR(refv, scores["Gbsa_Com_Energy"], 1e-4);
    EXPECT_NEAR(zero, scores["Gbsa_Com_Energy-Gbsa_Lig_Energy-Gbsa_Pro_Energy"], 1e-4);
}

}
