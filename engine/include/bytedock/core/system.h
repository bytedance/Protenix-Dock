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

#include <string>
#include <vector>

#include "bytedock/lib/dtype.h"

namespace bytedock {

/**
 * TODO(shishaochen): Separate force field parameters from interaction classes. Then,
 * precalculation of pairwise items can be enabled.
 */

struct harmonic_bond {
    index_t ids[2];  // i, j

    param_t k;
    param_t length;
};

struct harmonic_angle {
    index_t ids[3];  // i, j, k

    param_t k;
    param_t radian;  // In program
    param_t degree;  // In file
};

// Copy from `byteff.parametrize.ffparams.baseffparams.MAX_PERIODICITY_PROPER`
#define FF_DIHEDRAL_NDIMS 6

struct fourier_dihedral {
    index_t ids[4];  // i, j, k, l

    // For a dimension `dd`, involve its contribution only when `periodicity[dd]>0`
    index_t periodicity[FF_DIHEDRAL_NDIMS];
    param_t phase[FF_DIHEDRAL_NDIMS];  // Radian in program
    param_t k[FF_DIHEDRAL_NDIMS];
    param_t phase_raw[FF_DIHEDRAL_NDIMS];  // Degree in file
};

struct atom_pair {
    index_t ids[2];
};

struct lj_vdw {
    param_t sigma;
    param_t epsilon;
};

struct pi_ring {
    index_t ids[6];
    index_t size;
};

struct torsion_index {
    index_t ids[4];
};

struct torsion_strain_penalty {
    param_t coefs[4];
    std::vector<torsion_index> matched;
};

const param_t kCoulombFactor = 332.0637141491396;  // For pairwise Coulumb interaction

struct force_field_params {
    std::vector<harmonic_bond> bonds;
    std::vector<harmonic_angle> angles;
    std::vector<fourier_dihedral> impropers;
    std::vector<fourier_dihedral> propers;
    std::vector<atom_pair> pairs;  // 1-4 pairs
    std::vector<atom_pair> others;
    std::vector<index_t> vdw_types;
    std::vector<lj_vdw> vdw_params;
    std::vector<param_t> partial_charges;
    std::vector<index_t> hydrophobic_atoms;
    std::vector<index_t> cation_atoms;
    std::vector<index_t> anion_atoms;
    std::vector<pi_ring> pi5_rings;
    std::vector<pi_ring> pi6_rings;
    std::vector<torsion_strain_penalty> bad_torsions;
    std::vector<atom_pair> hbonddon_charged;
    std::vector<atom_pair> hbonddon_neutral;
    std::vector<index_t> hbondacc_charged;
    std::vector<index_t> hbondacc_neutral;
    std::vector<atom_pair> rotatable_bonds;
    std::vector<atom_pair> all_bonds;  // `bond_index` is superset `bonds` in protein
    std::vector<index_t> atomic_numbers;
    param_t molecular_weight;
};

struct atom_position {
    param_t xyz[3];
};
typedef std::vector<atom_position> molecule_pose;

struct rigid_group {
    /**
     * Two atom indexs defining a torsion axis.
     * The first atom belongs to the position-fixed part of the receptor.
     * The second atom belongs this terminal group.
     */
    index_t axis[2];
    std::vector<index_t> atoms;
};

// receptor_variable=torsion
struct rotable_geometry {
    index_t num_atoms;
    index_t num_rotatable_bonds;
    std::vector<rigid_group> terminals;
    std::vector<index_t> bond_mappings;  // Bond order of Cpp engine to that in Python
};

// Atom indexes are those after position permutation
struct fragment {
    index_t id;  // Unique identifier
    std::vector<fragment> children;  // Fragment identifiers of children

    // Atom range
    index_t start;  // Included
    index_t end;  // Excluded

    // Rotatable bond
    index_t parent;
    index_t mine;
};

struct fragment_tree {
    index_t num_atoms;
    std::vector<fragment> roots;
    std::vector<fragment *> fragments;
    std::vector<index_t> permutation_index;  // Atoms in a fragement will be neighbors
    std::vector<index_t> reverse_index;
};

std::string to_string(const force_field_params& ffp, const bool ligand,
                      const std::string& prefix = "");
std::string to_string(const molecule_pose& ffp, const std::string& prefix = "");
std::string to_string(const rigid_group& terminal, const std::string& prefix = "");
std::string to_string(const rotable_geometry& geo, const std::string& prefix = "");
std::string to_string(const fragment& frag, const std::string& prefix = "");
std::string to_string(const fragment_tree& tree, const std::string& prefix = "");

}
