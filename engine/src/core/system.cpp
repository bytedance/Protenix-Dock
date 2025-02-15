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

#include "bytedock/core/system.h"

#include <queue>

#include "bytedock/lib/printer.h"

namespace bytedock {

static std::string to_string(const std::vector<harmonic_bond>& bonds,
                             const std::string& prefix) {
    std::ostringstream oss;
    oss << "{" << std::endl;
    oss << prefix + "  \"atom_ids\": [" << std::endl;
    const std::string indent(prefix + "    ");
    for (auto& bb : bonds) {
        oss << indent << "[" << bb.ids[0] << ", " << bb.ids[1] << "]," << std::endl;
    }
    oss << prefix + "  ]," << std::endl;
    oss << prefix + "  \"k\": [";
    for (auto& bb : bonds) oss << bb.k << ", ";
    oss << prefix + "]," << std::endl;
    oss << prefix + "  \"length\": [";
    for (auto& bb : bonds) oss << bb.length << ", ";
    oss << prefix + "]," << std::endl;
    oss << prefix << "}";
    return oss.str();
}

static std::string to_string(const std::vector<harmonic_angle>& angles,
                             const std::string& prefix) {
    std::ostringstream oss;
    oss << "{" << std::endl;
    oss << prefix + "  \"atom_ids\": [" << std::endl;
    const std::string indent(prefix + "    ");
    for (auto& aa : angles) {
        oss << indent << "[" << aa.ids[0] << ", " << aa.ids[1] << ", " << aa.ids[2]
            << "]," << std::endl;
    }
    oss << prefix + "  ]," << std::endl;
    oss << prefix + "  \"k\": [";
    for (auto& aa : angles) oss << aa.k << ", ";
    oss << prefix + "]," << std::endl;
    oss << prefix + "  \"degree\": [";
    for (auto& aa : angles) oss << aa.degree << ", ";
    oss << prefix + "]," << std::endl;
    oss << prefix << "}";
    return oss.str();
}

static std::string to_string(const std::vector<fourier_dihedral>& diherals,
                             const std::string& prefix) {
    std::ostringstream oss;
    oss << "{" << std::endl;
    oss << prefix + "  \"atom_ids\": [" << std::endl;
    const std::string indent(prefix + "    ");
    for (auto& dd : diherals) {
        oss << indent << "[" << dd.ids[0] << ", " << dd.ids[1] << ", "
            << dd.ids[2] << ", " << dd.ids[3] << "]," << std::endl;
    }
    oss << prefix + "  ]," << std::endl;
    oss << prefix + "  \"periodicity\": [" << std::endl;
    for (auto& dd : diherals) {
        oss << indent << "[";
        for (size_t i = 0; i < FF_DIHEDRAL_NDIMS; i++) {
            if (dd.periodicity[i] == 0) break;
            oss << dd.periodicity[i] << ", ";
        }
        oss << "]," << std::endl;
    }
    oss << prefix + "  ]," << std::endl;
    oss << prefix + "  \"k\": [" << std::endl;
    for (auto& dd : diherals) {
        oss << indent << "[";
        for (size_t i = 0; i < FF_DIHEDRAL_NDIMS; i++) {
            if (dd.periodicity[i] == 0) break;
            oss << dd.k[i] << ", ";
        }
        oss << "]," << std::endl;
    }
    oss << prefix + "  ]," << std::endl;
    oss << prefix + "  \"phase\": [" << std::endl;
    for (auto& dd : diherals) {
        oss << indent << "[";
        for (size_t i = 0; i < FF_DIHEDRAL_NDIMS; i++) {
            if (dd.periodicity[i] == 0) break;
            oss << dd.phase_raw[i] << ", ";
        }
        oss << "]," << std::endl;
    }
    oss << prefix + "  ]," << std::endl;
    oss << prefix << "}";
    return oss.str();
}

static std::string format_sigmas(const std::vector<lj_vdw>& values,
                                 const std::string& prefix) {
    std::ostringstream oss;
    size_t nrows = values.size() / NELEMENTS_PER_LINE;
    for (size_t i = 0; i < nrows; i++) {
        oss << prefix;
        for (size_t j = 0; j < NELEMENTS_PER_LINE; j++) {
            oss << values[i*NELEMENTS_PER_LINE+j].sigma << ", ";
        }
        oss << std::endl;
    }
    size_t remained = values.size() % NELEMENTS_PER_LINE;
    if (remained > 0) {
        nrows *= NELEMENTS_PER_LINE;  // Used as `offset`
        oss << prefix;
        for (size_t j = 0; j < remained; j++) oss << values[nrows+j].sigma << ", ";
        oss << std::endl;
    }
    return oss.str();
}

static std::string format_epsilons(const std::vector<lj_vdw>& values,
                                   const std::string& prefix) {
    std::ostringstream oss;
    size_t nrows = values.size() / NELEMENTS_PER_LINE;
    for (size_t i = 0; i < nrows; i++) {
        oss << prefix;
        for (size_t j = 0; j < NELEMENTS_PER_LINE; j++) {
            oss << values[i*NELEMENTS_PER_LINE+j].epsilon << ", ";
        }
        oss << std::endl;
    }
    size_t remained = values.size() % NELEMENTS_PER_LINE;
    if (remained > 0) {
        nrows *= NELEMENTS_PER_LINE;  // Used as `offset`
        oss << prefix;
        for (size_t j = 0; j < remained; j++) oss << values[nrows+j].epsilon << ", ";
        oss << std::endl;
    }
    return oss.str();
}

static std::string to_string(const std::vector<lj_vdw>& vdw_params,
                             const std::string& prefix) {
    std::ostringstream oss;
    const std::string indent(prefix + "    ");
    oss << "{" << std::endl;
    oss << prefix << "  \"sigma\": [" << std::endl;
    oss << format_sigmas(vdw_params, indent);
    oss << prefix << "  ]," << std::endl;
    oss << prefix << "  \"epsilon\": [" << std::endl;
    oss << format_epsilons(vdw_params, indent);
    oss << prefix << "  ]," << std::endl;
    oss << prefix << "}";
    return oss.str();
}

static std::string to_string(const std::vector<atom_pair>& pairs,
                             const std::string& prefix) {
    std::ostringstream oss;
    oss << "[" << std::endl;
    for (auto& item : pairs) {
        oss << prefix << "  [" << item.ids[0] << ", " << item.ids[1]
            << "]," << std::endl;
    }
    oss << prefix + "]";
    return oss.str();
}

static std::string to_string(const std::vector<pi_ring>& rings,
                             const std::string& prefix) {
    std::ostringstream oss;
    oss << "[" << std::endl;
    for (auto& item : rings) {
        oss << prefix << "  [" << item.ids[0];
        for (size_t i = 1; i < item.size; i++) oss << ", " << item.ids[i];
        oss << "]," << std::endl;
    }
    oss << prefix + "]";
    return oss.str();
}

static std::string to_string(const std::vector<torsion_index>& torsions,
                             const std::string& prefix) {
    std::ostringstream oss;
    const index_t* one;
    oss << "[" << std::endl;
    for (size_t i = 0; i < torsions.size(); i++) {
        one = torsions[i].ids;
        if (one[0] == one[1]) break;
        oss << prefix << "  [" << one[0] << ", " << one[1]
            << ", " << one[2] << ", " << one[3] << "]";
        if (i < torsions.size() - 1) oss << ",";
        oss << std::endl;
    }
    oss << prefix << "]";
    return oss.str();
}

static std::string to_string(const std::vector<torsion_strain_penalty>& penalties,
                             const std::string& prefix) {
    std::ostringstream oss;
    oss << "{" << std::endl;
    oss << prefix << "  \"coefficents\": [" << std::endl;
    for (auto& pp : penalties) {
        oss << "[" << pp.coefs[0] << ", " << pp.coefs[1]
            << ", " << pp.coefs[2] << ", " << pp.coefs[3] << "], ";
    }
    oss << prefix << "  ]," << std::endl;
    oss << prefix << "  \"atom_ids\": [" << std::endl;
    std::string indent(prefix + "    ");
    for (auto& pp : penalties) {
        oss << indent << to_string(pp.matched, indent) << ", " << std::endl;
    }
    oss << prefix << "  ]" << std::endl;
    oss << prefix << "}";
    return oss.str();
}

static std::string to_string(const std::vector<std::vector<index_t> >& groups,
                             const std::string& prefix) {
    std::ostringstream oss;
    oss << "[" << std::endl;
    for (auto& gg : groups) {
        oss << prefix << "  [";
        for (auto& aa : gg) oss << aa << ", ";
        oss << "]," << std::endl;
    }
    oss << prefix + "]";
    return oss.str();
}

std::string to_string(const force_field_params& ffp, const bool ligand, const std::string& prefix) {
    std::ostringstream oss;
    const std::string indent(prefix + "  ");
    oss << "{" << std::endl;
    oss << prefix << "  \"harmonic_bonds\": " << to_string(ffp.bonds, indent)
        << "," << std::endl;
    oss << prefix << "  \"harmonic_angles\": " << to_string(ffp.angles, indent)
        << "," << std::endl;
    oss << prefix << "  \"improper_torsions\": " << to_string(ffp.impropers, indent)
        << "," << std::endl;
    oss << prefix << "  \"proper_torsions\": " << to_string(ffp.propers, indent)
        << "," << std::endl;
    oss << prefix << "  \"pairs_14\": " << to_string(ffp.pairs, indent)
        << "," << std::endl;
    oss << prefix << "  \"pairs_intra\": " << to_string(ffp.others, indent)
        << "," << std::endl;
    if (ligand) {
        oss << prefix << "  \"vdw_types\": [" << std::endl
            << format(ffp.vdw_types, indent)
            << "]," << std::endl;
    } else {
        oss << prefix << "  \"vdw_params\": " << to_string(ffp.vdw_params, indent)
            << "," << std::endl;
    }
    oss << prefix << "  \"partial_charges\": [" << std::endl
        << format(ffp.partial_charges, indent)
        << prefix << "  ]," << std::endl;
    oss << prefix << "  \"hydrophobic_atoms\": [" << std::endl
        << format(ffp.hydrophobic_atoms, indent)
        << prefix << "  ]," << std::endl;
    oss << prefix << "  \"cation_atoms\": [" << std::endl
        << format(ffp.cation_atoms, indent)
        << prefix << "  ]," << std::endl;
    oss << prefix << "  \"anion_atoms\": [" << std::endl
        << format(ffp.anion_atoms, indent)
        << prefix << "  ]," << std::endl;
    oss << prefix << "  \"pi5_rings\": " << to_string(ffp.pi5_rings, indent)
        << "," << std::endl;
    oss << prefix << "  \"pi6_rings\": " << to_string(ffp.pi6_rings, indent)
        << "," << std::endl;
    if (ffp.bad_torsions.size() > 0) {  // Only ligand has penalties
        oss << prefix << "  \"bad_torsions\": " << to_string(ffp.bad_torsions, indent)
            << "," << std::endl;
    }
    oss << prefix << "  \"hbonddon_charged\": "
        << to_string(ffp.hbonddon_charged, indent)
        << "," << std::endl;
    oss << prefix << "  \"hbonddon_neutral\": "
        << to_string(ffp.hbonddon_neutral, indent)
        << "," << std::endl;
    oss << prefix << "  \"hbondacc_charged\": [" << std::endl
        << format(ffp.hbondacc_charged, indent)
        << prefix << "  ]," << std::endl;
    oss << prefix << "  \"hbondacc_neutral\": [" << std::endl
        << format(ffp.hbondacc_neutral, indent)
        << prefix << "  ]," << std::endl;
    if (ffp.rotatable_bonds.size() > 0) {
        oss << prefix << "  \"rotatable_bonds\": "
            << to_string(ffp.rotatable_bonds, indent)
            << "," << std::endl;
    }
    oss << prefix << "  \"atomic_numbers\": [" << std::endl
        << format(ffp.atomic_numbers, indent)
        << prefix << "  ]," << std::endl;
    if (ffp.molecular_weight > 0_r) {
        oss << prefix << "  \"molecular_weight\": " << ffp.molecular_weight
            << "," << std::endl;
    }
    oss << prefix << "  \"bond_index\": "
        << to_string(ffp.all_bonds, indent)
        << std::endl;
    oss << prefix << "}";
    return oss.str();
}

std::string to_string(const molecule_pose& pose, const std::string& prefix) {
    if (pose.size() == 0) return "[]";
    std::ostringstream oss;
    oss << "[" << std::endl;
    size_t last = pose.size() - 1;
    for (size_t i = 0; i < last; i++) {
        auto& item = pose[i];
        oss << prefix << "  ["
            << item.xyz[0] << ", " << item.xyz[1] << ", " << item.xyz[2]
            << "]," << std::endl;
    }
    auto& item = pose[last];
    oss << prefix << "  ["
        << item.xyz[0] << ", " << item.xyz[1] << ", " << item.xyz[2]
        << "]" << std::endl;
    oss << prefix << "]";
    return oss.str();
}

std::string to_string(const rigid_group& terminal, const std::string& prefix) {
    std::ostringstream oss;
    oss << "{" << std::endl;
    oss << prefix << "  \"axis\": ["
        << terminal.axis[0] << ", " << terminal.axis[1]
        << "]," << std::endl;
    oss << prefix << "  \"atoms\": [";
    for (auto& num: terminal.atoms) oss << num << ", ";
    oss << "]," << std::endl;
    oss << prefix << "}";
    return oss.str();
}

static std::string to_string(const std::vector<rigid_group> terminals,
                             const std::string& prefix = "") {
    std::ostringstream oss;
    oss << "[" << std::endl;
    const std::string indent(prefix + "  ");
    for (auto& one : terminals) {
        oss << indent << to_string(one, indent) << "," << std::endl;
    }
    oss << prefix << "]";
    return oss.str();
}

std::string to_string(const rotable_geometry& geo, const std::string& prefix) {
    std::ostringstream oss;
    const std::string indent(prefix + "  ");
    oss << "{" << std::endl;
    oss << prefix << "  \"num_atoms\": " << geo.num_atoms << "," << std::endl;
    oss << prefix << "  \"num_rotatable_bonds\": " << geo.num_rotatable_bonds
        << "," << std::endl;
    oss << prefix << "  \"terminals\": " << to_string(geo.terminals, indent)
        << "," << std::endl;
    oss << prefix << "  \"rot_bond_orders\": [" << std::endl
        << format<index_t, 10>(geo.bond_mappings, indent)
        << prefix << "  ]" << std::endl;
    oss << prefix << "}";
    return oss.str();
}

std::string to_string(const fragment& frag, const std::string& prefix) {
    std::ostringstream oss;
    oss << "{" << std::endl;
    std::string indent(prefix + "  ");
    oss << indent << "\"id\": " << frag.id << "," << std::endl;
    oss << indent << "\"index_range\": ["
        << frag.start << ", " << frag.end << "]," << std::endl;
    oss << indent << "\"rotatable_bond\": ["
        << frag.parent << ", " << frag.mine << "]," << std::endl;
    if (frag.children.size() > 0) {
        oss << indent << "\"children\": [";
        size_t i;
        for (i = 0; i < frag.children.size() - 1; ++i) {
            oss << frag.children[i].id << ", ";
        }
        oss << frag.children[i].id << "]" << std::endl;
    } else {
        oss << indent << "\"children\": []" << std::endl;
    }
    oss << prefix << "}";
    return oss.str();
}

// Traverse the fragment tree in a layer-first way
std::string to_string(const std::vector<fragment>& frag_list,
                      const std::string& prefix) {
    std::ostringstream oss;
    oss << "[" << std::endl;
    std::string indent(prefix + "  ");
    std::queue<const fragment *> q;
    for (auto& item: frag_list) q.push(&item);
    const fragment* current;
    while (q.size() > 0) {
        if (oss.tellp() > 2) oss << "," << std::endl;
        current = q.front();
        q.pop();
        oss << indent << to_string(*current, indent);
        for (auto& child : current->children) q.push(&child);
    }
    oss << std::endl;
    oss << prefix << "]";
    return oss.str();
}

std::string to_string(const fragment_tree& tree, const std::string& prefix) {
    std::ostringstream oss;
    oss << "{" << std::endl;
    std::string indent(prefix + "  ");
    oss << indent << "\"num_atoms\": " << tree.num_atoms << "," << std::endl;
    oss << indent << "\"permutation_index\": [" << std::endl
        << format<index_t, 10>(tree.permutation_index, indent)
        << indent << "]," << std::endl;
    oss << indent << "\"reverse_index\": [" << std::endl
        << format<index_t, 10>(tree.reverse_index, indent)
        << indent << "]," << std::endl;
    oss << indent << "\"roots\": " << to_string(tree.roots, indent) << std::endl;
    oss << prefix << "}";
    return oss.str();
}

}
