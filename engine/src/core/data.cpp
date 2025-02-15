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

#include "bytedock/core/data.h"

#include <sstream>
#include <unordered_map>

#include <boost/lexical_cast.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "bytedock/core/constant.h"
#include "bytedock/ext/logging.h"
#include "bytedock/lib/error.h"

namespace bytedock {

namespace pt = boost::property_tree;

static void fill(const pt::ptree& node, std::vector<harmonic_bond>& dest) {
    size_t dim = node.get_child("FF_Bonds_atomidx").size();
    dest.resize(dim);
    size_t idx = 0;
    for (auto& kv : node.get_child("FF_Bonds_atomidx")) {
        dim = 0;
        for (auto& ij : kv.second) {
            dest[idx].ids[dim] = ij.second.get_value<index_t>();
            dim++;
        }
        idx++;
    }
    idx = 0;
    for (auto& kv : node.get_child("FF_Bonds_k")) {
        dest[idx].k = kv.second.get_value<param_t>();
        idx++;
    }
    idx = 0;
    for (auto& kv : node.get_child("FF_Bonds_length")) {
        dest[idx].length = kv.second.get_value<param_t>();
        idx++;
    }
}

static void fill(const pt::ptree& node, std::vector<harmonic_angle>& dest) {
    size_t dim = node.get_child("FF_Angles_atomidx").size();
    dest.resize(dim);
    size_t idx = 0;
    for (auto& kv : node.get_child("FF_Angles_atomidx")) {
        dim = 0;
        for (auto& ijk : kv.second) {
            dest[idx].ids[dim] = ijk.second.get_value<index_t>();
            dim++;
        }
        idx++;
    }
    idx = 0;
    for (auto& kv : node.get_child("FF_Angles_k")) {
        dest[idx].k = kv.second.get_value<param_t>();
        idx++;
    }
    idx = 0;
    for (auto& kv : node.get_child("FF_Angles_angle")) {
        dest[idx].degree = kv.second.get_value<param_t>();
        dest[idx].radian = MATH_TO_RADIAN(dest[idx].degree);
        idx++;
    }
}

static void fill(const pt::ptree& node, const std::string& prefix,
                 std::vector<fourier_dihedral>& dest) {
    size_t dim = node.get_child(prefix + "_atomidx").size();
    dest.resize(dim);
    size_t idx = 0;
    for (auto& kv : node.get_child(prefix + "_atomidx")) {
        dim = 0;
        for (auto& ijkl : kv.second) {
            dest[idx].ids[dim] = ijkl.second.get_value<index_t>();
            dim++;
        }
        idx++;
    }
    idx = 0;
    for (auto& kv : node.get_child(prefix + "_periodicity")) {
        dim = 0;
        for (auto& elem : kv.second) {
            dest[idx].periodicity[dim] = elem.second.get_value<index_t>();
            dim++;
        }
        while (dim < FF_DIHEDRAL_NDIMS) {
            dest[idx].periodicity[dim] = 0;
            dim++;
        }
        idx++;
    }
    idx = 0;
    for (auto& kv : node.get_child(prefix + "_k")) {
        dim = 0;
        for (auto& elem : kv.second) {
            dest[idx].k[dim] = elem.second.get_value<param_t>();
            dim++;
        }
        idx++;
    }
    idx = 0;
    for (auto& kv : node.get_child(prefix + "_phase")) {
        dim = 0;
        for (auto& elem : kv.second) {
            dest[idx].phase_raw[dim] = elem.second.get_value<param_t>();
            dest[idx].phase[dim] = MATH_TO_RADIAN(dest[idx].phase_raw[dim]);
            dim++;
        }
        idx++;
    }
}

static void fill(const pt::ptree& node, std::vector<lj_vdw>& dest) {
    size_t idx = node.get_child("FF_vdW_sigma").size();
    dest.resize(idx);
    idx = 0;
    for (auto& kv : node.get_child("FF_vdW_sigma")) {
        dest[idx].sigma = kv.second.get_value<param_t>();
        idx++;
    }
    idx = 0;
    for (auto& kv : node.get_child("FF_vdW_epsilon")) {
        dest[idx].epsilon = kv.second.get_value<param_t>();
        idx++;
    }
}

inline void fill(const pt::ptree& node, std::vector<index_t>& dest) {
    dest.resize(node.size());
    size_t idx = 0;
    for (auto& one : node) {
        dest[idx] = one.second.get_value<index_t>();
        idx++;
    }
}

inline void fill(const pt::ptree& node, std::vector<param_t>& dest) {
    dest.resize(node.size());
    size_t idx = 0;
    for (auto& one : node) {
        dest[idx] = one.second.get_value<param_t>();
        idx++;
    }
}

static void fill(const pt::ptree& node, std::vector<atom_pair>& dest) {
    dest.resize(node.size());
    size_t i = 0, j;
    for (auto& pair : node) {
        j = 0;
        for (auto& atom : pair.second) {
            dest[i].ids[j] = atom.second.get_value<index_t>();
            j++;
        }
        i++;
    }
}

static void fill(const pt::ptree& node, std::vector<pi_ring>& dest) {
    dest.resize(node.size());
    size_t i = 0, j;
    for (auto& pair : node) {
        j = 0;
        for (auto& atom : pair.second) {
            dest[i].ids[j] = atom.second.get_value<index_t>();
            j++;
        }
        dest[i].size = j;
        i++;
    }
}

static void fill(const pt::ptree& node, std::vector<torsion_index>& dest) {
    size_t idx = node.size(), dim;
    dest.resize(idx);
    idx = 0;
    for (auto& outer : node) {
        dim = 0;
        for (auto& inner : outer.second) {
            dest[idx].ids[dim] = inner.second.get_value<index_t>();
            dim++;
        }
        idx++;
    }
}

static void fill(const pt::ptree& node, std::vector<torsion_strain_penalty>& dest) {
    size_t idx = node.get_child("tstrain_params_pose_selection").size();
    if (idx == 0) return;
    dest.resize(idx);
    size_t dim;
    idx = 0;
    for (auto& outer : node.get_child("tstrain_params_pose_selection")) {
        dim = 0;
        for (auto& inner : outer.second) {
            dest[idx].coefs[dim] = inner.second.get_value<param_t>();
            dim++;
        }
        idx++;
    }
    idx = 0;
    for (auto& outer : node.get_child("tstrain_atomidx_pose_selection")) {
        fill(outer.second, dest[idx].matched);
        idx++;
    }
}

inline void fill(const pt::ptree& node, std::vector<std::vector<index_t> >& dest) {
    dest.resize(node.size());
    size_t i = 0;
    for (auto& parent : node) {
        fill(parent.second, dest[i]);
        i++;
    }
}

inline void nolmaliz_vdw_types(std::vector<index_t>& vdw_types) {
    for (size_t i = 0; i < vdw_types.size(); i++) {
        vdw_types[i] -= KVdwOffset;
        if (vdw_types[i] > kNumVdwTypes || vdw_types[i] < 0) {
            std::ostringstream oss;
            oss << "VDW type [" << vdw_types[i] << "] of atom#"
                << i << "is illegal!";
            throw failed_conf_error(oss.str());
        }
    }
}

inline void fill(const pt::ptree& node, const bool ligand, force_field_params& dest) {
    fill(node, dest.bonds);
    fill(node, dest.angles);
    fill(node, "FF_ImproperTorsions", dest.impropers);
    fill(node, "FF_ProperTorsions", dest.propers);
    if (ligand) {
        fill(node.get_child("FF_vdW_paraidx"), dest.vdw_types);
        nolmaliz_vdw_types(dest.vdw_types);
    } else {
        fill(node, dest.vdw_params);
    }
    fill(node.get_child("partial_charges"), dest.partial_charges);
    fill(node.get_child("FF_Nonbonded14_atomidx"), dest.pairs);
    fill(node.get_child("FF_NonbondedAll_atomidx"), dest.others);
    fill(node.get_child("hydrophobic_atomidx"), dest.hydrophobic_atoms);
    fill(node.get_child("cation_atomidx"), dest.cation_atoms);
    fill(node.get_child("anion_atomidx"), dest.anion_atoms);
    fill(node.get_child("piring5_atomidx"), dest.pi5_rings);
    fill(node.get_child("piring6_atomidx"), dest.pi6_rings);
    fill(node.get_child("hbonddon_charged_atomidx"), dest.hbonddon_charged);
    fill(node.get_child("hbonddon_neut_atomidx"), dest.hbonddon_neutral);
    fill(node.get_child("hbondacc_charged_atomidx"), dest.hbondacc_charged);
    fill(node.get_child("hbondacc_neut_atomidx"), dest.hbondacc_neutral);
    if (ligand) {
        fill(node, dest.bad_torsions);
        fill(node.get_child("rotatable_bond_index"), dest.rotatable_bonds);
    }
    dest.molecular_weight = node.get<param_t>("molecular_weight", -1.);
    fill(node.get_child("atomic_numbers"), dest.atomic_numbers);
    fill(node.get_child("bond_index"), dest.all_bonds);
}

static void fill(const pt::ptree& node, molecule_pose& dest) {
    const size_t natoms = node.size();
    dest.resize(natoms);
    size_t idx = 0, dim;
    for (auto& ii : node) {
        dim = 0;
        for (auto& jj : ii.second) {
            dest[idx].xyz[dim] = jj.second.get_value<param_t>();
            dim++;
        }
        idx++;
    }
}

static void fill(const pt::ptree& node, rotable_geometry& dest) {
    dest.num_atoms = node.get<index_t>("num_atoms");
    dest.terminals.clear();
    size_t offset, idx;
    for (auto& kv : node.get_child("data")) {
        auto& batch = kv.second;
        offset = dest.terminals.size();
        dest.terminals.resize(offset + batch.get_child("rot_bond_1").size());
        idx = 0;
        for (auto& one : batch.get_child("rot_bond_1")) {
            dest.terminals[offset + idx].axis[0] = one.second.get_value<index_t>();
            idx++;
        }
        idx = 0;
        for (auto& one : batch.get_child("rot_bond_2")) {
            dest.terminals[offset + idx].axis[1] = one.second.get_value<index_t>();
            idx++;
        }
        idx = 0;
        for (auto& ii : batch.get_child("fragments")) {
            auto& group = dest.terminals[offset + idx].atoms;
            for (auto& jj : ii.second) group.push_back(jj.second.get_value<index_t>());
            idx++;
        }
    }
    dest.bond_mappings.resize(dest.terminals.size());
    dest.num_rotatable_bonds = node.get<index_t>("num_rotatable_bonds", 0);
    idx = 0;
    for (auto& kv : node.get_child("data")) {
        auto& batch = kv.second;
        for (auto& one : batch.get_child("rot_bond_order")) {
            dest.bond_mappings[idx] = one.second.get_value<index_t>();
            offset = dest.bond_mappings[idx] + 1;  // Used as `tmp`
            if (offset > dest.num_rotatable_bonds) dest.num_rotatable_bonds = offset;
            ++idx;
        }
    }
}

void torsional_receptor::parse(std::shared_ptr<std::istream> fin) {
    step_timer t(PARSE_RECEPTOR_DATA_STEP);
    pt::ptree properties;
    pt::read_json(*fin, properties);
    fill(properties.get_child("geometry"), geo_);
    LOG_DEBUG << "\"torsional_receptor.geometry\": " << to_string(geo_);
    fill(properties.get_child("ffdata"), false, ffp_);
    LOG_DEBUG << "\"torsional_receptor.ff_data\": " << to_string(ffp_, false);
    fill(properties.get_child("xyz"), coords_);
    LOG_DEBUG << "\"torsional_receptor.xyz\": " << to_string(coords_);
    if (geo_.num_atoms != coords_.size()) {
        std::ostringstream oss;
        oss << "Atom count mismatches! Wanted [" << geo_.num_atoms << "] vs. Parsed ["
            << coords_.size() << "]";
        throw failed_conf_error(oss.str());
    }

    // Precalculate rotation-related fields
    index_t axis_atom;
    leaves_.resize(geo_.terminals.size());
    nmovable_ = 0;
    for (size_t i = 0; i < leaves_.size(); i++) {
        auto& leaf = leaves_[i];
        const auto& rgroup = geo_.terminals[i];
        leaf.axis_origin = coords_[rgroup.axis[0]];
        axis_atom = rgroup.axis[1];
        minus_3d(coords_[axis_atom].xyz, leaf.axis_origin.xyz, leaf.axis_vector.xyz);
        normalize_3d(leaf.axis_vector.xyz);
        for (const auto& j : rgroup.atoms) {
            if (j != axis_atom) leaf.movable_atoms.push_back(j);
        }
        nmovable_ += leaf.movable_atoms.size();
    }
}

molecule_pose torsional_receptor::apply_parameters(
    const std::vector<param_t>& torsions
) const {
    molecule_pose copied = coords_;
    matrix_3x3 rot;
    vector_3d tmp;
    for(size_t i = 0; i < torsions.size(); i++) {
        auto& leaf = leaves_[i];
        get_rotation_matrix(leaf.axis_vector, torsions[i], rot);

        // Apply rotation matrix to all atoms of the terminal
        for (const auto& j : leaf.movable_atoms) {
            atom_position& dest_pos = copied[j];
            minus_3d(dest_pos.xyz, leaf.axis_origin.xyz, tmp.xyz);
            multiply_3x3_3d(rot.data, tmp.xyz, dest_pos.xyz);
            add_inplace_3d(dest_pos.xyz, leaf.axis_origin.xyz);
        }
    }
    return copied;
}

inline void fill(const pt::ptree& node, std::vector<molecule_pose>& dest) {
    dest.resize(node.size());
    size_t idx = 0;
    for (auto& ii : node) {
        fill(ii.second, dest[idx]);
        idx++;
    }
}

static void pre_order_traverse(
    fragment& node, std::unordered_map<std::string, std::pair<index_t, index_t> >& repo
) {
    std::ostringstream oss;
    for (auto& child : node.children) {
        oss.str("");
        oss << node.id << "->" << child.id;
        auto& bond = repo[oss.str()];
        child.parent = bond.first;
        child.mine = bond.second;
        pre_order_traverse(child, repo);
    }
}

static void fill(const pt::ptree& node, fragment_tree& dest) {
    dest.num_atoms = node.get<index_t>("num_atoms");

    // Parse fragment tree 
    auto& frag_split_index = node.get_child("frag_split_index");
    dest.fragments.resize(frag_split_index.size() + 1);
    index_t frag_id, child_id;
    fragment* out_frag;
    for (auto& in_level : node.get_child("frag_traverse_levels")) {
        for (auto& in_frag : in_level.second) {
            frag_id = boost::lexical_cast<index_t>(in_frag.first);
            auto& parents = in_frag.second.get_child("parent");
            auto& in_children = in_frag.second.get_child("children");
            if (parents.size() == 0) {  // Since no parent, introduce myself
                child_id = dest.roots.size();
                dest.roots.push_back({});
                dest.roots[child_id].id = frag_id;
                dest.fragments[frag_id] = dest.roots.data() + child_id;
            }
            out_frag = dest.fragments[frag_id];
            out_frag->children.resize(in_children.size());
            child_id = 0;
            for (auto& in_child : in_children) {  // Introduce my children
                frag_id = in_child.second.get_value<index_t>();
                out_frag->children[child_id].id = frag_id;
                dest.fragments[frag_id] = out_frag->children.data() + child_id;
                ++child_id;
            }
        }
    }

    // Fill in atom range
    std::vector<index_t> atom_range(dest.fragments.size() + 1);
    atom_range[0] = 0;
    atom_range[dest.fragments.size()] = dest.num_atoms;
    frag_id = 0;
    child_id = 0;  // Used as `prev_id`
    for (auto& one : frag_split_index) {
        out_frag = dest.fragments[frag_id];
        out_frag->start = child_id;
        child_id = one.second.get_value<index_t>();
        out_frag->end = child_id;
        ++frag_id;
    }
    out_frag = dest.fragments[frag_split_index.size()];
    out_frag->start = child_id;
    out_frag->end = dest.num_atoms;

    // Parse permutation index & reverse one
    auto& permute_index = node.get_child("permute_index");
    dest.permutation_index.resize(permute_index.size());
    dest.reverse_index.resize(permute_index.size());
    frag_id = 0;  // Used as `offset`
    for (auto& one : permute_index) {
        child_id = one.second.get_value<index_t>();  // Used as `value`
        dest.permutation_index[frag_id] = child_id;
        dest.reverse_index[child_id] = frag_id;
        frag_id++;
    }

    // Fill in bonds
    std::unordered_map<std::string, std::pair<index_t, index_t> > bond_map;
    std::vector<index_t> bond(2);
    for (auto& outer : node.get_child("edge_to_bond")) {
        child_id = 0;  // Used as offset
        for (auto& inner : outer.second) {
            bond[child_id] = dest.reverse_index[inner.second.get_value<index_t>()];
            ++child_id;
        }
        bond_map[outer.first] = std::make_pair(bond[0], bond[1]);
    }
    for (auto& child : dest.roots) pre_order_traverse(child, bond_map);
}

static std::string to_string(const std::vector<molecule_pose>& poses,
                             const std::string& prefix = "") {
    std::ostringstream oss;
    oss << "[" << std::endl;
    const std::string indent(prefix + "  ");
    for (auto& item : poses) {
        oss << indent << to_string(item, indent) << "," << std::endl;
    }
    oss << prefix << "]" << std::endl;
    return oss.str();
}

void free_ligand::parse(std::shared_ptr<std::istream> fin, const bool geometry) {
    step_timer t(PARSE_LIGAND_DATA_STEP);
    pt::ptree properties;
    pt::read_json(*fin, properties);
    if (geometry) {
        fill(properties.get_child("geometry"), geo_);
        LOG_DEBUG << "\"free_ligand.geometry\": " << to_string(geo_);
    }
    if (properties.count("mapped_smiles") > 0) {
        smiles_ = properties.get<std::string>("mapped_smiles");
    }
    fill(properties.get_child("ffdata"), true, ffp_);
    LOG_DEBUG << "\"free_ligand.ff_data\": " << to_string(ffp_, true);
    fill(properties.get_child("xyz"), poses_);
    LOG_DEBUG << "\"free_ligand.nposes\": " << poses_.size();
    LOG_DEBUG << "\"free_ligand.xyz\": " << to_string(poses_);
    for (size_t i = 0; i < ffp_.vdw_types.size(); i++) {
        if (ffp_.vdw_types[i] > kNumVdwTypes) {
            std::ostringstream oss;
            oss << "VDW type [" << ffp_.vdw_types[i] << "] of atom#"
                << i << "is illegal!";
            throw failed_conf_error(oss.str());
        }
    }
}

inline molecule_pose rearrange_atoms(const molecule_pose& input,
                                     const std::vector<index_t>& permutation) {
    molecule_pose output(input.size());
    for (size_t i = 0; i < permutation.size(); ++i) output[i] = input[permutation[i]];
    return output;
}

// Parse rotatable bonds into axes
static void pre_order_traverse(const fragment& node, const molecule_pose& ligand_xyz,
                               std::vector<vector_3d>& rb_axes) {
    minus_3d(ligand_xyz[node.mine].xyz, ligand_xyz[node.parent].xyz,
             rb_axes[node.id].xyz);
    normalize_3d(rb_axes[node.id].xyz);
    for (auto& child : node.children) pre_order_traverse(child, ligand_xyz, rb_axes);
}

// Update atom coordates: world -> local
static void post_order_traverse(const fragment& node, molecule_pose& ligand_xyz) {
    for (auto& child : node.children) {
        post_order_traverse(child, ligand_xyz);
        minus_inplace_3d(ligand_xyz[child.mine].xyz, ligand_xyz[node.mine].xyz);
    }
    size_t i = node.start;
    while (i < node.mine) {
        minus_inplace_3d(ligand_xyz[i].xyz, ligand_xyz[node.mine].xyz);
        ++i;
    }
    ++i;  // Handled by function caller / in parent fragment
    while (i < node.end) {
        minus_inplace_3d(ligand_xyz[i].xyz, ligand_xyz[node.mine].xyz);
        ++i;
    }
}

// Update atom coordinates: local -> torsion -> world
static void pre_order_traverse(const fragment& node,
                               const matrix_3x3& parent_mat,
                               std::vector<param_t>::const_iterator& iter,
                               std::vector<vector_3d>& rb_axes,
                               molecule_pose& ligand_xyz) {
    param_t theta = *iter;
    ++iter;
    matrix_3x3 curr_mat;
    get_rotation_matrix(rb_axes[node.id], theta, curr_mat);
    auto fused_mat = multiply_matrix_3x3(curr_mat, parent_mat);

    // `node.mine` is rotated in function caller
    vector_3d tmp;
    param_t *origin = ligand_xyz[node.mine].xyz;
    for (size_t i = node.start; i < node.end; ++i) {
        if (i != node.mine) {
            multiply_3x3_3d(fused_mat.data, ligand_xyz[i].xyz, tmp.xyz);
            add_3d(tmp.xyz, origin, ligand_xyz[i].xyz);
        }
    }

    // Rotate children's axes and bond atoms
    for (auto& child : node.children) {
        multiply_3x3_3d(fused_mat.data, rb_axes[child.id].xyz, tmp.xyz);
        rb_axes[child.id] = tmp;
        multiply_3x3_3d(fused_mat.data, ligand_xyz[child.mine].xyz, tmp.xyz);
        add_3d(tmp.xyz, origin, ligand_xyz[child.mine].xyz);
    }

    // Handle children fragments
    for (auto& child : node.children) {
        pre_order_traverse(child, fused_mat, iter, rb_axes, ligand_xyz);
    }
}

molecule_pose free_ligand::apply_parameters(
    const molecule_pose& ligand_xyz, const atom_position& new_center,
    const matrix_3x3& orientation, const std::vector<param_t>& torsions
) const {
    auto copied = rearrange_atoms(ligand_xyz, geo_.permutation_index);  // Move

    // Convert rotatable bonds to axes
    std::vector<vector_3d> rb_axes(geo_.fragments.size());
    for (auto& outer : geo_.roots) {
        // Root fragment have no rotatable bonds but orientation
        for (auto& inner : outer.children) pre_order_traverse(inner, copied, rb_axes);
    }

    // Convert atom coordinates to relative ones for rotation & torsions
    atom_position old_center = calculate_geometric_center(copied);
    for (auto& outer : geo_.roots) {
        for (auto& inner : outer.children) {
            post_order_traverse(inner, copied);
            minus_inplace_3d(copied[inner.mine].xyz, old_center.xyz);
        }
        for (size_t i = outer.start; i < outer.end; ++i) {
            minus_inplace_3d(copied[i].xyz, old_center.xyz);
        }
    }

    // Apply transition & orientation to root node
    vector_3d tmp;
    for (auto& outer : geo_.roots) {
        for (size_t i = outer.start; i < outer.end; ++i) {
            multiply_3x3_3d(orientation.data, copied[i].xyz, tmp.xyz);
            add_3d(tmp.xyz, new_center.xyz, copied[i].xyz);
        }
        // Rotate children's axes and bond atoms
        for (auto& inner : outer.children) {
            multiply_3x3_3d(orientation.data, rb_axes[inner.id].xyz, tmp.xyz);
            rb_axes[inner.id] = tmp;
            multiply_3x3_3d(orientation.data, copied[inner.mine].xyz, tmp.xyz);
            add_3d(tmp.xyz, new_center.xyz, copied[inner.mine].xyz);
        }
    }

    // Handle children fragments
    std::vector<param_t>::const_iterator iter = torsions.begin();
    for (auto& outer : geo_.roots) {
        for (auto& inner : outer.children) {
            pre_order_traverse(inner, orientation, iter, rb_axes, copied);
        }
    }
    return rearrange_atoms(copied, geo_.reverse_index);
}

}
