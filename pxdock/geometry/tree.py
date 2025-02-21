# Copyright (C) 2025 ByteDance and/or its affiliates

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import copy
import warnings
from typing import Iterable

import networkx as nx
import numpy as np


def create_graph(bond_index: list[tuple[int, int]], num_atoms: int) -> nx.Graph:
    G = nx.Graph()
    G.add_nodes_from(range(num_atoms))
    G.add_edges_from(bond_index)
    return G


def find_tree_longest_path(tree: nx.Graph) -> list[int]:
    farthest_node_1 = list(nx.bfs_tree(tree, source=list(tree.nodes())[0]))[-1]
    farthest_node_2 = list(nx.bfs_tree(tree, source=farthest_node_1))[-1]
    longest_path = nx.shortest_path(
        tree, source=farthest_node_1, target=farthest_node_2
    )

    return longest_path


def find_tree_center(tree):
    """
    find the longest path of the tree, and return the center node of the longest path.
    """
    longest_path = find_tree_longest_path(tree)
    # if longest_path is even, there exist two center nodes, return one of them.
    center = longest_path[len(longest_path) // 2]
    return center


def next_level_traversal(tree, current_level, visited: set) -> set:
    visited = set(visited)
    current_level = set(current_level)
    next_level = set()
    for node in current_level:
        neighbors = set(int(n) for n in tree.neighbors(node))
        children = neighbors - visited
        n = len(next_level)
        next_level.update(children)
        assert len(next_level) == n + len(
            children
        ), "Something must be wrong. The input may not be a tree."
    return next_level


def traverse_next_level(tree: nx.graph, current_level: list, visited: set):
    visited = set(visited)
    current_level = set(current_level)
    next_level = {}
    for node in current_level:
        neighbors = set(int(n) for n in tree.neighbors(node))
        children = neighbors - visited
        parent = neighbors - children
        assert len(parent) <= 1, "Something must be wrong. The input may not be a tree."
        next_level[node] = {"children": sorted(list(children)), "parent": list(parent)}
    return next_level


def level_order_traversal(tree, start_level):
    levels = []
    current_level = list(set(start_level))
    visited = set([])

    while current_level:
        visited.update(set(current_level))
        next_level = traverse_next_level(tree, current_level, visited)
        if next_level:
            levels.append(next_level)
        current_level = sum([val["children"] for val in next_level.values()], [])
    return levels


class RigidFragmentTree:
    def __init__(
        self,
        bond_index: list[list[int]],
        is_rotatable: list[int],
        num_atoms: int = None,
    ):
        """
        Representation of a molecule using a tree constructed by its rigid fragments and rotatable bonds connecting each other

        Args:
            bond_index: List[List[int]] with size [nbond, 2]. Each row contains atom indices of a bond.
            is_rotatable: List[int] with ize [nbond]. Each element indicates whether the bond is rotatable.
            num_atoms: int. The number of atoms in the molecule.
        """
        is_rotatable = np.array(is_rotatable).astype(dtype=bool).reshape(-1)
        bond_index = np.array(bond_index).astype(int)
        num_bonds = len(bond_index)
        num_atoms = num_atoms or bond_index.max() + 1

        assert num_atoms > bond_index.max()
        assert num_bonds == len(is_rotatable)

        self.bond_index, self.is_rotatable = bond_index, is_rotatable

        rbond_index = bond_index[is_rotatable]

        # molecule graph
        self.G_mol = create_graph(bond_index, num_atoms)

        # fragment quotient graph
        self.G, self.edge_to_bond = self.quotient_graph_of_rigid_fragments(
            self.G_mol, edge_cut=rbond_index
        )
        self.ordered_edges = [e for e in self.G.edges]
        self.edge_to_order = {
            f"{e[0]}->{e[1]}": i for i, e in enumerate(self.ordered_edges)
        }
        self.edge_to_order.update(
            {f"{e[1]}->{e[0]}": i for i, e in enumerate(self.ordered_edges)}
        )

        # basic info
        self.number_of_nodes = self.G.number_of_nodes()
        self.number_of_edges = len(self.ordered_edges)
        self.number_of_atoms = num_atoms

        self._root = None

    def find_tree_center(self, tree: nx.Graph):
        if tree.number_of_nodes() <= 2:
            # choose the larger component
            max_size = -1
            for i in tree.nodes():
                n_size = len(tree.nodes[i]["atoms"])
                if n_size > max_size:
                    root = i
                    max_size = n_size
            return root
        else:
            # choose the node in the center of the longest path
            return find_tree_center(tree)

    @property
    def root(self):
        """Find the tree center of G_frag. Use it as the root frame."""
        if self._root is None:
            if not nx.is_tree(self.G):
                components = list(nx.connected_components(self.G))
                trees = [self.G.subgraph(cp) for cp in components]
                self._root = [self.find_tree_center(tree) for tree in trees]
            else:
                self._root = self.find_tree_center(self.G)
        return self._root

    def levels(self, start_level=None):
        """Get the level-order traversal of the tree"""
        if start_level is None:
            if isinstance(self.root, int):
                start_level = [self.root]
            else:
                start_level = self.root
        return level_order_traversal(self.G, start_level=start_level)

    def node(self, i: int):
        return {"index": i, "atoms": self.G.nodes[i]["atoms"]}

    def fragment_size(self, i: int) -> int:
        return len(self.G.nodes[i]["atoms"])

    def edge(self, i: int, j: int) -> dict:
        key = f"{i}->{j}"
        if key not in self.edge_to_bond:
            warnings.warn(f"edge ({i}, {j}) is not in the fragment quotient graph.")
            return

        return {
            "index": (i, j),
            "atoms": self.edge_to_bond[key],
            "order": self.edge_to_order[key],
        }

    def neighbors(self, i: int) -> list:
        return [n for n in self.G.neighbors(i)]

    @property
    def nodes(self) -> list:
        return [self.node(i) for i in range(self.number_of_nodes)]

    @property
    def edges(self) -> list:
        return [self.edge(e[0], e[1]) for e in self.ordered_edges]

    def fragment(self, i: int):
        return self.G.nodes[i]["atoms"]

    @property
    def fragments(self) -> list:
        return [self.G.nodes[i]["atoms"] for i in range(self.number_of_nodes)]

    @staticmethod
    def quotient_graph_of_rigid_fragments(
        G: nx.Graph, edge_cut: Iterable[Iterable[int]]
    ):
        """
        Generate the fragment quotient graph of the rigid fragments.
        The nodes of the new graph are the rigid fragments, and the edges are those rotatable bonds.
        """

        # remove edges in edge_cut
        G0 = copy.deepcopy(G)
        G0.remove_edges_from(edge_cut)

        # get all connected fragments
        fragments = [
            set([int(x) for x in frag]) for frag in nx.connected_components(G0)
        ]

        # create new quotient graph Q
        # each node in Q is a rigid fragment
        # each edge in Q represents a rotatable bond
        Q = nx.Graph()
        n_frags = len(fragments)
        for i in range(n_frags):
            Q.add_node(i, atoms=sorted(list(fragments[i])))

        directed_edge_to_rot_bond = {}
        for i in range(n_frags):
            for j in range(i + 1, n_frags):
                rot_bond = None
                for edge in edge_cut:
                    edge = [int(x) for x in edge]
                    if edge[0] in fragments[i] and edge[1] in fragments[j]:
                        rot_bond = [edge[0], edge[1]]
                    if edge[1] in fragments[i] and edge[0] in fragments[j]:
                        rot_bond = [edge[1], edge[0]]

                    if rot_bond is not None:
                        # edge is detected
                        break
                if rot_bond is not None:
                    Q.add_edge(i, j)
                    directed_edge_to_rot_bond[f"{i}->{j}"] = [rot_bond[0], rot_bond[1]]
                    directed_edge_to_rot_bond[f"{j}->{i}"] = [rot_bond[1], rot_bond[0]]

        # make sure Q is a Tree
        assert nx.is_forest(
            Q
        ), "The generated graph of rigid fragments is not a forest!"

        return Q, directed_edge_to_rot_bond
