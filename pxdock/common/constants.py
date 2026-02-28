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

import re

CHG_FACTOR = 332.0637141491396
SMALL_NUMBER = 1e-20

RND_VDW_THRESHOLD = 100
RND_COUNT_THRESHOLD = 2000

ELEMENTS_TO_ATOMIC_NUMBERS = {
    "H": 1,
    "He": 2,
    "Li": 3,
    "Be": 4,
    "B": 5,
    "C": 6,
    "N": 7,
    "O": 8,
    "F": 9,
    "Ne": 10,
    "Na": 11,
    "Mg": 12,
    "Al": 13,
    "Si": 14,
    "P": 15,
    "S": 16,
    "Cl": 17,
    "Ar": 18,
    "K": 19,
    "Ca": 20,
    "Sc": 21,
    "Ti": 22,
    "V": 23,
    "Cr": 24,
    "Mn": 25,
    "Fe": 26,
    "Co": 27,
    "Ni": 28,
    "Cu": 29,
    "Zn": 30,
    "Ga": 31,
    "Ge": 32,
    "As": 33,
    "Se": 34,
    "Br": 35,
    "Kr": 36,
    "Rb": 37,
    "Sr": 38,
    "Y": 39,
    "Zr": 40,
    "Nb": 41,
    "Mo": 42,
    "Tc": 43,
    "Ru": 44,
    "Rh": 45,
    "Pd": 46,
    "Ag": 47,
    "Cd": 48,
    "In": 49,
    "Sn": 50,
    "Sb": 51,
    "Te": 52,
    "I": 53,
    "Xe": 54,
    "Cs": 55,
    "Ba": 56,
    "La": 57,
    "Ce": 58,
    "Pr": 59,
    "Nd": 60,
    "Pm": 61,
    "Sm": 62,
    "Eu": 63,
    "Gd": 64,
    "Tb": 65,
    "Dy": 66,
    "Ho": 67,
    "Er": 68,
    "Tm": 69,
    "Yb": 70,
    "Lu": 71,
    "Hf": 72,
    "Ta": 73,
    "W": 74,
    "Re": 75,
    "Os": 76,
    "Ir": 77,
    "Pt": 78,
    "Au": 79,
    "Hg": 80,
    "Tl": 81,
    "Pb": 82,
    "Bi": 83,
    "Th": 90,
    "Pa": 91,
    "U": 92,
    "Np": 93,
    "Pu": 94,
    "Am": 95,
    "Cm": 96,
    "Bk": 97,
    "Cf": 98,
    "Es": 99,
    "Fm": 100,
    "Md": 101,
    "No": 102,
    "Lr": 103,
    "Rf": 104,
    "Db": 105,
    "Sg": 106,
    "Bh": 107,
    "Hs": 108,
    "Mt": 109,
    "Ds": 110,
    "Rg": 111,
    "Cn": 112,
    "Nh": 113,
    "Fl": 114,
    "Mc": 115,
    "Lv": 116,
    "Ts": 117,
    "Og": 118,
}

DEFAULT_RADIUS = [
    1.5,  # NULL
    1.2,  # H
    1.5,  # He
    1.5,  # Li
    1.5,  # Be
    1.5,  # B
    1.7,  # C
    1.55,  # N
    1.5,  # O
    1.5,  # F
    1.5,  # Ne
    1.5,  # Na
    1.5,  # Mg
    1.5,  # Al
    2.1,  # Si
    1.85,  # P
    1.8,  # S
    1.7,  # Cl
]  # in Angstrom
ATOMIC_NUMBERS_TO_ELEMENTS = {i: ai for ai, i in ELEMENTS_TO_ATOMIC_NUMBERS.items()}

DEFAULT_SCREEN = [
    0.8,  # NULL
    0.85,  # H
    0.8,  # He
    0.8,  # Li
    0.8,  # Be
    0.8,  # B
    0.72,  # C
    0.79,  # N
    0.85,  # O
    0.88,  # F
    0.8,  # Ne
    0.8,  # Na
    0.8,  # Mg
    0.8,  # Al
    0.8,  # Si
    0.86,  # P
    0.96,  # S
]

# 用不着的均填充的1.5
VDW_VINA = [
    1.5,  # NULL
    1.0,  # H
    1.4,  # He #
    1.8,  # Li #
    1.53,  # Be #
    1.92,  # B #
    2.0,  # C
    1.75,  # N
    1.60,  # O
    1.545,  # F
    1.5,  # Ne #
    1.75,  # Na
    0.65,  # Mg
    1.5,  # Al #
    2.3,  # Si
    2.1,  # P
    2.0,  # S
    2.045,  # Cl
    1.5,  # Ar #
    1.5,  # K #
    0.99,  # Ca
    1.5,  # Sc
    1.5,  # Ti
    1.5,  # V
    1.5,  # Cr
    0.65,  # Mn
    0.65,  # Fe
    0.65,  # Co #
    0.65,  # Ni #
    0.65,  # Cu
    0.65,  # Zn
    1.5,  # Ga #
    1.5,  # Ge #
    1.5,  # As #
    1.5,  # Se #
    2.165,  # Br
]


def mapped_smiles_to_atoms(mapped_smiles: str) -> list[str]:
    matches = re.findall(r"\[([A-Za-z@\+\-]+):(\d+)", mapped_smiles)
    sorted_atoms = [
        match[0].rstrip("+-@").capitalize()
        for match in sorted(matches, key=lambda x: int(x[1]))
    ]
    return sorted_atoms


def mapped_smiles_to_atomic_numbers(mapped_smiles: str) -> list[int]:
    sorted_atoms = mapped_smiles_to_atoms(mapped_smiles)
    return [ELEMENTS_TO_ATOMIC_NUMBERS[x] for x in sorted_atoms]


if __name__ == "__main__":
    mapped_smiles = "[C:1]([c:2]1[c:3]([H:32])[c:4]([N:5]2[C:6]([H:33])([H:34])[C:7]([H:35])([H:36])[C:8]([H:37])([H:38])[C@:9]([c:11]3[n:12][n:13]([H:28])[c:14]([C:16](=[O:17])[N:18]([H:43])[H:44])[c:15]3[H:42])([H:39])[C:10]2([H:40])[H:41])[c:19]2[c:20]([n+:21]1[H:47])[c:22]([F:23])[c:24]([F:27])[c:25]([H:45])[c:26]2[H:46])([H:29])([H:30])[H:31]"
    atoms = mapped_smiles_to_atoms(mapped_smiles)
    print(atoms)
    print(len(atoms))
