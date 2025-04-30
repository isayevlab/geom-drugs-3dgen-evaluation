"""
Valency lookup table for elements found in the GEOM-Drugs dataset, including hydrogen counts.

This dictionary defines chemically plausible valency configurations for atoms
based on their formal charges. Each element is mapped to a nested dictionary
where:
  - keys are formal charges (int),
  - values are lists of allowed valencies (number of bonds including implicit Hs).

The table was constructed from the cleaned and validated GEOM-Drugs dataset and
is used for evaluating molecular stability in generative 3D molecule models.

Usage:
    valencies = geom_drugs_h_valencies[element][formal_charge]
"""

geom_drugs_h_valencies = {
    "Br": {0: [ 1], 1: [2]},
    "C": {0: [ 4 ], -1: [ 3 ], 1: [3]},
    "N": {0: [ 3 ], 1: [ 4 ], -1: [ 2 ], -2: [ 1 ] },
    "H": { 0: [ 1 ]},
    "S": { 0: [ 2, 6, 3], 1: [ 3 ], 2: [ 4], 3: [ 5, 2 ], -1: [ 1]},
    "O": { 0: [ 2], -1: [ 1 ], 1: [ 3] },
    "F": { 0: [ 1 ]},
    "Cl": { 0: [ 1 ], 1: [ 2 ] },
    "P": { 0: [ 5, 3 ], 1: [ 4 ] },
    "I": { 0: [ 1 ], 1: [ 2 ], 2: [ 3 ]},
    "Si": { 0: [ 4 ], 1: [ 5 ] },
    "B": { -1: [ 4 ], 0: [ 3 ] },
    "Bi": { 2: [ 5 ],  0: [ 3]
    }
}


geom_drugs_h_bad_valencies = {
    "H": {0: 1, 1: 0, -1: 0},
    "C": {0: [3, 4], 1: 3, -1: 3},
    "N": {0: [2, 3], 1: [2, 3, 4], -1: 2},
    "O": {0: 2, 1: 3, -1: 1},
    "F": {0: 1, -1: 0},
    "B": 3,
    "Al": 3,
    "Si": 4,
    "P": {0: [3, 5], 1: 4},
    "S": {0: [2, 6], 1: [2, 3], 2: 4, 3: 5, -1: 3},
    "Cl": 1,
    "As": 3,
    "Br": {0: 1, 1: 2},
    "I": 1,
    "Hg": [1, 2],
    "Bi": [3, 5],
    "Se": [2, 4, 6],
}


geom_drugs_h_tuple_valencies = {
    "Br": {
        0: [(0, 1)],
        1: [(0, 2)],
    },
    "C": {
        0: [(0, 4), (2, 2), (2, 1), (3, 0)],
        -1: [(0, 3), (2, 1), (3, 0)],
        1: [(0, 3), (2, 1), (3, 0)],
    },
    "N": {
        0: [(0, 3), (2, 0), (2, 1), (3, 0)],
        1: [(0, 4), (2, 0), (2, 1), (2, 2), (3, 0)],
        -1: [(0, 2), (2, 0)],
        -2: [(0, 1)],
    },
    "H": {
        0: [(0, 1)],
    },
    "S": {
        0: [(0, 2), (0, 3), (0, 6), (2, 0)],
        1: [(0, 3), (2, 0), (2, 1), (3, 0)],
        2: [(0, 4), (2, 1), (2, 2)],
        3: [(0, 2), (0, 5)],
        -1: [(0, 1)],
    },
    "O": {
        0: [(0, 2), (2, 0)],
        -1: [(0, 1)],
        1: [(0, 3)],
    },
    "F": {
        0: [(0, 1)],
    },
    "Cl": {
        0: [(0, 1)],
        1: [(0, 2)],
    },
    "P": {
        0: [(0, 3), (0, 5)],
        1: [(0, 4)],
    },
    "I": {
        0: [(0, 1)],
        1: [(0, 2)],
        2: [(0, 3)],
    },
    "Si": {
        0: [(0, 4)],
        1: [(0, 5)],
    },
    "B": {
        -1: [(0, 4)],
        0: [(0, 3)],
    },
    "Bi": {
        0: [(0, 3)],
        2: [(0, 5)],
    }
}
