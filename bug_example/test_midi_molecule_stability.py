import torch

from midi_molecule_stability import Molecule, check_stability

class _DatasetInfo:
    def __init__(self, atom_decoder):
        self.atom_decoder = atom_decoder

def _zeros_xyz(n):
    return torch.zeros((n, 3), dtype=torch.float)

def test_ch2_n_ch2():
    """
    Structure: CH2–N–CH2 (no hydrogen on N).
    Objective: Demonstrate that a neutral nitrogen with valence=2, and carbon with valence=3 passes
    stability per the current allowed_bonds ('N': {0: [2, 3]}).
    """
    # Decoder order: 0='C', 1='H', 2='N'
    atom_decoder = ['C', 'H', 'N']
    dsinfo = _DatasetInfo(atom_decoder)

    # Indices:
    # 0 = C(left), 1 = H(L1), 2 = H(L2),
    # 3 = N,       4 = C(right), 5 = H(R1), 6 = H(R2)
    atom_types = torch.tensor([0, 1, 1, 2, 0, 1, 1], dtype=torch.long)
    charges    = torch.zeros(7, dtype=torch.long)

    # Initialize 7x7 bond matrix (LongTensor)
    bonds = torch.zeros((7, 7), dtype=torch.long)

    # Left CH2: C(0)–H(1), C(0)–H(2)
    bonds[0, 1] = bonds[1, 0] = 1
    bonds[0, 2] = bonds[2, 0] = 1

    # Backbone: C(0)–N(3)–C(4)
    bonds[0, 3] = bonds[3, 0] = 1
    bonds[3, 4] = bonds[4, 3] = 1

    # Right CH2: C(4)–H(5), C(4)–H(6)
    bonds[4, 5] = bonds[5, 4] = 1
    bonds[4, 6] = bonds[6, 4] = 1

    positions = _zeros_xyz(7)

    # Build Molecule via your original class
    mol = Molecule(
        atom_types=atom_types,
        bond_types=bonds,
        positions=positions,
        charges=charges,
        atom_decoder=atom_decoder,
    )

    # Run original stability check
    stable, n_stable, n_atoms = check_stability(mol, dsinfo)

    # Compute integer valencies as used by the function (no aromatic bonds here)
    valencies = bonds.sum(dim=-1)
    n_val = int(valencies[3].item())  # N is index 3

    # Assertions for the demonstration:
    # - EXACTLY 7 atoms
    assert n_atoms == 7
    # - Neutral N has valence 2 in this construction (two single bonds to carbons)
    assert n_val == 2
    # - Per allowed_bonds: 'N': {0: [2, 3]} → valence 2 is allowed for neutral N → stable
    assert bool(stable.item()) is True, (
        "Neutral N with valence=2 passes the current stability rules for this CH2–N–CH2 example."
    )

    def test_ch2_n_ch2():
        """
        Structure: CH2–N–CH2 (no hydrogen on N).
        Objective: Demonstrate that a neutral nitrogen with valence=2, and carbon with valence=3 passes
        stability per the current allowed_bonds ('N': {0: [2, 3]}).
        """
        # Decoder order: 0='C', 1='H', 2='N'
        atom_decoder = ['C', 'H', 'N']
        dsinfo = _DatasetInfo(atom_decoder)

        # Indices:
        # 0 = C(left), 1 = H(L1), 2 = H(L2),
        # 3 = N,       4 = C(right), 5 = H(R1), 6 = H(R2)
        atom_types = torch.tensor([0, 1, 1, 2, 0, 1, 1], dtype=torch.long)
        charges = torch.zeros(7, dtype=torch.long)

        # Initialize 7x7 bond matrix (LongTensor)
        bonds = torch.zeros((7, 7), dtype=torch.long)

        # Left CH2: C(0)–H(1), C(0)–H(2)
        bonds[0, 1] = bonds[1, 0] = 1
        bonds[0, 2] = bonds[2, 0] = 1

        # Backbone: C(0)–N(3)–C(4)
        bonds[0, 3] = bonds[3, 0] = 4
        bonds[3, 4] = bonds[4, 3] = 4

        # Right CH2: C(4)–H(5), C(4)–H(6)
        bonds[4, 5] = bonds[5, 4] = 1
        bonds[4, 6] = bonds[6, 4] = 1

        positions = _zeros_xyz(7)

        # Build Molecule via your original class
        mol = Molecule(
            atom_types=atom_types,
            bond_types=bonds,
            positions=positions,
            charges=charges,
            atom_decoder=atom_decoder,
        )

        # Run original stability check
        stable, n_stable, n_atoms = check_stability(mol, dsinfo)

        # Compute integer valencies as used by the function (no aromatic bonds here)
        valencies = bonds.sum(dim=-1)
        n_val = int(valencies[3].item())  # N is index 3

        # Assertions for the demonstration:
        # - EXACTLY 7 atoms
        assert n_atoms == 7
        # - Neutral N has valence 2 in this construction (two single bonds to carbons)
        assert n_val == 2
        # - Per allowed_bonds: 'N': {0: [2, 3]} → valence 2 is allowed for neutral N → stable
        assert bool(stable.item()) is True, (
            "Neutral N with valence=2 passes the current stability rules for this CH2–N–CH2 example."
        )


def test_ch2_n_ch2():
    """
    Structure: CH2–N–CH2 (no hydrogen on N).
    Objective: Demonstrate that a neutral nitrogen with valence=2, and carbon with valence=3 passes
    stability per the current allowed_bonds ('N': {0: [2, 3]}).
    """
    # Decoder order: 0='C', 1='H', 2='N'
    atom_decoder = ['C', 'H', 'N']
    dsinfo = _DatasetInfo(atom_decoder)

    # Indices:
    # 0 = C(left), 1 = H(L1), 2 = H(L2),
    # 3 = N,       4 = C(right), 5 = H(R1), 6 = H(R2)
    atom_types = torch.tensor([0, 1, 1, 2, 0, 1, 1], dtype=torch.long)
    charges    = torch.zeros(7, dtype=torch.long)

    # Initialize 7x7 bond matrix (LongTensor)
    bonds = torch.zeros((7, 7), dtype=torch.long)

    # Left CH2: C(0)–H(1), C(0)–H(2)
    bonds[0, 1] = bonds[1, 0] = 1
    bonds[0, 2] = bonds[2, 0] = 1

    # Backbone: C(0)–N(3)–C(4)
    bonds[0, 3] = bonds[3, 0] = 1
    bonds[3, 4] = bonds[4, 3] = 1

    # Right CH2: C(4)–H(5), C(4)–H(6)
    bonds[4, 5] = bonds[5, 4] = 1
    bonds[4, 6] = bonds[6, 4] = 1

    positions = _zeros_xyz(7)

    # Build Molecule via your original class
    mol = Molecule(
        atom_types=atom_types,
        bond_types=bonds,
        positions=positions,
        charges=charges,
        atom_decoder=atom_decoder,
    )

    # Run original stability check
    stable, n_stable, n_atoms = check_stability(mol, dsinfo)

    # Compute integer valencies as used by the function (no aromatic bonds here)
    valencies = bonds.sum(dim=-1)
    n_val = int(valencies[3].item())  # N is index 3

    # Assertions for the demonstration:
    # - EXACTLY 7 atoms
    assert n_atoms == 7
    # - Neutral N has valence 2 in this construction (two single bonds to carbons)
    assert n_val == 2
    # - Per allowed_bonds: 'N': {0: [2, 3]} → valence 2 is allowed for neutral N → stable
    assert bool(stable.item()) is True, (
        "Neutral N with valence=2 passes the current stability rules for this CH2–N–CH2 example."
    )


def test_ch3_ch3():
    """
    Structure: CH3:CH3
    We encode the C–C bond as '4' to indicate aromatic-like (intended 1.5).
    Intended valency at each carbon: 3 (C–H) + 1.5 (C–C_aromatic) = 4.5  -> NOT allowed for neutral carbon {3,4}
    Bugged valency (because 1.5 truncated to 1): 3 + 1 = 4 -> allowed -> function returns STABLE
    This test proves the aromatic replacement bug flips the classification.
    """
    # Decoder: 0='C', 1='H'
    atom_decoder = ['C', 'H']
    dsinfo = _DatasetInfo(atom_decoder)

    # Indices:
    # 0=C(L), 1=H(L1), 2=H(L2), 3=H(L3),
    # 4=C(R), 5=H(R1), 6=H(R2), 7=H(R3)
    atom_types = torch.tensor([0, 1, 1, 1, 0, 1, 1, 1], dtype=torch.long)
    charges    = torch.zeros(8, dtype=torch.long)

    bonds = torch.zeros((8, 8), dtype=torch.long)

    # Left CH3
    bonds[0, 1] = bonds[1, 0] = 1
    bonds[0, 2] = bonds[2, 0] = 1
    bonds[0, 3] = bonds[3, 0] = 1

    # Right CH3
    bonds[4, 5] = bonds[5, 4] = 1
    bonds[4, 6] = bonds[6, 4] = 1
    bonds[4, 7] = bonds[7, 4] = 1

    # C–C "aromatic" placeholder bond encoded as 4
    bonds[0, 4] = bonds[4, 0] = 4

    bonds_copy = bonds.clone()

    positions = _zeros_xyz(8)

    mol = Molecule(
        atom_types=atom_types,
        bond_types=bonds,
        positions=positions,
        charges=charges,
        atom_decoder=atom_decoder,
    )

    # Run original (bugged) stability check
    stable_bug, _, _ = check_stability(mol, dsinfo)

    # Bugged integer path: 4 -> 1, so each carbon sees valency 4 (3 H + 1 C) which IS allowed
    assert bool(stable_bug.item()) is True, (
        "Bugged function misclassifies CH3–(aromatic)–CH3 as stable because 1.5 was truncated to 1."
    )

    bonds_copy[bonds_copy == 4] = 1.5

    assert 1.5 not in bonds_copy

