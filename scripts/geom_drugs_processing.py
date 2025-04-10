import json
import os
import pickle
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import SanitizeFlags

# Covalent radii dictionary
covalent_radii = {
    1: 0.31, 6: 0.76, 7: 0.71, 8: 0.66, 9: 0.57, 15: 1.07, 16: 1.05, 17: 1.02, 35: 1.20, 53: 1.39
}


def get_reference_distance_matrix(adjacency_matrix, numbers):
    """Compute reference bond distance matrix based on covalent radii."""
    adjacency_mask = (adjacency_matrix > 0).astype(int)
    rads = np.array([covalent_radii.get(i, 1.5) for i in numbers])  # Default to 1.5 Å if unknown
    return rads[:, np.newaxis] + rads[np.newaxis, :]


def calculate_distance_map(coordinates):
    """Compute distance matrix from 3D coordinates."""
    diff = coordinates[:, :, np.newaxis, :] - coordinates[:, np.newaxis, :, :]
    return np.linalg.norm(diff, axis=-1)


def check_topology(adjacency_matrix, numbers, coordinates, tolerance=0.4):
    """Validate molecular topology based on expected bond distances."""
    num_conformers = coordinates.shape[0]
    adjacency_mask = (adjacency_matrix > 0).astype(int)
    ref_dist = get_reference_distance_matrix(adjacency_matrix, numbers) * adjacency_mask
    data_dist = calculate_distance_map(coordinates) * adjacency_mask

    diffs = np.abs(data_dist - ref_dist[np.newaxis, :, :]) <= (
            ref_dist[np.newaxis, :, :] * tolerance)
    return diffs.all(axis=(1, 2))


def process_molecule(mol):
    """Sanitize, kekulize, and validate a molecule."""
    try:
        mol = Chem.Mol(mol)  # Create a new instance to avoid modifying original data
        Chem.SanitizeMol(mol, sanitizeOps=SanitizeFlags.SANITIZE_ALL)
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Exception:
        return None  # Invalid molecule

    # Ensure it has a single connectivity component
    if len(Chem.GetMolFrags(mol)) > 1:
        return None

    return mol


def validate_topology(mol):
    """Ensure all conformers of a molecule have consistent topology."""
    adjacency_matrix = Chem.GetAdjacencyMatrix(mol)
    numbers = np.array([atom.GetAtomicNum() for atom in mol.GetAtoms()])

    conformers = mol.GetConformers()
    if not conformers:
        return False

    coordinates = np.array([conf.GetPositions() for conf in conformers])
    return check_topology(adjacency_matrix, numbers, coordinates).all()


def process_geom_drugs(input_folder, output_folder):
    """Process all GEOM-Drugs split files, applying the required filters."""
    os.makedirs(output_folder, exist_ok=True)
    valency_dict = {}

    for split in ["test_data.pickle", "train_data.pickle", "val_data.pickle"]:
        input_path = Path(input_folder) / split
        output_path = Path(output_folder) / split

        if not input_path.exists():
            print(f"Skipping {input_path}, file not found.")
            continue

        with open(input_path, "rb") as f:
            data = pickle.load(f)

        initial_size = len(data)  # Track initial molecule count
        conformer_count = 0

        i = 0  # Use manual indexing to allow in-place modification
        while i < len(data):
            smiles, mols = data[i]

            # Convert input SMILES to an RDKit molecule & compute canonical SMILES
            reference_mol = Chem.MolFromSmiles(smiles)
            if reference_mol is None:
                del data[i]
                continue

            # Process & sanitize molecules
            sanitized_mols = [process_molecule(mol) for mol in mols]
            sanitized_mols = [mol for mol in sanitized_mols if
                              mol is not None]  # Remove None entries

            if not sanitized_mols:
                del data[i]
                continue

            # Check topology
            if not all(validate_topology(mol) for mol in sanitized_mols):
                del data[i]
                continue

            # Compute valency dictionary
            for mol in sanitized_mols:
                for atom in mol.GetAtoms():
                    element = atom.GetSymbol()
                    charge = str(atom.GetFormalCharge())
                    valency = atom.GetExplicitValence()

                    if element not in valency_dict:
                        valency_dict[element] = {}

                    if charge not in valency_dict[element]:
                        valency_dict[element][charge] = []

                    if valency not in valency_dict[element][charge]:
                        valency_dict[element][charge].append(valency)

            # Track total number of conformers
            conformer_count += sum(mol.GetNumConformers() for mol in sanitized_mols)

            # If molecule passes all filters, keep it
            data[i] = (smiles, sanitized_mols)
            i += 1  # Increment index only if molecule is kept

        final_size = len(data)  # Track final molecule count
        removed_molecules = initial_size - final_size
        removal_percentage = (removed_molecules / initial_size) * 100 if initial_size > 0 else 0

        print(f"Processed {split}: {final_size} molecules saved.")
        print(f"Removed {removed_molecules} molecules ({removal_percentage:.2f}% of dataset).")

        print(f"Total number of conformers in dataset: {conformer_count}")

        # Save modified dataset
        with open(output_path, "wb") as f:
            pickle.dump(data, f)

    # Save valency dictionary
    valency_output = Path(output_folder) / "valency_dict.json"
    with open(valency_output, "w") as f:
        json.dump(valency_dict, f, indent=4)

    print(f"Valency dictionary saved to {valency_output}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Process GEOM-Drugs dataset.")
    parser.add_argument("--input_folder", type=str, required=True,
                        help="Path to midi_split folder.")
    parser.add_argument("--output_folder", type=str, required=True,
                        help="Path to save filtered results.")

    args = parser.parse_args()
    process_geom_drugs(args.input_folder, args.output_folder)
