import json
import os
import pickle
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import SanitizeFlags

# Covalent radii dictionary
covalent_radii = {
    1: 0.31, 6: 0.76, 7: 0.71, 8: 0.66, 9: 0.57, 15: 1.07, 16: 1.05, 17: 1.02,
    35: 1.20, 53: 1.39
}


def get_reference_distance_matrix(adjacency_matrix, numbers):
    adjacency_mask = (adjacency_matrix > 0).astype(int)
    rads = np.array([covalent_radii.get(i, 1.5) for i in numbers])
    return rads[:, np.newaxis] + rads[np.newaxis, :]


def calculate_distance_map(coordinates):
    diff = coordinates[:, :, np.newaxis, :] - coordinates[:, np.newaxis, :, :]
    return np.linalg.norm(diff, axis=-1)


def check_topology(adjacency_matrix, numbers, coordinates, tolerance=0.4):
    adjacency_mask = (adjacency_matrix > 0).astype(int)
    ref_dist = get_reference_distance_matrix(adjacency_matrix, numbers) * adjacency_mask
    data_dist = calculate_distance_map(coordinates) * adjacency_mask
    diffs = np.abs(data_dist - ref_dist[np.newaxis, :, :]) <= (
        ref_dist[np.newaxis, :, :] * tolerance)
    return diffs.all(axis=(1, 2))


def process_molecule(mol):
    try:
        mol = Chem.Mol(mol)
        Chem.SanitizeMol(mol, sanitizeOps=SanitizeFlags.SANITIZE_ALL)
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Exception:
        return None
    if len(Chem.GetMolFrags(mol)) > 1:
        return None
    return mol


def validate_topology(mol):
    adjacency_matrix = Chem.GetAdjacencyMatrix(mol)
    numbers = np.array([atom.GetAtomicNum() for atom in mol.GetAtoms()])
    conformers = mol.GetConformers()
    if not conformers:
        return False
    coordinates = np.array([conf.GetPositions() for conf in conformers])
    return check_topology(adjacency_matrix, numbers, coordinates).all()


def process_geom_drugs(input_folder, output_folder):
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

        initial_size = len(data)
        conformer_count = 0
        initial_conformer_count = 0
        conformer_counts = []

        i = 0
        while i < len(data):
            smiles, mols = data[i]

            reference_mol = Chem.MolFromSmiles(smiles)
            if reference_mol is None:
                del data[i]
                continue

            initial_conformer_count += len(mols)

            sanitized_mols = [process_molecule(mol) for mol in mols]
            sanitized_mols = [mol for mol in sanitized_mols if mol is not None]

            if not sanitized_mols or not all(validate_topology(mol) for mol in sanitized_mols):
                del data[i]
                continue

            for mol in sanitized_mols:
                for atom in mol.GetAtoms():
                    element = atom.GetSymbol()
                    charge = str(atom.GetFormalCharge())
                    valency = atom.GetExplicitValence()
                    valency_dict.setdefault(element, {}).setdefault(charge, [])
                    if valency not in valency_dict[element][charge]:
                        valency_dict[element][charge].append(valency)
            conformer_counts.append(len(sanitized_mols))
            conformer_count += len(sanitized_mols)
            data[i] = (smiles, sanitized_mols)
            i += 1

        final_size = len(data)
        removed_molecules = initial_size - final_size
        dropped_conformers = initial_conformer_count - sum(conformer_counts)

        print(f"\nProcessed {split}:")
        print(f"  Molecules saved: {final_size}")
        print(f"  Molecules removed: {removed_molecules} ({(removed_molecules / initial_size * 100):.2f}%)")
        print(f"  Kept conformers: {sum(conformer_counts)}")
        print(f"  Dropped conformers: {dropped_conformers} ({(dropped_conformers / initial_conformer_count * 100):.2f}%)")

        if conformer_counts:
            print(f"  Conformer count per molecule: min={min(conformer_counts)}, max={max(conformer_counts)}, mean={np.mean(conformer_counts):.2f}")

        with open(output_path, "wb") as f:
            pickle.dump(data, f)

    valency_output = Path(output_folder) / "valency_dict.json"
    with open(valency_output, "w") as f:
        json.dump(valency_dict, f, indent=4)

    print(f"\nValency dictionary saved to {valency_output}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Process GEOM-Drugs dataset.")
    parser.add_argument("--input_folder", type=str, required=True, help="Path to midi_split folder.")
    parser.add_argument("--output_folder", type=str, required=True, help="Path to save filtered results.")
    args = parser.parse_args()

    process_geom_drugs(args.input_folder, args.output_folder)
