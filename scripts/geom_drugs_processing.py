"""
GEOM-Drugs Dataset Processing and Validation Module

This module implements the complete preprocessing pipeline for the GEOM-Drugs dataset,
including molecular sanitization, topology validation, and valency analysis. It addresses
chemical accuracy issues in the original dataset by applying rigorous validation criteria.

Key Processing Steps:
1. Molecular sanitization using RDKit with full kekulization
2. Topology validation using distance-based bond length checking
3. Fragment filtering to remove disconnected molecules
4. Valency dictionary construction for stability analysis

Topology Validation Strategy:
The module uses covalent radii and tolerance-based distance checking to validate
that molecular geometries are chemically reasonable. This catches issues like:
- Unrealistic bond lengths that would indicate poor conformer generation
- Overlapping atoms or impossible geometries
- Inconsistent connectivity across conformers

Distance Tolerance Logic:
- Uses element-specific covalent radii as reference distances
- Applies 40% tolerance by default (configurable)
- Validates all conformers must have consistent topology
- Removes molecules with any invalid conformers
"""

import json
import os
import pickle
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import SanitizeFlags

# Covalent radii dictionary for topology validation
# Values in Angstroms, used to determine expected bond lengths
covalent_radii = {
    1: 0.31,   # Hydrogen
    6: 0.76,   # Carbon  
    7: 0.71,   # Nitrogen
    8: 0.66,   # Oxygen
    9: 0.57,   # Fluorine
    15: 1.07,  # Phosphorus
    16: 1.05,  # Sulfur
    17: 1.02,  # Chlorine
    35: 1.20,  # Bromine
    53: 1.39   # Iodine
}


def get_reference_distance_matrix(adjacency_matrix, numbers):
    """
    Calculate expected bond lengths based on covalent radii.

    """
    adjacency_mask = (adjacency_matrix > 0).astype(int)
    # Get covalent radii for each atom, default to 1.5 Å for unknown elements
    rads = np.array([covalent_radii.get(i, 1.5) for i in numbers])
    # Sum of covalent radii gives expected bond length
    return rads[:, np.newaxis] + rads[np.newaxis, :]


def calculate_distance_map(coordinates):
    # Broadcasting to compute all pairwise distances efficiently
    diff = coordinates[:, :, np.newaxis, :] - coordinates[:, np.newaxis, :, :]
    return np.linalg.norm(diff, axis=-1)


def check_topology(adjacency_matrix, numbers, coordinates, tolerance=0.4):
    """
    Validate molecular topology using distance-based bond length checking.
    
    Validation Logic:
    1. Calculate expected bond lengths from covalent radii
    2. Measure actual bond lengths in conformer geometries
    3. Check if actual lengths are within tolerance of expected values
    4. Require ALL bonded pairs in ALL conformers to pass validation
    
    Args:
        adjacency_matrix (np.ndarray): Molecular connectivity [n_atoms, n_atoms]
        numbers (np.ndarray): Atomic numbers [n_atoms]
        coordinates (np.ndarray): Conformer coordinates [n_conformers, n_atoms, 3]
        tolerance (float): Relative tolerance for bond length validation (default: 0.4 = 40%)
        
    Returns:
        np.ndarray: Boolean array indicating which conformers pass validation
    """
    # Create mask for bonded atom pairs only
    adjacency_mask = (adjacency_matrix > 0).astype(int)
    
    # Get expected distances and mask to bonded pairs only
    ref_dist = get_reference_distance_matrix(adjacency_matrix, numbers) * adjacency_mask
    
    # Calculate actual distances and mask to bonded pairs only
    data_dist = calculate_distance_map(coordinates) * adjacency_mask
    
    # Check if actual distances are within tolerance of expected distances
    # tolerance is relative: actual_dist <= expected_dist * (1 + tolerance)
    diffs = np.abs(data_dist - ref_dist[np.newaxis, :, :]) <= (
        ref_dist[np.newaxis, :, :] * tolerance)
    
    # All bonded pairs in all conformers must pass validation
    return diffs.all(axis=(1, 2))


def process_molecule(mol):
    """
    Sanitize and kekulize a molecular structure with strict validation.
    
    Sanitization Steps Applied:
    - SANITIZE_ALL: Complete RDKit sanitization pipeline
    - Kekulization with clearAromaticFlags=True: Ensures proper aromatic/Kekulé form
    
    Fragment Filtering:
    - Rejects molecules with multiple disconnected fragments
    - Ensures analysis focuses on single molecular entities
    
    Args:
        mol (rdkit.Chem.Mol): Input RDKit molecule object
        
    Returns:
        rdkit.Chem.Mol or None: Sanitized molecule or None if processing fails
    """
    try:
        # Create a copy to avoid modifying the original
        mol = Chem.Mol(mol)
        
        # Apply full RDKit sanitization pipeline
        Chem.SanitizeMol(mol, sanitizeOps=SanitizeFlags.SANITIZE_ALL)
        
        # Kekulize with aromatic flag clearing for consistent valency analysis
        # This ensures aromatic systems are properly handled in stability calculations
        Chem.Kekulize(mol, clearAromaticFlags=True)
        
    except Exception:
        # Any sanitization failure results in molecule rejection
        return None
        
    # Reject molecules with multiple fragments (disconnected components)
    if len(Chem.GetMolFrags(mol)) > 1:
        return None
        
    return mol


def validate_topology(mol):
    """
    Validate that all conformers of a molecule have reasonable geometries.
    """
    # Extract molecular connectivity
    adjacency_matrix = Chem.GetAdjacencyMatrix(mol)
    
    # Get atomic numbers for covalent radius calculations
    numbers = np.array([atom.GetAtomicNum() for atom in mol.GetAtoms()])
    
    # Extract conformer coordinates
    conformers = mol.GetConformers()
    if not conformers:
        return False  # No conformers = invalid
        
    coordinates = np.array([conf.GetPositions() for conf in conformers])
    
    # Apply topology validation to all conformers
    return check_topology(adjacency_matrix, numbers, coordinates).all()


def process_geom_drugs(input_folder, output_folder):
    """
    Complete processing pipeline for GEOM-Drugs dataset.
    
    Processing Pipeline:
    1. Load raw pickle files from MiDi GEOM-Drugs splits
    2. Validate SMILES strings can be parsed by RDKit
    3. Apply molecular sanitization and kekulization
    4. Perform topology validation on all conformers
    5. Filter out molecules with invalid geometries or connectivity
    6. Build comprehensive valency dictionary for stability analysis
    7. Save cleaned data and analysis results
    
    Args:
        input_folder (str): Path to directory containing raw GEOM-Drugs pickle files
        output_folder (str): Path to save processed dataset and analysis files
        
    Output Files:
        - {train,val,test}_data.pickle: Cleaned molecular datasets
        - valency_dict.json: Comprehensive valency lookup table
        
        
    """
    os.makedirs(output_folder, exist_ok=True)
    valency_dict = {}

    # Process each data split independently
    for split in ["test_data.pickle", "train_data.pickle", "val_data.pickle"]:
        input_path = Path(input_folder) / split
        output_path = Path(output_folder) / split

        if not input_path.exists():
            print(f"Skipping {input_path}, file not found.")
            continue

        # Load raw dataset
        with open(input_path, "rb") as f:
            data = pickle.load(f)

        # Initialize statistics tracking
        initial_size = len(data)
        conformer_count = 0
        initial_conformer_count = 0
        conformer_counts = []

        # Process molecules with in-place filtering
        i = 0
        while i < len(data):
            smiles, mols = data[i]

            # Validate SMILES can be parsed (sanity check)
            reference_mol = Chem.MolFromSmiles(smiles)
            if reference_mol is None:
                del data[i]
                continue

            initial_conformer_count += len(mols)

            # Apply molecular sanitization to all conformers
            sanitized_mols = [process_molecule(mol) for mol in mols]
            sanitized_mols = [mol for mol in sanitized_mols if mol is not None]

            # Require at least one valid conformer and all must pass topology validation
            if not sanitized_mols or not all(validate_topology(mol) for mol in sanitized_mols):
                del data[i]
                continue

            # Build valency dictionary from validated molecules
            # This creates an empirical valency table from the cleaned dataset
            for mol in sanitized_mols:
                for atom in mol.GetAtoms():
                    element = atom.GetSymbol()
                    charge = str(atom.GetFormalCharge())
                    valency = atom.GetExplicitValence()
                    valency_dict.setdefault(element, {}).setdefault(charge, [])
                    if valency not in valency_dict[element][charge]:
                        valency_dict[element][charge].append(valency)
                        
            # Update statistics and store cleaned data
            conformer_counts.append(len(sanitized_mols))
            conformer_count += len(sanitized_mols)
            data[i] = (smiles, sanitized_mols)
            i += 1

        # Calculate and report processing statistics
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

        # Save cleaned dataset
        with open(output_path, "wb") as f:
            pickle.dump(data, f)

    # Save comprehensive valency dictionary for stability analysis
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
