import argparse
from collections import defaultdict

import numpy as np
from rdkit import Chem

from geom_utils.pair_geometry import compute_bond_angles_diff
from geom_utils.pair_geometry import compute_bond_lengths_diff
from geom_utils.pair_geometry import compute_torsion_angles_diff
from geom_utils.utils import bond_type_to_symbol, compute_differences


def print_results(results):
    """
    Print average difference and weight for each pattern, and the total weighted average.
    """
    periodic_table = Chem.GetPeriodicTable()
    sorted_results = sorted(results.items(), key=lambda x: np.mean(x[1][1]), reverse=True)

    total_weighted_diff = 0.0
    total_weight = 0.0

    print("Final Results Across All Files (sorted by weight):")

    for key, (avg_diff_list, weight_list) in sorted_results:
        mean_diff = np.mean(avg_diff_list)
        mean_weight = np.mean(weight_list)

        atoms = key[::2]
        bonds = key[1::2]
        atom_symbols = [periodic_table.GetElementSymbol(a) for a in atoms]
        bond_symbols = [bond_type_to_symbol(b) for b in bonds]

        bond_repr = ""
        for i, atom in enumerate(atom_symbols):
            bond_repr += atom
            if i < len(bond_symbols):
                bond_repr += bond_symbols[i]

        if mean_weight > 0.005:
            print(f"{bond_repr}: Mean = {mean_diff:.4f}, Weight = {mean_weight:.4f}")

        total_weighted_diff += mean_diff * mean_weight
        total_weight += mean_weight

    if total_weight > 0:
        print(f"\nTotal Weighted Mean Difference: {total_weighted_diff / total_weight:.4f}")
    else:
        print("\nNo valid results to report.")


def run_analysis(pairs, analysis_type="bond_length"):
    accumulated_results = defaultdict(lambda: ([], []))

    is_not_none = lambda x: (x[0] is not None) and (x[1] is not None)
    pairs = list(filter(is_not_none, pairs))

    if analysis_type == "bond_length":
        results = compute_differences(pairs, compute_bond_lengths_diff)
    elif analysis_type == "bond_angle":
        results = compute_differences(pairs, compute_bond_angles_diff)
    elif analysis_type == "torsion_angle":
        results = compute_differences(pairs, compute_torsion_angles_diff)
    else:
        raise ValueError(f"Unknown analysis type: {analysis_type}")

    accumulate_results(accumulated_results, results)
    return accumulated_results


def accumulate_results(accumulated_results, new_results):
    for key, (avg_diff_list, _, weight_list) in new_results.items():
        accumulated_results[key][0].append(avg_diff_list)
        accumulated_results[key][1].append(weight_list)


def summarize_results(results):
    periodic_table = Chem.GetPeriodicTable()
    total_weighted_diffs = []

    for key, (avg_diff_list, weight_list) in results.items():
        weighted_diffs = np.array(avg_diff_list) * np.array(weight_list)
        total_weighted_diffs.append(np.sum(weighted_diffs))

    return np.sum(total_weighted_diffs) / sum(np.sum(results[k][1]) for k in results)


def run_subsets_analysis(pairs, analysis_type, n_subsets):
    fold_size = len(pairs) // n_subsets
    scores = []

    for i in range(n_subsets):
        fold_pairs = pairs[i * fold_size: (i + 1) * fold_size] if i < n_subsets - 1 else pairs[
                                                                                         i * fold_size:]
        results = run_analysis(fold_pairs, analysis_type=analysis_type)
        score = summarize_results(results)
        scores.append(score)

    return scores


def print_scores(name, scores):
    mean = np.mean(scores)
    std = np.std(scores)
    print(f"{name}: {mean:.4f} ± {std:.4f}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--init_sdf", type=str, required=True,
                        help="Path to initial structure SDF file.")
    parser.add_argument("--opt_sdf", type=str, required=True,
                        help="Path to optimized structure SDF file.")
    parser.add_argument("--hydrogens", action="store_true", default=True,
                        help="Include hydrogens in geometry calculations.")
    parser.add_argument("--n_subsets", type=int, default=1,
                        help="Number of subsets to compute std over.")
    args = parser.parse_args()

    # Read molecules from SDF files
    suppl1 = Chem.SDMolSupplier(args.init_sdf, sanitize=False, removeHs=not args.hydrogens)
    suppl2 = Chem.SDMolSupplier(args.opt_sdf, sanitize=False, removeHs=not args.hydrogens)

    pairs = [(m1, m2) for m1, m2 in zip(suppl1, suppl2) if m1 is not None and m2 is not None]

    if args.n_subsets == 1:
        print("Bond Lengths")
        lengths = run_analysis(pairs, analysis_type="bond_length")
        print_results(lengths)

        print("\nBond Angles")
        angles = run_analysis(pairs, analysis_type="bond_angle")
        print_results(angles)

        print("\nTorsions")
        torsions = run_analysis(pairs, analysis_type="torsion_angle")
        print_results(torsions)
    else:
        print(f"Running analysis in {args.n_subsets} subsets...\n")

        bond_length_scores = run_subsets_analysis(pairs, "bond_length", args.n_subsets)
        bond_angle_scores = run_subsets_analysis(pairs, "bond_angle", args.n_subsets)
        torsion_scores = run_subsets_analysis(pairs, "torsion_angle", args.n_subsets)

        print_scores("Bond Lengths", bond_length_scores)
        print_scores("Bond Angles", bond_angle_scores)
        print_scores("Torsions", torsion_scores)
