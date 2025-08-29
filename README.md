# GEOM-Drugs Revisited  
**Toward More Chemically Accurate Benchmarks for 3D Molecule Generation**

This repository accompanies the paper **"GEOM-Drugs Revisited"**, which introduces a corrected evaluation pipeline for 3D molecular generative models based on GEOM-Drugs. It includes scripts for preprocessing, stability assessment, geometric analysis, and energy-based benchmarking using GFN2-xTB.

---

## 📦 Dependencies

Make sure the following Python packages are installed:

- `rdkit`
- `numpy`
- `tqdm`
- `torch`

External tools:
- [`xtb`](https://github.com/grimme-lab/xtb) (for geometry optimization via GFN2-xTB)

Install dependencies:
```bash
conda config --add channels conda-forge
conda install xtb
pip install numpy tqdm rdkit torch
```

---

## 🗂 Directory Structure

```
.
├──scripts
    ├──geom_drugs_processing.py
    ├── compute_molecule_stability.py
    ├── energy_benchmark/
       ├── compute_pair_geometry.py
       ├── rmsd_energy.py
       └── xtb_optimization.py
└── geom_utils/
   ├── molecule_stability.py
   ├── pair_geometry.py
   ├── geom_drugs_valency_table.py
   └── utils.py
```

---

## 🧪 Script Descriptions

### 1. `geom_drugs_processing.py`

**Purpose:** Sanitize and filter the MiDi GEOM-Drugs split files. Removes disconnected molecules and invalid valencies. Also builds a valency lookup dictionary.

This script takes as input a filepath to the GEOM-Drugs splits from [MiDi](https://github.com/cvignac/MiDi); these can be downloaded using the following command:

```console
wget -r -np -nH --cut-dirs=2 --reject "index.html*" https://bits.csb.pitt.edu/files/geom_raw/
```

**Usage:**
```bash
python geom_drugs_processing.py \
  --input_folder path/to/midi_split \
  --output_folder path/to/cleaned_split
```

**Outputs:**
- Filtered `train/val/test_data.pickle`
- `valency_dict.json` summarizing observed valency configurations

---

### 2. `compute_molecule_stability.py`

**Purpose:** Compute molecular and atom stability of molecules using correct valency tables.

**Usage:**
```bash
python compute_molecule_stability.py \
  --sdf path/to/structures.sdf \
  --aromatic auto \
  --n_subsets 5
```

**Arguments:**
- `--sdf`: Path to an SDF file
- `--aromatic`: `true`, `false`, or `auto` to control aromatic bond handling
- `--n_subsets`: If set, computes mean ± std across subsets

---

### 3. `energy_benchmark/xtb_optimization.py`

**Purpose:** Optimize molecules with GFN2-xTB and extract energy/RMSD values.

**Usage:**
```bash
MKL_NUM_THREADS=16 OMP_NUM_THREADS=16 python energy_benchmark/xtb_optimization.py \
  --input_sdf path/to/generated.sdf \
  --output_sdf path/to/optimized_output.sdf \
  --init_sdf path/to/saved_initial_structures.sdf
```

**Note:** Requires `xtb` to be installed and available in your `PATH`.

---

### 4. `energy_benchmark/compute_pair_geometry.py`

**Purpose:** Compare initial and optimized structures by computing changes in bond lengths, bond angles, and torsion angles.

**Usage:**
```bash
python energy_benchmark/compute_pair_geometry.py \
  --init_sdf path/to/initial.sdf \
  --opt_sdf path/to/optimized.sdf \
  --n_subsets 5
```

---

### 5. `energy_benchmark/rmsd_energy.py`

**Purpose:** Compute GFN2-xTB energy gains, MMFF energy drops, and RMSD across molecule pairs.

**Usage:**
```bash
python energy_benchmark/rmsd_energy.py \
  --init_sdf path/to/initial.sdf \
  --opt_sdf path/to/optimized.sdf \
  --n_subsets 5
```

---

## 📄 Alternative Data Formats

For improved reproducibility and easier maintenance, valency tables are also available in JSON format with schema validation:

### Valency Tables (JSON Format)
```
valency_tables/
├── geom_drugs_h_tuple_valencies.json     # Main aromatic-aware valency table
├── valency_schema.json                   # JSON schema for validation
├── legacy_valencies.json                 # Legacy compatibility table
└── legacy_valency_schema.json           # Schema for legacy format
```

### Usage with JSON Tables
```python
from geom_utils.valency_loader import load_valency_table

# Load main tuple-based valency table
valency_table = load_valency_table('tuple', validate_schema=True)

# Load legacy table for compatibility
legacy_table = load_valency_table('legacy', validate_schema=True)

# Use with stability computation
from geom_utils.molecule_stability import compute_molecules_stability
stability = compute_molecules_stability(molecules, allowed_bonds=valency_table)
```

## 📖 Citation
```
@article{nikitin2025geom,
  title={GEOM-Drugs Revisited: Toward More Chemically Accurate Benchmarks for 3D Molecule Generation},
  author={Nikitin, Filipp and Dunn, Ian and Koes, David Ryan and Isayev, Olexandr},
  journal={arXiv preprint arXiv:2505.00169},
  year={2025}
}
```

