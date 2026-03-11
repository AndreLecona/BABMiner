# BABMiner: β-α-β Protein Motif Tool

A Python tool for identifying and extracting β-α-β (beta-alpha-beta) protein motifs from PDB files.

## Features

- **Motif detection**: Identifies beta-alpha-beta structural patterns in proteins
- **Parallel processing**: Process multiple PDB files concurrently
- **Fragment extraction**: Automatically extracts and saves motif fragments:
  - PDB structure files (`fragments_pdb/`)
  - FASTA sequences (`fragments_fasta/`)
  - Metadata CSV (`metadata.csv`)
- **Rich analysis**: Secondary structure, hydrogen bonding, strand orientation, handedness, sequence coverage

## Installation

```bash
git clone https://github.com/yourusername/BABMiner.git
cd BABMiner
pip install -e .
```

## Quick Start
### Analyze a single PDB file

```bash
BABMiner protein.pdb --outdir results
```

### Analyze a multiple PDB files in parallel

```bash
BABMiner /path/to/pdbs --outdir results --processes 8
```

## Command-Line Options

```
positional arguments:
  input                 PDB file or directory containing PDB files

optional arguments:
  -o, --outdir DIR      Output directory (required)
  -p, --processes N     Parallel processes (default: auto-detect)
  --angle DEGREES       Angle threshold (default: 40)
  --min_hb N           Minimum H-bonds (default: 2)
  -v, --verbose        Enable verbose output
  -h, --help           Show help
```

## Output

Running BABMiner creates:

```
results/
├── fragments_pdb/          # PDB structure files
│   ├── 1abc_1.pdb
│   └── ...
├── fragments_fasta/        # Amino acid sequences
│   ├── 1abc_1.fasta
│   └── ...
└── metadata.csv            # Motif metadata
```

### CSV Columns

| Column | Description |
|--------|-------------|
| `identifier` | Motif ID (e.g., `1abc_1`) |
| `chain` | Chain identifier |
| `motif_type` | Network distance |
| `beta1`, `beta2` | Residue ranges |
| `orientation` | parallel/antiparallel |
| `handedness` | right/left/unknown |
| `coverage` | Sequence fraction covered |

## Python API

```python
from BABMiner import bab_finder
from BABMiner.utils import load_pdb_structure

# Load structure
structure, model = load_pdb_structure('protein.pdb')
from Bio.PDB import DSSP
dssp = DSSP(model, 'protein.pdb')

# Find motifs
motifs, coverage = bab_finder(dssp, model, angle_threshold=40, min_bonds=2)

for motif in motifs:
    print(f"Chain {motif['chain']}: {motif['orientation']}")
```

## Algorithm

1. Load PDB and compute DSSP secondary structure
2. Identify contiguous β-strands (≥3 residues)
3. Calculate hydrogen bonds and angles between strands
4. Build graph of connected strands
5. Find β-α-β sequence patterns
6. Filter and compute orientation/handedness

## License

GNU General Public License v3 (see LICENSE file)

## Contributing

Contributions welcome! Please fork, add features, and submit pull requests.

## Citation

If you use BABMiner, please cite:

```bibtex
@software{babminer2026,
  title={BABMiner: β-α-β protein motif tool},
  author={Andre Lecona Buttelli},
  year={2026},
  url={https://github.com/yourusername/BABMiner}
}
```
