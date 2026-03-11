# Project Structure & Module Guide

## Directory Layout

```
BABMiner/
├── find_babs/                  # Main package
│   ├── __init__.py            # Package initialization
│   ├── __main__.py            # Entry point for python -m find_babs
│   ├── cli.py                 # Command-line interface
│   ├── constants.py           # Constants, defaults, codes
│   ├── io.py                  # File I/O (PDB loading, CSV writing)
│   ├── utils.py               # Core algorithm functions (consolidated)
│   ├── finder.py              # Main motif finding orchestration
│   └── processor.py           # Per-file analysis workflow
├── setup.py                    # Package setup & installation
├── requirements.txt            # Python dependencies
├── README.md                   # Documentation
├── LICENSE                     # GNU GPLv3 License
├── .gitignore                  # Git ignore rules
└── STRUCTURE.md                # This file
```

## Module Responsibilities

| Module | Purpose |
|--------|---------|
| **`utils.py`** | Core algorithm functions: structure loading, beta segment detection, hydrogen bonding, network graph building, orientation/handedness analysis, BAB candidate identification |
| **`finder.py`** | Main orchestration: coordinates all analysis steps and returns final BAB motifs |
| **`processor.py`** | Per-file workflow: loads PDB, runs finder, extracts fragments (PDB/FASTA), assembles metadata |
| **`cli.py`** | Argument parsing and command-line interface with parallelization |
| **`io.py`** | File operations: gather PDB files, write CSV output, create directories |
| **`constants.py`** | Global constants, default parameters, secondary structure codes |

## Data Flow

```
User Input (CLI)
    ↓
gather_pdb_files() [io.py]
    ↓
For each PDB file (parallel or sequential):
    ├─ load_pdb_structure() [utils.py]
    ├─ DSSP computation
    ├─ bab_finder() [finder.py]
    │   ├─ find_beta_segments() [utils.py]
    │   ├─ calculate_beta_vectors() [utils.py]
    │   ├─ compute_bond_info() [utils.py]
    │   ├─ build_beta_sheet_network() [utils.py]
    │   ├─ find_bab_candidates() [utils.py]
    │   ├─ determine_orientation() [utils.py]
    │   └─ determine_handedness() [utils.py]
    └─ analyze_pdb_file() [processor.py]
        ├─ Dice.extract() → fragments_pdb/
        ├─ write FASTA → fragments_fasta/
        └─ collect metadata → results[]
    ↓
write_csv() [io.py] → metadata.csv
    ↓
Output
```

## Adding New Features

### To add a new analysis metric:
1. Implement calculation in `utils.py` (appropriate section)
2. Update `filter_true_babs()` in `finder.py` to compute and store
3. Add column to `CSV_COLUMNS` in `constants.py`
4. Update metadata assembly in `processor.py`

### To modify motif detection parameters:
1. Update defaults in `constants.py`
2. Add/modify CLI argument in `cli.py`
3. Pass through `processor.py` → `analyze_pdb_file()` → `bab_finder()`

### To change DSSP/structure handling:
1. Modify functions in `utils.py` (STRUCTURE HANDLING section)
2. Update error handling in `processor.py`

## Testing Checklist

- [ ] Single PDB file analysis
- [ ] Directory with multiple PDB files
- [ ] Parallel processing (-p flag)
- [ ] Custom angle threshold (--angle)
- [ ] Custom min H-bonds (--min_hb)
- [ ] CSV output format verification
- [ ] Fragment extraction (PDB & FASTA)
- [ ] Error handling for invalid files

## Development Workflow

```bash
# Clone and enter directory
git clone https://github.com/yourusername/BABMiner.git
cd BABMiner

# Editable install with dependencies
pip install -e .

# Test CLI
python -m find_babs --help
babminer --help
babminer test.pdb -o results
```

## Common Customizations

### Changing default parameters:
Edit `find_babs/constants.py`:
```python
DEFAULT_ANGLE_THRESHOLD = 45.0  # Change from 40
DEFAULT_MIN_BONDS = 3           # Change from 2
```

### Adding new DSSP-based filters:
Update logic in `finder.py` and `utils.py` (appropriate section)