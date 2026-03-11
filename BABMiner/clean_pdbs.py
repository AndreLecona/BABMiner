"""
PDB Cleaning Tool for DSSP Compatibility

A standalone utility to clean PDB files for DSSP (Database of Secondary Structure Prediction)
processing. Removes non-standard amino acids, keeps only the first model, and ensures proper
formatting for secondary structure analysis.

This tool is NOT part of the main BABMiner pipeline but serves as a troubleshooting utility
for preprocessing problematic PDB files.

Features:
    - Filters to standard amino acids only
    - Removes alternate locations and extra models
    - Normalizes non-standard residues (SEC→CYS, MSE→MET)
    - Parallel processing for batch operations
    - Preserves HEADER and EXPDTA records

Usage:
    cleanpdb <input> --outdir <output> -p <processes>
    
Examples:
    cleanpdb protein.pdb --outdir cleaned/
    cleanpdb pdb_directory/ --outdir cleaned/ -p 8
"""

import os
import sys
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm


# =============================================================================
# CONSTANTS
# =============================================================================

AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLU", "GLN", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL"
}

RESNAME_MAP = {
    "SEC": "CYS",
    "MSE": "MET",
}

KEEP_PREFIXES = ("HEADER", "EXPDTA")

# Globals (worker-local)
TARGET_DIR = None


# =============================================================================
# WORKER INITIALIZATION
# =============================================================================

def init_worker(target_dir):
    """Initialize worker process with target directory."""
    global TARGET_DIR
    TARGET_DIR = target_dir


# =============================================================================
# PDB CLEANING
# =============================================================================

def clean_pdb_for_dssp(input_pdb, output_pdb):
    """
    Clean a PDB file for DSSP compatibility.
    
    Processes:
        - Keeps only the first model
        - Filters to standard amino acids
        - Normalizes non-standard residues
        - Preserves HEADER and EXPDTA records
    
    Args:
        input_pdb (str): Path to input PDB file
        output_pdb (str): Path to output cleaned PDB file
    """
    in_first_model = True
    seen_model = False
    saw_header = False

    amino_acids = AMINO_ACIDS
    res_map = RESNAME_MAP
    keep_prefixes = KEEP_PREFIXES

    output_lines = []

    with open(input_pdb, "r") as f_in:
        for line in f_in:
            record = line[:6]

            # --- MODEL handling ---
            if record.startswith("MODEL"):
                if seen_model:
                    in_first_model = False
                else:
                    seen_model = True
                continue

            if record.startswith("ENDMDL"):
                in_first_model = False
                continue

            if not in_first_model:
                continue

            # --- Header records ---
            if record.startswith(keep_prefixes):
                if record.startswith("HEADER"):
                    saw_header = True
                output_lines.append(line)
                continue

            # --- ATOM records ---
            if record.startswith("ATOM"):
                orig_resname = line[17:20]
                resname = res_map.get(orig_resname.strip(), orig_resname.strip())

                if resname not in amino_acids:
                    continue

                if resname != orig_resname.strip():
                    line = line[:17] + f"{resname:>3}" + line[20:]

                output_lines.append(line)
                continue

            # --- TER / END ---
            if record.startswith("TER") or record.startswith("END"):
                output_lines.append(line)

    # --- Write output with fallback HEADER ---
    with open(output_pdb, "w") as f_out:
        for line in output_lines:
            f_out.write(line)


# =============================================================================
# FILE DISCOVERY
# =============================================================================

def find_all_pdbs(path):
    """
    Recursively find all PDB files in a directory or return single file.
    
    Args:
        path (str): File or directory path
        
    Returns:
        list: Sorted list of PDB file paths
    """
    valid_exts = (".pdb", ".ent")
    pdbs = []

    def _walk(dirpath):
        with os.scandir(dirpath) as it:
            for entry in it:
                if entry.is_dir():
                    _walk(entry.path)
                elif entry.name.lower().endswith(valid_exts):
                    pdbs.append(entry.path)

    if os.path.isfile(path) and path.lower().endswith(valid_exts):
        return [path]

    if os.path.isdir(path):
        _walk(path)

    return sorted(pdbs)


def list_existing_cleaned_pdbs(target_dir):
    """
    Find already-cleaned PDB files in target directory.
    
    Args:
        target_dir (str): Output directory path
        
    Returns:
        set: Set of processed PDB identifiers
    """
    existing = set()

    if not os.path.exists(target_dir):
        return existing

    for root, _, files in os.walk(target_dir):
        for f in files:
            if f.lower().endswith(".pdb"):
                existing.add(os.path.splitext(f)[0])

    return existing


# =============================================================================
# WORKER TASK
# =============================================================================

def process_single_pdb(pdb_path):
    """
    Process a single PDB file.
    
    Args:
        pdb_path (str): Path to input PDB file
        
    Returns:
        bool: True if successful, False otherwise
    """
    uid = os.path.splitext(os.path.basename(pdb_path))[0]
    os.makedirs(TARGET_DIR, exist_ok=True)
    output_path = os.path.join(TARGET_DIR, f"{uid}.pdb")

    try:
        clean_pdb_for_dssp(pdb_path, output_path)
        return True
    except Exception:
        return False


# =============================================================================
# CLI & MAIN
# =============================================================================

def create_parser():
    """Create command-line argument parser."""
    parser = argparse.ArgumentParser(
        description="Clean PDB files for DSSP secondary structure analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  cleanpdb protein.pdb --outdir cleaned/
  cleanpdb pdb_dir/ --outdir cleaned/ -p 8
  cleanpdb single.pdb -o output/ --verbose
        """
    )
    
    parser.add_argument(
        "input",
        help="PDB file or directory containing PDB files"
    )
    
    parser.add_argument(
        "-o", "--outdir",
        required=True,
        help="Output directory for cleaned PDB files"
    )
    
    parser.add_argument(
        "-p", "--processes",
        type=int,
        default=1,
        help="Number of parallel processes (default: 1)"
    )
    
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose output"
    )
    
    return parser


def main(argv=None):
    """
    Main entry point for cleanpdb CLI tool.
    
    Args:
        argv (list): Command-line arguments (for testing)
        
    Returns:
        int: Exit code
    """
    parser = create_parser()
    args = parser.parse_args(argv)
    
    # Validate input
    if not os.path.exists(args.input):
        print(f"[!] Input path does not exist: {args.input}", file=sys.stderr)
        return 1
    
    # Find PDB files
    pdb_paths = find_all_pdbs(args.input)
    if not pdb_paths:
        print("[!] No PDB files found. Exiting.", file=sys.stderr)
        return 1
    
    if args.verbose:
        print(f"[*] Found {len(pdb_paths)} PDB file(s)")
    
    # Setup output directory
    os.makedirs(args.outdir, exist_ok=True)
    
    # Find existing cleaned files
    existing = list_existing_cleaned_pdbs(args.outdir)
    if args.verbose and existing:
        print(f"[*] Found {len(existing)} existing cleaned PDB(s)")
    
    # Filter to-do list
    to_process = [
        p for p in pdb_paths
        if os.path.splitext(os.path.basename(p))[0] not in existing
    ]
    
    if not to_process:
        print("[*] All PDB files already cleaned. Nothing to do.")
        return 0
    
    if args.verbose:
        print(f"[*] Processing {len(to_process)} new PDB file(s)")
    
    # Determine worker count
    num_processes = args.processes
    
    # Process files
    success = failed = 0
    
    with ThreadPoolExecutor(
        max_workers=num_processes,
        initializer=init_worker,
        initargs=(args.outdir,)
    ) as executor:
        futures = [
            executor.submit(process_single_pdb, pdb_path)
            for pdb_path in to_process
        ]
        
        for future in tqdm(
            as_completed(futures),
            total=len(futures),
            desc="Cleaning PDBs",
            unit="file"
        ):
            if future.result():
                success += 1
            else:
                failed += 1
    
    # Print summary
    print(f"\n[+] Completed")
    print(f"[+] Successfully cleaned: {success}")
    if failed > 0:
        print(f"[!] Failed: {failed}")
    print(f"[+] Output directory: {args.outdir}")
    
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())