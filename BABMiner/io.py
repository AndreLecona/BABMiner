"""
File I/O operations for PDB and output files.
"""

import os
import csv
from glob import glob


def get_pdb_code_from_file(pdb_file):
    """Extract PDB code from file header or use filename."""
    filename_code = os.path.splitext(os.path.basename(pdb_file))[0]
    try:
        with open(pdb_file, "r") as fh:
            for line in fh:
                if line.startswith("HEADER"):
                    code = line[62:66].strip()
                    return code if code else filename_code
    except Exception:
        pass
    return filename_code


def get_exp_method_from_file(pdb_file):
    """Extract experimental method from PDB file header."""
    try:
        with open(pdb_file, "r") as fh:
            for line in fh:
                if line.startswith("EXPDTA"):
                    return line[6:].strip()
    except Exception:
        pass
    return "UNKNOWN"


def mkdirs_for_outdir(outdir):
    """Create output subdirectories for fragments."""
    p_pdb = os.path.join(outdir, "fragments_pdb")
    p_fa = os.path.join(outdir, "fragments_fasta")
    os.makedirs(p_pdb, exist_ok=True)
    os.makedirs(p_fa, exist_ok=True)
    return p_pdb, p_fa


def gather_pdb_files(input_path):
    """Gather all PDB files from input path (file or directory)."""
    if os.path.isfile(input_path):
        return [input_path]
    
    pdb_files = glob(os.path.join(input_path, "**", "*.pdb"), recursive=True)
    ent_files = glob(os.path.join(input_path, "**", "*.ent"), recursive=True)
    return sorted(pdb_files + ent_files)


def write_csv(outdir, rows, columns):
    """Write metadata results to CSV file."""
    csv_path = os.path.join(outdir, "metadata.csv")
    os.makedirs(outdir, exist_ok=True)
    
    with open(csv_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns)
        writer.writeheader()
        for r in rows:
            row = {k: r.get(k, "") for k in columns}
            writer.writerow(row)
    
    return csv_path
