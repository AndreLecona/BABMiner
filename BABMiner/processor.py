"""
Per-file analysis and fragment extraction.
"""

import os
from Bio.PDB import Dice, DSSP
from .utils import load_pdb_structure
from .finder import bab_finder
from .io import get_pdb_code_from_file, get_exp_method_from_file


def analyze_pdb_file(pdb_file, dir_pdb, dir_fa, angle_threshold=40, min_hb=2):
    """
    Analyze single PDB file and extract BAB fragments.
    
    Creates fragment FASTA and PDB files in output directories.
    
    Args:
        pdb_file: Path to PDB file
        dir_pdb: Output directory for PDB fragments
        dir_fa: Output directory for FASTA fragments
        angle_threshold: Angle threshold for motif detection
        min_hb: Minimum hydrogen bonds
        
    Returns:
        list: List of metadata dictionaries for found motifs
    """
    uid = os.path.basename(pdb_file)
    uid_noext = uid[:-4] if uid.lower().endswith(".pdb") else uid
    
    # Load structure and compute DSSP
    try:
        structure, model = load_pdb_structure(pdb_file)
        dssp = DSSP(model, pdb_file)
    except Exception as e:
        print(f"[!] Error processing {pdb_file}: {e}")
        return []
    
    # Find motifs
    babs, coverage = bab_finder(dssp, model, angle_threshold=angle_threshold, min_bonds=min_hb)
    if not babs:
        return []
    
    # Extract fragments and metadata
    results = []
    pdb_code = get_pdb_code_from_file(pdb_file).lower()
    
    for idx, bab in enumerate(babs, 1):
        identifier = f"{pdb_code}_{idx}"
        
        # Fragment range
        beta1_start_res = bab['beta1_range'][0]
        beta2_end_res = bab['beta2_range'][1]
        chain_id = bab['chain']
        
        # Save PDB fragment
        pdb_path = os.path.join(dir_pdb, f"{identifier}.pdb")
        try:
            Dice.extract(
                structure,
                chain_id=chain_id,
                start=beta1_start_res,
                end=beta2_end_res,
                filename=pdb_path
            )
        except Exception as e:
            print(f"[!] Could not extract PDB fragment {identifier}: {e}")
        
        # Save FASTA
        fa_path = os.path.join(dir_fa, f"{identifier}.fasta")
        try:
            with open(fa_path, "w") as fh:
                fh.write(f">{identifier}\n")
                fh.write(bab['sequence'] + "\n")
        except Exception as e:
            print(f"[!] Could not write FASTA {identifier}: {e}")
        
        # Build metadata
        meta = {
            "identifier": identifier,
            "uid": uid_noext,
            "chain": chain_id,
            "motif_type": bab['distance'],
            "beta1": f"{beta1_start_res}-{bab['beta1_range'][1]}",
            "beta2": f"{bab['beta2_range'][0]}-{beta2_end_res}",
            "h&l_len": bab['segment_length'],
            "length": beta2_end_res - beta1_start_res + 1,
            "res_bonded": f"[{';'.join(str(x) for x in bab['bonds'])}]",
            "angles": f"[{';'.join(f'{x:.2f}' for x in bab['angles'])}]",
            "orientation": bab['orientation'],
            "handedness": bab['handedness'] or "",
            "angle": f"{bab['angle']:.2f}" if bab['angle'] is not None else "",
            "dssp": bab['dssp'],
            "exp_method": get_exp_method_from_file(pdb_file),
            "coverage": f"{coverage:.3f}"
        }
        
        results.append(meta)
    
    return results
