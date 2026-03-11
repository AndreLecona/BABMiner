"""
Command-line interface for BAB motif finder.
"""

import os
import sys
import argparse
import traceback
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm

from .constants import CSV_COLUMNS, DEFAULT_ANGLE_THRESHOLD, DEFAULT_MIN_BONDS
from .io import gather_pdb_files, mkdirs_for_outdir, write_csv
from .processor import analyze_pdb_file


def create_parser():
    """Create and return argument parser."""
    parser = argparse.ArgumentParser(
        description="β-α-β (beta-alpha-beta) motif finder - Identify and extract BAB motifs from PDB files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single file
  python -m find_babs path/to/protein.pdb --outdir results
  
  # Directory with parallelization
  python -m find_babs /path/to/pdbs --outdir results --processes 8
  
  # Custom parameters
  python -m find_babs data/ --outdir output --angle 45 --min_hb 3 -p 4
        """
    )
    
    parser.add_argument(
        "input",
        help="PDB file or directory containing PDB files"
    )
    
    parser.add_argument(
        "-o", "--outdir",
        required=True,
        help="Output directory (creates: fragments_pdb/, fragments_fasta/, metadata.csv)"
    )
    
    parser.add_argument(
        "-p", "--processes",
        type=int,
        default=None,
        help="Number of parallel processes (default: auto-detect CPU count)"
    )
    
    parser.add_argument(
        "--angle",
        type=float,
        default=DEFAULT_ANGLE_THRESHOLD,
        help=f"Maximum angle threshold for strand orientation in degrees (default: {DEFAULT_ANGLE_THRESHOLD})"
    )
    
    parser.add_argument(
        "--min_hb",
        type=int,
        default=DEFAULT_MIN_BONDS,
        help=f"Minimum hydrogen bonds between strands (default: {DEFAULT_MIN_BONDS})"
    )
    
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose output"
    )
    
    return parser


def main(argv=None):
    """Main CLI entry point."""
    parser = create_parser()
    args = parser.parse_args(argv)
    
    # Validate input
    if not os.path.exists(args.input):
        print(f"[!] Input path does not exist: {args.input}", file=sys.stderr)
        return 1
    
    # Gather PDB files
    pdb_files = gather_pdb_files(args.input)
    if not pdb_files:
        print("[!] No PDB files found. Exiting.", file=sys.stderr)
        return 1
    
    if args.verbose:
        print(f"[*] Found {len(pdb_files)} PDB file(s)")
    
    # Setup output directories
    os.makedirs(args.outdir, exist_ok=True)
    dir_pdb, dir_fa = mkdirs_for_outdir(args.outdir)
    
    if args.verbose:
        print(f"[*] Output directory: {args.outdir}")
    
    # Determine number of processes
    num_processes = args.processes
    if num_processes is None:
        num_processes = min(32, (os.cpu_count() or 4) * 2)
    
    if args.verbose:
        print(f"[*] Using {num_processes} process(es)")
    
    # Process files
    all_results = []
    
    if num_processes > 1 and len(pdb_files) > 1:
        # Parallel processing
        with ThreadPoolExecutor(max_workers=num_processes) as exe:
            futures = {
                exe.submit(analyze_pdb_file, pdb, dir_pdb, dir_fa, args.angle, args.min_hb): pdb 
                for pdb in pdb_files
            }
            
            with tqdm(total=len(futures), desc="Processing PDBs", unit="pdb") as pbar:
                for fut in as_completed(futures):
                    pdb = futures[fut]
                    try:
                        res = fut.result()
                        if res:
                            all_results.extend(res)
                    except Exception as e:
                        print(f"[!] Error processing {pdb}: {e}", file=sys.stderr)
                        if args.verbose:
                            traceback.print_exc()
                    pbar.update(1)
    else:
        # Sequential processing
        iterator = tqdm(pdb_files, desc="Processing PDBs", unit="pdb") if len(pdb_files) > 1 else pdb_files
        for pdb in iterator:
            try:
                res = analyze_pdb_file(pdb, dir_pdb, dir_fa, args.angle, args.min_hb)
                if res:
                    all_results.extend(res)
            except Exception as e:
                print(f"[!] Error processing {pdb}: {e}", file=sys.stderr)
                if args.verbose:
                    traceback.print_exc()
    
    # Write results
    csv_path = write_csv(args.outdir, all_results, CSV_COLUMNS)
    
    print(f"[+] Wrote metadata CSV: {csv_path}")
    print(f"[+] PDB fragments saved in: {dir_pdb}")
    print(f"[+] FASTA fragments saved in: {dir_fa}")
    print(f"[+] Total motifs found: {len(all_results)}")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
