"""
Constants and configuration for BAB motif finding.
"""

# Amino acid conversion mapping
THREE_TO_ONE = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F", 
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L", 
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R", 
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y"
}

# CSV output columns
CSV_COLUMNS = [
    "identifier",
    "uid",
    "chain",
    "motif_type",
    "beta1",
    "beta2",
    "h&l_len",
    "length",
    "res_bonded",
    "angles",
    "orientation",
    "handedness",
    "angle",
    "dssp",
    "exp_method",
    "coverage"
]

# Default parameters
DEFAULT_ANGLE_THRESHOLD = 40.0
DEFAULT_MIN_BONDS = 2
DEFAULT_PROCESSES = None  # Will use CPU count

# DSSP secondary structure codes
BETA_CODE = 'E'
HELIX_CODE = 'H'
MIN_BETA_LENGTH = 3
