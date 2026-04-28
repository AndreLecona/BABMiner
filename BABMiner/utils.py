"""
Unified utility functions for BAB motif finding.

This module consolidates functions from:
- Structure handling
- Beta segment detection and analysis
- Hydrogen bond analysis
- Network graph building
- Motif orientation and handedness
- BAB candidate identification and filtering
"""

from Bio.PDB import PDBParser
import numpy as np
import networkx as nx
from .constants import BETA_CODE, MIN_BETA_LENGTH, HELIX_CODE


# =============================================================================
# STRUCTURE HANDLING
# =============================================================================

def load_pdb_structure(pdb_file):
    """
    Load PDB file and return structure and model.
    
    Args:
        pdb_file: Path to PDB file
        
    Returns:
        tuple: (structure, model) or (None, None) on error
        
    Raises:
        Exception: If PDB parsing fails
    """
    parser = PDBParser()
    structure = parser.get_structure('pdb', pdb_file)
    model = structure[0]
    return structure, model


# =============================================================================
# BETA SEGMENT DETECTION & ANALYSIS
# =============================================================================

def find_beta_segments(dssp):
    """
    Identify contiguous beta strand regions from DSSP output.
    
    Args:
        dssp: DSSP output dictionary
        
    Returns:
        list: List of beta strand segments, each containing residue keys
    """
    beta_segments = []
    current = []
    
    for key in sorted(dssp.keys()):
        if dssp[key][2] == BETA_CODE:
            current.append(key)
        else:
            if len(current) >= MIN_BETA_LENGTH:
                beta_segments.append(current)
            current = []
    
    if len(current) >= MIN_BETA_LENGTH:
        beta_segments.append(current)
    
    return beta_segments


def map_residues_to_segments(beta_segments):
    """
    Create mapping from residue key to segment index.
    
    Args:
        beta_segments: List of beta segments
        
    Returns:
        dict: Mapping of residue key -> segment index
    """
    res_to_seg = {}
    for idx, seg in enumerate(beta_segments):
        for res in seg:
            res_to_seg[res] = idx
    return res_to_seg


def calculate_beta_vectors(model, beta_segments):
    """
    Compute direction vectors for each beta strand.
    
    Args:
        model: BioPython model object
        beta_segments: List of beta segments
        
    Returns:
        list: Direction vectors (numpy arrays) for each segment
    """
    vectors = []
    
    for seg in beta_segments:
        try:
            chain_id = seg[0][0]
            chain = model[chain_id]
            start_res = chain[seg[0][1]]
            end_res = chain[seg[-1][1]]
            
            if 'CA' in start_res and 'CA' in end_res:
                start_ca = start_res['CA'].get_coord()
                end_ca = end_res['CA'].get_coord()
                vec = end_ca - start_ca
                vectors.append(vec)
            else:
                vectors.append(np.array([0, 0, 0]))
        except Exception:
            vectors.append(np.array([0, 0, 0]))
    
    return vectors


def get_beta_range(beta_segments, seg_idx):
    """
    Extract residue number range for a beta segment.
    
    Args:
        beta_segments: List of beta segments
        seg_idx: Index of segment
        
    Returns:
        tuple: (start_res_num, end_res_num)
    """
    seg = beta_segments[seg_idx]
    start_res_num = seg[0][1][1]
    end_res_num = seg[-1][1][1]
    return start_res_num, end_res_num


def get_ca_coordinate(model, chain_id, res_key):
    """
    Get CA atom coordinate for a specific residue.
    
    Args:
        model: BioPython model object
        chain_id: Chain identifier
        res_key: Residue key tuple (chain_id, res_id)
        
    Returns:
        numpy.array: Coordinate array or None if not found
    """
    try:
        chain = model[chain_id]
        res = chain[res_key[1]]
        if 'CA' in res:
            return res['CA'].get_coord()
    except Exception:
        pass
    return None


# =============================================================================
# BETA SHEET NETWORK BUILDING
# =============================================================================

def build_beta_sheet_network(beta_segments, bond_info, angle_threshold=40, min_bonds=2):
    """
    Build graph of connected beta strands based on hydrogen bonding.
    
    Args:
        beta_segments: List of beta segments
        bond_info: Bond information dictionary from bonding.py
        angle_threshold: Maximum angle for parallel/antiparallel strands (degrees)
        min_bonds: Minimum hydrogen bonds to connect strands
        
    Returns:
        networkx.Graph: Network of connected beta segments
    """
    G = nx.Graph()
    G.add_nodes_from(range(len(beta_segments)))
    
    for i in range(len(beta_segments)):
        for j in range(i + 1, len(beta_segments)):
            count, angle, _ = bond_info[(i, j)]
            
            # Connect if angle is parallel (< threshold) or antiparallel (> 180-threshold)
            # AND sufficient hydrogen bonds
            if (angle < angle_threshold or abs(angle - 180) < angle_threshold) and count >= min_bonds:
                G.add_edge(i, j)
    
    return G


# =============================================================================
# MOTIF ANALYSIS (Orientation & Handedness)
# =============================================================================

def determine_orientation(angles):
    """
    Determine strand orientation (parallel or antiparallel).
    
    Args:
        angles: List of angles between beta strands
        
    Returns:
        int: 0 for parallel, 1 for antiparallel, None if no angles
    """
    if not angles:
        return None
    return sum(angle > 90 for angle in angles) % 2


def calculate_signed_angle(ca_beta1_start, ca_beta1_end, ca_beta2_end, ca_helix_list):
    """
    Calculate signed angles for handedness determination.
    
    Args:
        ca_beta1_start: CA coordinate of first beta strand start
        ca_beta1_end: CA coordinate of first beta strand end
        ca_beta2_end: CA coordinate of second beta strand end
        ca_helix_list: List of CA coordinates in helix region
        
    Returns:
        list: Signed angles in degrees
    """
    angles = []
    vec_b1 = ca_beta1_end - ca_beta1_start
    vec_ref = ca_beta2_end - ca_beta1_end
    vec_nom = np.cross(vec_b1, vec_ref)
    plane_normal = np.cross(vec_nom, vec_ref)
    
    for ca_helix in ca_helix_list:
        vec_helix = ca_helix - ca_beta1_end
        dot_pn_vh = np.dot(plane_normal, vec_helix)
        dot_pn_pn = np.dot(plane_normal, plane_normal)
        
        if dot_pn_pn > 0:
            proj = vec_helix - (dot_pn_vh / dot_pn_pn) * plane_normal
        else:
            proj = vec_helix
        
        cross_prod = np.cross(vec_ref, proj)
        dot_prod = np.dot(vec_ref, proj)
        signed_angle = np.degrees(np.arctan2(np.dot(plane_normal, cross_prod), dot_prod))
        angles.append(signed_angle)
    
    return angles


def determine_handedness(angle_list):
    """
    Determine handedness (right or left) based on weighted angles.
    
    Args:
        angle_list: List of signed angles
        
    Returns:
        str: "right", "left", "unknown", or None if empty
    """
    length = len(angle_list)
    if length == 0:
        return None
    
    center = (length - 1) / 2
    rights = 0.0
    lefts = 0.0
    
    for i, a in enumerate(angle_list):
        angle = a if a >= 0 else 360 + a
        # Larger weight near ends, smaller near center
        weight = abs(i - center) + 1
        
        if angle <= 180:
            rights += weight
        else:
            lefts += weight
    
    if rights > lefts:
        return "right"
    elif lefts > rights:
        return "left"
    else:
        return "unknown"


# =============================================================================
# HYDROGEN BOND ANALYSIS
# =============================================================================

def count_hydrogen_bonds(dssp, beta_segments, i, j):
    """
    Count hydrogen bonds between two beta segments.
    
    Args:
        dssp: DSSP output dictionary
        beta_segments: List of beta segments
        i: Index of first segment
        j: Index of second segment
        
    Returns:
        tuple: (bond_count, list of bonding residue pairs)
    """
    count = 0
    bonding_residues = []
    
    for res in beta_segments[i]:
        dssp_idx = dssp[res][0]
        for bp_idx in [6, 8]:  # DSSP indices for hydrogen bond partners
            bp = dssp[res][bp_idx]
            if bp != 0:
                target_idx = dssp_idx + bp
                for r2 in beta_segments[j]:
                    if dssp[r2][0] == target_idx:
                        count += 1
                        bonding_residues.append((res, r2))
    
    return count, bonding_residues


def calculate_angle(v1, v2):
    """
    Calculate angle between two vectors in degrees.
    
    Args:
        v1: First vector (numpy array)
        v2: Second vector (numpy array)
        
    Returns:
        float: Angle in degrees (0-180)
    """
    norm1 = np.linalg.norm(v1)
    norm2 = np.linalg.norm(v2)
    
    if norm1 > 0 and norm2 > 0:
        cos_angle = np.dot(v1, v2) / (norm1 * norm2)
        cos_angle = np.clip(cos_angle, -1, 1)
        return np.degrees(np.arccos(cos_angle))
    
    return 0


def compute_bond_info(dssp, beta_segments, vectors):
    """
    Calculate hydrogen bond counts and angles for all segment pairs.
    
    Args:
        dssp: DSSP output dictionary
        beta_segments: List of beta segments
        vectors: Beta strand direction vectors
        
    Returns:
        dict: Bond information keyed by (i, j) tuple with values
              (bond_count, angle, bonding_residues)
    """
    bond_info = {}
    
    for i in range(len(beta_segments)):
        for j in range(i + 1, len(beta_segments)):
            count, bonding_residues = count_hydrogen_bonds(dssp, beta_segments, i, j)
            angle = calculate_angle(vectors[i], vectors[j])
            bond_info[(i, j)] = (count, angle, bonding_residues)
    
    return bond_info


# =============================================================================
# BAB CANDIDATE IDENTIFICATION & FILTERING
# =============================================================================

def find_bab_candidates(ss_list, dssp, res_to_seg):
    """
    Identify all beta-alpha-beta sequence patterns.
    
    Args:
        ss_list: List of secondary structure codes in order
        dssp: DSSP output dictionary
        res_to_seg: Residue to segment mapping
        
    Returns:
        list: List of (seg_start, seg_end, pos_start, pos_end) tuples
    """
    bab_list = []
    keys = sorted(dssp.keys())
    
    for start in range(len(ss_list)):
        if ss_list[start] == BETA_CODE and keys[start] in res_to_seg:
            seg_start = res_to_seg[keys[start]]
            
            for end in range(start + 2, len(ss_list)):
                if (ss_list[end] == BETA_CODE and 
                    keys[end] in res_to_seg and 
                    keys[start][0] == keys[end][0]):  # Same chain
                    
                    seg_end = res_to_seg[keys[end]]
                    
                    if seg_start != seg_end:
                        middle = ss_list[start + 1:end]
                        # Alpha helix present and no other beta strands
                        if BETA_CODE not in middle and HELIX_CODE in middle:
                            bab_list.append((seg_start, seg_end, start, end))
    
    return bab_list


def extract_path_info(G, bond_info, seg_start, seg_end):
    """
    Get network path between segments and collect angle/bond info.
    
    Args:
        G: NetworkX graph of beta sheet
        bond_info: Bond information dictionary
        seg_start: Starting segment index
        seg_end: Ending segment index
        
    Returns:
        tuple: (path, distance, angles, bonds) where distance is path length - 1
    """
    path = nx.shortest_path(G, seg_start, seg_end)
    distance = len(path) - 1
    angles, bonds = [], []
    
    for i in range(len(path) - 1):
        idx_a, idx_b = path[i], path[i + 1]
        pair = (idx_a, idx_b) if idx_a < idx_b else (idx_b, idx_a)
        count, angle, _ = bond_info[pair]
        angles.append(angle)
        bonds.append(count)
    
    return path, distance, angles, bonds

def extract_sequence_from_dssp(dssp, pos_start, pos_end):
    """
    Extract amino acid sequence from DSSP positions.
    
    Args:
        dssp: DSSP output dictionary
        pos_start: Starting position in DSSP order
        pos_end: Ending position in DSSP order
        
    Returns:
        str: Single-letter amino acid sequence
    """
    keys = sorted(dssp.keys())
    res_map = {key[1][1]: dssp[key][1] for key in keys}

    seq = []
    for resnum in range(pos_start, pos_end + 1):
        aa = res_map.get(resnum, 'X')
        if aa in ('!', 'X'):
            aa = 'X'
        seq.append(aa)

    return ''.join(seq)
