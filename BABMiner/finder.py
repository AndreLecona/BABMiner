"""
Main BAB motif finding logic.
"""

import networkx as nx
from .utils import (
    find_beta_segments, map_residues_to_segments, calculate_beta_vectors, 
    get_beta_range, get_ca_coordinate, compute_bond_info, build_beta_sheet_network,
    find_bab_candidates, extract_path_info, extract_sequence_from_dssp,
    determine_orientation, calculate_signed_angle, determine_handedness
)
from .constants import HELIX_CODE


def find_coverage(true_bab, dssp):
    """
    Calculate fraction of protein sequence covered by BAB motifs.
    
    Args:
        true_bab: List of BAB motif dictionaries
        dssp: DSSP output dictionary
        
    Returns:
        float: Fraction of sequence covered (0-1)
    """
    keys = sorted(dssp.keys())
    protein_length = len(keys)
    
    covered_indices = set()
    for bab in true_bab:
        start = bab['pos_start']
        end = bab['pos_end']
        for i in range(start, end + 1):
            if 0 <= i < protein_length:
                covered_indices.add(i)
    
    return len(covered_indices) / protein_length if protein_length > 0 else 0.0


def filter_true_babs(bab_list, G, beta_segments, bond_info, model, dssp):
    """
    Filter BAB candidates to only those with valid network paths.
    Compute handedness and orientation information.
    
    Args:
        bab_list: Candidates from find_bab_candidates()
        G: Beta sheet network graph
        beta_segments: List of beta segments
        bond_info: Bond information dictionary
        model: BioPython model object
        dssp: DSSP output dictionary
        
    Returns:
        list: Filtered list of BAB motif dictionaries
    """
    true_bab = []
    
    for seg_start, seg_end, pos_start, pos_end in bab_list:
        if not nx.has_path(G, seg_start, seg_end):
            continue
        
        path, distance, angles, bonds = extract_path_info(G, bond_info, seg_start, seg_end)
        chain = beta_segments[seg_start][0][0]
        beta1_range = get_beta_range(beta_segments, seg_start)
        beta2_range = get_beta_range(beta_segments, seg_end)
        segment_length = pos_end - pos_start - 1
        
        orientation = determine_orientation(angles)
        bab_sequence = extract_sequence_from_dssp(dssp, pos_start, pos_end)
        orientation_str = "parallel" if orientation == 0 else "antiparallel"
        
        handedness = None
        angle = None
        
        if orientation == 0:  # Only compute handedness for parallel orientations
            keys = sorted(dssp.keys())
            ca_beta1_start = get_ca_coordinate(model, chain, beta_segments[seg_start][0])
            ca_beta1_end = get_ca_coordinate(model, chain, beta_segments[seg_start][-1])
            ca_beta2_end = get_ca_coordinate(model, chain, beta_segments[seg_end][-1])
            
            helix_keys = []
            for i in range(pos_start + 1, pos_end):
                if keys[i] in dssp and dssp[keys[i]][2] == HELIX_CODE:
                    helix_keys.append(keys[i])
            
            if helix_keys:
                ca_list_coords = []
                for key in helix_keys:
                    ca = get_ca_coordinate(model, chain, key)
                    if ca is not None:
                        ca_list_coords.append(ca)
                
                if (ca_list_coords and all([ca_beta1_start is not None, 
                                           ca_beta1_end is not None, 
                                           ca_beta2_end is not None])):
                    angle_list = calculate_signed_angle(ca_beta1_start, ca_beta1_end, 
                                                       ca_beta2_end, ca_list_coords)
                    handedness = determine_handedness(angle_list)
                    angle = sum(angle_list) / len(angle_list) if angle_list else None
        
        true_bab.append({
            'chain': chain,
            'distance': distance,
            'beta1_range': beta1_range,
            'beta2_range': beta2_range,
            'segment_length': segment_length,
            'bonds': bonds,
            'angles': angles,
            'orientation': orientation_str,
            'handedness': handedness,
            'angle': angle,
            'seg_start': seg_start,
            'seg_end': seg_end,
            'pos_start': pos_start,
            'pos_end': pos_end,
            'sequence': bab_sequence,
            'dssp': ''.join([dssp[sorted(dssp.keys())[i]][2] for i in range(pos_start, pos_end + 1)])
        })
    
    return true_bab


def bab_finder(dssp, model, angle_threshold=40, min_bonds=2):
    """
    Main function: Analyze structure and extract all BAB motifs.
    
    Args:
        dssp: DSSP secondary structure assignment
        model: BioPython model object
        angle_threshold: Maximum angle threshold for parallel/antiparallel (degrees)
        min_bonds: Minimum hydrogen bonds to connect segments
        
    Returns:
        tuple: (list of BAB motif dicts, sequence coverage fraction)
    """
    beta_segments = find_beta_segments(dssp)
    if len(beta_segments) < 2:
        return [], 0.0
    
    res_to_seg = map_residues_to_segments(beta_segments)
    vectors = calculate_beta_vectors(model, beta_segments)
    bond_info = compute_bond_info(dssp, beta_segments, vectors)
    
    G = build_beta_sheet_network(beta_segments, bond_info, 
                                angle_threshold=angle_threshold, 
                                min_bonds=min_bonds)
    
    ss_list = [dssp[key][2] for key in sorted(dssp.keys())]
    bab_candidates = find_bab_candidates(ss_list, dssp, res_to_seg)
    true_bab = filter_true_babs(bab_candidates, G, beta_segments, bond_info, model, dssp)
    coverage = find_coverage(true_bab, dssp)
    
    return true_bab, coverage
