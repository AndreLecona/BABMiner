"""
BABMiner: β-α-β Motif Mining Tool
A CLI tool for identifying and extracting beta-alpha-beta protein motifs from PDB files.
"""

__version__ = "1.0.0"
__author__ = "[YOUR_NAME]"

from .finder import bab_finder

__all__ = ["bab_finder"]
