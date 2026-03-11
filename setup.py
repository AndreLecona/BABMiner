"""
Setup configuration for BABMiner package.
"""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="BABMiner",
    version="1.0.0",
    author="Andre Lecona Buttelli",
    author_email="andrelecon@elsi.jp",
    description="A CLI tool for mining β-α-β protein motifs from PDB files",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/AndreLecona/BABMiner",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.7",
    install_requires=[
        "biopython>=1.79",
        "numpy>=1.19",
        "networkx>=2.5",
        "tqdm>=4.50",
    ],
    entry_points={
        "console_scripts": [
            "babminer=BABMiner.cli:main",
        ],
    },
)
