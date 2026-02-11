"""
Calculates features for peptide sequences using Macrel library.
"""
import pandas as pd
import numpy as np
from macrel.macrel_features import compute_all
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import torch
from transformers import AutoTokenizer, AutoModelForTokenClassification, pipeline
import math
import peptides
import loader


def macrel_descriptors_from_seq(sequence: str) -> np.ndarray:
    """
    Extracts 22 features using the compute_all function from macrel.
    """
    features = compute_all(sequence)

    # Convert to float32 numpy array for ONNX compatibility
    return np.array(features, dtype=np.float32).reshape(1, -1)


def peptides_descriptors_from_seqs(seqs: list[str]):
    """
    Extracts about 50 descriptors from package peptides.
    """
    descriptors_list_list = []
    for seq in seqs:
        descriptors = peptides.Peptide(seq).descriptors()  # compute descriptors
        # add row
        descriptors_list_list.append(descriptors)
    descriptors_df = pd.DataFrame(descriptors_list_list)
    return descriptors_df

def peptides_descriptors_from_fasta(filepath):
    """
    Extracts about 50 descriptors from package peptides.
    """
    sequences = loader.load_fasta(filepath)
    descriptors_df = peptides_descriptors_from_seqs(sequences)
    return descriptors_df

def calculate_hydrophobic_moment(sequence):
    angle_step = math.radians(100)  # 100 degrees for alpha-helix
    sum_sin = 0
    sum_cos = 0

    # Eisenberg Hydrophobicity Scale
    EISENBERG = {
        'A': 0.62, 'R': -2.53, 'N': -0.78, 'D': -0.90, 'C': 0.29,
        'Q': -0.85, 'E': -0.74, 'G': 0.48, 'H': -0.40, 'I': 1.38,
        'L': 1.06, 'K': -1.50, 'M': 0.64, 'F': 1.19, 'P': 0.12,
        'S': -0.18, 'T': -0.05, 'W': 0.81, 'Y': 0.26, 'V': 1.08
    }

    for n, res in enumerate(sequence):
        if res in EISENBERG:
            h = EISENBERG[res]
            sum_sin += h * math.sin(n * angle_step)
            sum_cos += h * math.cos(n * angle_step)

    moment = math.sqrt(sum_sin ** 2 + sum_cos ** 2) / len(sequence)
    return moment


def count_helices(sequence: str, verbose=False) -> float:
    analysed_seq = ProteinAnalysis(sequence)

    # returns a tuple: (Helix, Turn, Sheet)
    secondary_structure = analysed_seq.secondary_structure_fraction()
    res = secondary_structure[0]

    if verbose:
        print(f"Alpha-helix fraction: {secondary_structure[0] * 100:.2f}%")
    return res



