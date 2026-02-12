"""
Calculates features for peptide sequences using Macrel library.
"""
import pandas as pd
import numpy as np
from macrel.macrel_features import compute_all
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import peptides
import loader
from modlamp.descriptors import GlobalDescriptor, PeptideDescriptor


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


def alphahelices(sequence: str, verbose=False) -> float:
    analysed_seq = ProteinAnalysis(sequence)

    # returns a tuple: (Helix, Turn, Sheet)
    secondary_structure = analysed_seq.secondary_structure_fraction()
    res = secondary_structure[0]

    if verbose:
        print(f"Alpha-helix fraction: {secondary_structure[0] * 100:.2f}%")
    return res

def hydrophobic_moment(seq: str) -> float:
    # global uH value
    calc = PeptideDescriptor(seq, 'eisenberg')
    calc.calculate_moment(window=1000, angle=100, modality='mean')
    return float(calc.descriptor[0][0])

def hydrophobicity(seq: str) -> float:
    # average hydrophobicity (H)
    calc = PeptideDescriptor(seq, 'eisenberg')
    calc.calculate_global(window=1000, modality='mean')
    return float(calc.descriptor[0][0])

def charge(seq, amide=False) -> float:
    calc = GlobalDescriptor(seq)
    calc.calculate_charge(amide=amide)
    return float(calc.descriptor[0][0])

def amps_analysis(seqs: list[str], verbose=False) -> pd.DataFrame:
    hm = []
    h = []
    c = []
    ah = []
    for seq in seqs:
        hm.append(hydrophobic_moment(seq))
        h.append(hydrophobicity(seq))
        c.append(charge(seq))
        ah.append(alphahelices(seq))
    analyzed = {
        "Sequence" : seqs,
        "Hydrophobic moment" : hm,
        "Hydrophobicity" : h,
        "Charge" : c,
        "Alphahelices" : ah
    }
    if verbose:
        print(pd.DataFrame(analyzed))

    return pd.DataFrame(analyzed)




