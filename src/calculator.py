"""
Calculates features for peptide sequences using Macrel library.
"""
import pandas as pd
import numpy as np
from charset_normalizer.md import lru_cache
from macrel.macrel_features import compute_all
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import peptides
import loader
from modlamp.descriptors import GlobalDescriptor, PeptideDescriptor

from predictor import Predictor


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


@lru_cache(maxsize=10000)
def hydrophobic_moment(seq: str) -> float:
    # global uH value
    calc = PeptideDescriptor(seq, 'eisenberg')
    calc.calculate_moment(window=1000, angle=100, modality='mean')
    return float(calc.descriptor[0][0])


@lru_cache(maxsize=10000)
def hydrophobicity(seq: str) -> float:
    # average hydrophobicity (H)
    calc = PeptideDescriptor(seq, 'eisenberg')
    calc.calculate_global(window=1000, modality='mean')
    return float(calc.descriptor[0][0])


@lru_cache(maxsize=10000)
def charge(seq, amide=False) -> float:
    calc = GlobalDescriptor(seq)
    calc.calculate_charge(amide=amide)
    return float(calc.descriptor[0][0])


def amps_analysis(seqs: list[str], verbose=False) -> pd.DataFrame:
    hm = []     # hydrophobic moment
    h = []      # hydrophobicity
    c = []      # charge
    ah = []     # alphahelices
    for seq in seqs:
        hm.append(round(hydrophobic_moment(seq), 2))
        h.append(round(hydrophobicity(seq), 2))
        c.append(round(charge(seq), 2))
        ah.append(round(alphahelices(seq), 2))
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


class AMPKillerPredictor(Predictor):
    def __init__(self):
        import predictor
        self.macrel = predictor.MacrelPredictor()

    def calculate_and_predict_seqs(self, sequences: list[str]) -> list[float]:
        return [self.amp_killer_score(seq) for seq in sequences]

    def amp_killer_score(self, seq: str) -> float:
        """
        Calculates composite fitness based on Macrel and physical features.
        Bonuses for:
        - Moment (0.4-0.59),
        - Ratio (0.5-0.7),
        - Charge (+5 to +7),
        """

        # Ensure this returns the probability (0.0 to 1.0)
        macrel_prob = self.macrel.calculate_and_predict_seqs([seq])[0]

        # 2. Physical Descriptors (The nudges)
        pep = PeptideDescriptor(seq, 'eisenberg')
        pep.calculate_moment(window=1000, angle=100)
        uH = pep.descriptor[0][0]  # moment
        glob = GlobalDescriptor(seq)
        glob.calculate_charge(ph=7.4)
        z = glob.descriptor[0][0]  # charge
        glob.hydrophobic_ratio()
        h_ratio = glob.descriptor[0][0]  # ratio

        # --- WEIGHTING LOGIC ---
        # Multiplier of 2.0 ensures Macrel dominates the physical bonuses
        score = macrel_prob * 2

        # 1. Hydrophobic Moment (0.4 - 0.59)
        if 0.4 <= uH <= 0.59:
            score += 1
        elif uH > 0.59:
            score -= 1

        # 2. Hydrophobic Ratio (0.5 - 0.7)
        if 0.5 <= h_ratio <= 0.7:
            score += 1

        # 3. Net Charge (+5 to +7)
        if 5 <= z <= 7:
            score += 1
        #
        # # 4. Proline Centering
        # pro_pos = seq.find('P')
        # if pro_pos != -1:
        #     rel_pos = (pro_pos + 1) / len(seq)
        #     centering = 1.0 - abs(0.5 - rel_pos) * 2
        #     score += (centering * 0.2)
        #
        # # 5. R over K Preference
        # r_count = seq.count('R')
        # k_count = seq.count('K')
        # score += (r_count * 0.02) - (k_count * 0.05)

        return score
