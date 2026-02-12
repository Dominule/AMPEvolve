"""
This module takes a sequence with list of masked positions
and generates random possible completions from aabet (the alphabet of amino acids).
"""
import random
from modlamp.sequences import Helices
from modlamp.descriptors import GlobalDescriptor
from modlamp.sequences import Helices
from modlamp.descriptors import GlobalDescriptor

aabet = "ACDEFGHIKLMNPQRSTVWY"
aabet_without_C = "ADEFGHIKLMNPQRSTVWY"


# --- SETTINGS ---
used_alphabet = aabet
# ----------------

from modlamp.sequences import AmphipathicArc
from modlamp.descriptors import GlobalDescriptor, PeptideDescriptor

def generate_amphipatic_helices(n_candidates=100, min_len=12, max_len=25):
    """
    Generates potential AMPs, filtering for high amphipathicity and positive charge.
    """
    # 1. Generate Raw Candidates using an alpha-helix template (180 degree arc)
    # This class places Hydrophobic vs Polar AAs in a helical pattern.
    lib = AmphipathicArc(n_candidates, min_len, max_len)
    lib.generate_sequences(arcsize=180)
    raw_seqs = lib.sequences

    # 2. Calculate Hydrophobic Moment (uH) using PeptideDescriptor
    # We use the 'Eisenberg' scale as defined in your documentation.
    pep_desc = PeptideDescriptor(raw_seqs, scalename='Eisenberg')
    pep_desc.calculate_moment(window=1000, angle=100, modality='max')
    moments = pep_desc.descriptor # This is a numpy array

    # 3. Calculate Net Charge using GlobalDescriptor
    glob_desc = GlobalDescriptor(raw_seqs)
    glob_desc.calculate_charge(ph=7.4)
    charges = glob_desc.descriptor

    # 4. Filter for high amphipathicity and positive charge
    # Thresholds: Moment > 0.4 (Strongly amphipathic), Charge > 2.0 (Target bacteria)
    filtered_amps = []

    for i in range(len(raw_seqs)):
        uH = moments[i][0]
        charge = charges[i][0]

        if uH > 0.4 and charge > 2.0:
            filtered_amps.append(raw_seqs[i])

    return filtered_amps




def choose_aa(original_aa):
    possible_aas = [aa for aa in used_alphabet if aa != original_aa]
    return random.choice(possible_aas)


def generate_completions(sequence: str, masked_positions: list[int], num_completions: int) -> list[str]:
    """
    Generates completed peptides from source peptide.
    Changes only pointing positions (masked_positions).
    Args:
        sequence (str): source sequence
        masked_positions (list[int]):  Positions to be changed.
        num_completions (int): Number of generated peptides.
    """
    completions = []
    for _ in range(num_completions):
        completion = list(sequence)
        for pos in masked_positions:
            original_aa = completion[pos]
            new_aa = choose_aa(original_aa)
            completion[pos] = new_aa
        completions.append("".join(completion))
    return completions