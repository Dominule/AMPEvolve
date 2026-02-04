"""
This module takes a sequence with list of masked positions
and generates random possible completions from aabet (the alphabet of amino acids).
"""
import random

aabet = "ACDEFGHIKLMNPQRSTVWY"
aabet_without_C = "ADEFGHIKLMNPQRSTVWY"


# --- SETTINGS ---
used_alphabet = aabet
# ----------------


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