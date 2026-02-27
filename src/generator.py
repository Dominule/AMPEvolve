"""
This module generates various thingies...
"""
import random
from modlamp.sequences import Helices
from modlamp.descriptors import GlobalDescriptor
from modlamp.sequences import Helices
from modlamp.descriptors import GlobalDescriptor

aabet = "ACDEFGHIKLMNPQRSTVWY"
aabet_without_C = "ADEFGHIKLMNPQRSTVWY"


from modlamp.sequences import AmphipathicArc
from modlamp.descriptors import GlobalDescriptor, PeptideDescriptor
from modlamp.sequences import Kinked
from modlamp.descriptors import GlobalDescriptor, PeptideDescriptor


def generate_killer_kinks(n_candidates=10, min_len=18, max_len=25):
    """
    Generates kinked AMPs and identifies the hinge (Proline) position.
    """
    # 1. Generate sequences using the Kinked rule (basic residues every 3-4 AAs)
    lib = Kinked(n_candidates, min_len, max_len)
    lib.generate_sequences()
    raw_seqs = lib.sequences

    # 2. Compute descriptors for filtering
    pep = PeptideDescriptor(raw_seqs, 'eisenberg')
    pep.calculate_moment()

    glob = GlobalDescriptor(raw_seqs)
    glob.calculate_charge(ph=7.4)

    final_candidates = []

    for i, seq in enumerate(raw_seqs):
        uH = pep.descriptor[i][0]
        z = glob.descriptor[i][0]

        # Find the Proline (kink) index
        pro_pos = seq.find('P') + 1  # 1-based indexing for researchers

        # FILTER:
        # High moment (>0.4), Good charge (>3), and Centralized Kink (between 30%-70% of length)
        rel_pos = pro_pos / len(seq)
        if uH > 0.4 and z > 3.0 and (0.3 < rel_pos < 0.7):
            final_candidates.append({
                'sequence': seq,
                'moment': round(uH, 3),
                'charge': round(z, 2),
                'kink_at': pro_pos
            })

    return [candidate["sequence"] for candidate in final_candidates]


def generate_amphipatic_helices(n_candidates=100, min_len=18, max_len=25):
    """
    Generates potential helical AMPs, filtering for high amphipathicity and positive charge.
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




def choose_aa(original_aa, used_alphabet=aabet):
    possible_aas = [aa for aa in used_alphabet if aa != original_aa]
    return random.choice(possible_aas)


def generate_completions(sequence: str, masked_positions: list[int], num_completions: int, forget_cys=False) -> list[str]:
    """
    Generates completed peptides from source peptide.
    Changes only pointing positions (masked_positions).
    Args:
        sequence (str): source sequence
        masked_positions (list[int]):  Positions to be changed.
        num_completions (int): Number of generated peptides.
    """
    if forget_cys:
        used_alphabet = aabet_without_C
    else:
        used_alphabet = aabet

    completions = []
    for _ in range(num_completions):
        completion = list(sequence)
        for pos in masked_positions:
            original_aa = completion[pos]
            new_aa = choose_aa(original_aa, used_alphabet=used_alphabet)
            completion[pos] = new_aa
        completions.append("".join(completion))
    return completions


def generate_neighbours(sequence: str, num_neighbours: int, forget_cys=False) -> list[str]:
    """
    Generates peptides that differ by one amino acid from the source sequence.
    Args:
        sequence (str): source sequence
        num_neighbours (int): Number of generated peptides.
    """
    if forget_cys:
        used_alphabet = aabet_without_C
    else:
        used_alphabet = aabet

    neighbours = []
    for _ in range(num_neighbours):
        neighbour = list(sequence)
        pos = random.randint(0, len(sequence) - 1)
        original_aa = neighbour[pos]
        new_aa = choose_aa(original_aa, used_alphabet=used_alphabet)
        neighbour[pos] = new_aa
        neighbours.append("".join(neighbour))
    return neighbours


def generate_all_neighbors(peptide: str, alphabet: str = aabet_without_C) -> list[str]:
    neighbors = []
    for i in range(len(peptide)):
        neighbors.append(peptide[:i] + peptide[i + 1:])
        for aa in alphabet:
            if aa != peptide[i]:
                neighbors.append(peptide[:i] + aa + peptide[i + 1:])
                neighbors.append(peptide[:i] + aa + peptide[i:])
    return neighbors
