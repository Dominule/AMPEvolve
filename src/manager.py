"""
Loop for mutants:
1) Take one seq
2) set possible positions for masking and how many masks for epoch # conservative = 1
3) generate and predict
4) predict for original seq and add to the list
5) keep n highest

Lead-optimization - Hill-climbing:
1) take one sequence
2) set possible positions for masking and how many masks for epoch # conservative = 1
3) generate and predict
4) predict for original seq and add to the list
5) keep only the highest


Generate and predict:
3A) load sequences
3B) calculate features
3C) predict with macrel
"""

import pandas as pd
import loader
from generator import generate_completions
import calculator
import plotter
from predictor import MacrelPredictor
from hill_climber import climb_high

ROOT_PATH = "../"

# --- SETTINGS ---
sam_1 = "RIKWRVLLYRGHRFAKLGMKVIK"
dnv5 = "KRRWRNICGLFGKISL"
melittin = "GIGAVLKVLTTGLPALISWIKRKRQQ"
positions = [7, 9]
num_completions = 100   # number of generated sequences for hill_climber
storage_path_files = ROOT_PATH + "outputs"
storage_path_plots= ROOT_PATH + "outputs/plots/"
# ----------------


# generate sequences and predict
def generate_and_predict():
    """
    An example of generating a new sequences. Deprecated.
    """
    sequences = generate_completions(dnv5, positions, num_completions=num_completions)

    features = []
    for sequence in sequences:
        feat = calculator.macrel_descriptors_from_seq(sequence)
        features.append(feat)

    macrel = MacrelPredictor()
    predictions = macrel.predict_seqs(features)

    d = {"Sequence": sequences, "Macrel_prediction": predictions}
    df = pd.DataFrame(d)
    df.to_csv(storage_path_files + "dnv5_mutants.csv", index=True)

def execute():
    """
    Example of usage.
    """
    seqs = [dnv5, melittin, sam_1]
    for seq in seqs:
        calculator.count_helices(seq, verbose=True)
        result = calculator.calculate_hydrophobic_moment(seq)
        print(f"Hydrophobic Moment: {result:.3f}")
        plotter.show_helical_wheel(seq, storage_path_plots + f"{seq}.png")

    # calculation of descriptors
    # descriptors_df = calculator.peptides_descriptors_from_fasta("../inputs/sequences.fasta")
    # print(descriptors_df.head(10))


execute()



