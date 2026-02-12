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
import generator
import calculator
import plotter
from predictor import MacrelPredictor
from hill_climber import climb_high

ROOT_PATH = "../"

# --- SETTINGS ---
sam_1 = "RIKWRVLLYRGHRFAKLGMKVIK"
dnv5 = "KRRWRNICGLFGKISL"
melittin = "GIGAVLKVLTTGLPALISWIKRKRQQ"
penetratin = "RQIKIWFQNRRMKWKK"
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
    sequences = generator.generate_completions(dnv5, positions, num_completions=num_completions)

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
    ## AMP Information and helical_wheels
    # seqs = {dnv5 : "dnv5", melittin : "melittin", sam_1:"sam_1", penetratin : "penetratin"}
    # for seq, name in seqs.items():
    #     calculator.alphahelices(seq, verbose=True)
    #     result = calculator.calculate_hydrophobic_moment(seq)
    #     print(f"Hydrophobic Moment: {result:.3f}")
    #     plotter.show_helical_wheel(seq, storage_path_plots + f"{name}.png")

    ## Calculation of descriptors
    # descriptors_df = calculator.peptides_descriptors_from_fasta("../inputs/sequences.fasta")
    # print(descriptors_df.head(10))

    # generating
    seqs = generator.generate_amphipatic_helices()
    if not seqs:
        print("No sequences found.")
    print(seqs)

    df = calculator.amps_analysis(seqs, verbose=True)
    # plotter.save_helical_wheels(seqs, storage_path_plots +"generated/")




execute()



