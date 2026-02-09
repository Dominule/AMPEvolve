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
from calculator import get_peptide_features
from predictor import MacrelPredictor
from hill_climber import climb_high

ROOT_PATH = "../"

# --- SETTINGS ---
dnv5 = "KRRWRNICGLFGKISL"
positions = [7, 9]
num_completions = 100   # number of generated sequences
storage_path = ROOT_PATH + "outputs/dnv5_mutants.csv"
# ----------------


# generate sequences and predict
def generate_and_predict():
    # deprecated
    sequences = generate_completions(dnv5, positions, num_completions=num_completions)

    features = []
    for sequence in sequences:
        feat = get_peptide_features(sequence)
        features.append(feat)

    macrel = MacrelPredictor()
    predictions = macrel.predict_seqs(features)

    d = {"Sequence": sequences, "Macrel_prediction": predictions}
    df = pd.DataFrame(d)
    df.to_csv(storage_path, index=True)

def execute():
    # climb_high("KRRWRQVMGAFWKIKV", [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15], 18, 75, until_finished=True) # DNV5 tuning

    # file = "../inputs/selected_by_sam.fasta"
    # generated = loader.load_fasta(file)
    file = "../inputs/generated_AMPGen_3.csv"
    generated = loader.load_csv(file, column_name="Sequence")
    # generated = "MSFREMLSNSCPRSNGGG"
    macrel = MacrelPredictor()
    predictions = macrel.calculate_and_predict_seqs(generated)
    print(generated)
    print(predictions)

    # generated = climb_high(generated, [], until_finished=True, mask_all=True)
    # predictions = macrel.calculate_and_predict_seqs([generated])
    # print(generated)
    # print(predictions)


execute()



