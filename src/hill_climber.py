"""
Here is code for lead optimization.
"""
import pandas as pd

from generator import generate_completions
from calculator import get_peptide_features
from predictor import MacrelPredictor


def climb_high(seq: str, positions: list[int], num_completions: int, epochs=20) -> str:
    counter = 0
    while counter < epochs:

        # Mutants
        all_mutants = []
        for pos in positions:
            mutants = generate_completions(seq, [pos], num_completions)
            all_mutants = all_mutants + mutants
        seqs = all_mutants + [seq]

        # Features
        features = []
        for seq in seqs:
            feat = get_peptide_features(seq)
            features.append(feat)

        # Predictions
        macrel = MacrelPredictor()
        predictions = macrel.predict_seqs(features)
        d = {"Sequence": seqs, "Pred": predictions}
        df = pd.DataFrame(d)

        # Winner
        highest_value = df['Pred'].max()
        new_seq = df[df['Pred'] == highest_value]['Sequence'].values[0]
        print(new_seq)

        seq = new_seq
        counter += 1

    return seq