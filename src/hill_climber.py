"""
Here is code for lead optimization.
"""
import pandas as pd

from generator import generate_completions
import calculator
from predictor import MacrelPredictor


def climb_high(seq: str, positions: list[int], num_completions=10, epochs=10, mask_all=False, until_finished=False, verbose=False) -> str:
    """
    Optimizes peptide sequence with hill climbing.
        Args:
            seq (str): Peptide sequence (e.g., AWKGAGIISL).
            positions (list[int]): Positions to be modified.
            num_completions (int): Number of seqs to generate.
            epochs (int): Number of epochs. For each epoch, one position is modified.
            mask_all (bool): If true, ignores "positions" and masks all positions.
            until_finished (bool): If true, ignores "epochs" and generates until new sequence is equal to the previous one.
            verbose (bool): If true, prints new seq and its AMP proba every epoch.
        Returns:
            str: Optimalized sequence.
    """
    counter = 0
    while (counter < epochs) or until_finished:
        if mask_all:
            positions = range(0, len(seq))

        # Mutants
        all_mutants = []
        for pos in positions:
            mutants = generate_completions(seq, [pos], num_completions)
            all_mutants = all_mutants + mutants
        seqs = all_mutants + [seq]

        # Features
        features = []
        for seq in seqs:
            feat = calculator.macrel_descriptors_from_seq(seq)
            features.append(feat)

        # Predictions
        macrel = MacrelPredictor()
        predictions = macrel.predict_seqs(features)
        d = {"Sequence": seqs, "Pred": predictions}
        df = pd.DataFrame(d)

        # Winner
        highest_value = df['Pred'].max()
        new_seq = df[df['Pred'] == highest_value]['Sequence'].values[0]
        if verbose:
            print(new_seq)
            print(highest_value)

        if until_finished:
            if new_seq == seq:
                return seq

        seq = new_seq
        counter += 1

    return seq