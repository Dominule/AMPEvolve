"""
Combines calculations from calculator.py and predictor.py,
generates pd.DataFrame with seqs and additional information relevant for AMPs.
"""

import pandas as pd
import calculator
import predictor

def analyze_AMPs(seqs: list[str], storage_path) -> pd.DataFrame:
    df = calculator.amps_analysis(seqs)

    macrel = predictor.MacrelPredictor()
    predictions = macrel.calculate_and_predict_seqs(seqs)

    df["Macrel"] = predictions
    df.sort_values(by="Macrel", ascending=False, inplace=True)
    df.to_csv(storage_path, index=True)

    return df