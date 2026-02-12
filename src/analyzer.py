"""
Combines calculations from calculator.py and predictor.py,
generates pd.DataFrame with seqs and additional information relevant for AMPs.
"""

import calculator
import predictor

def analyze_AMPs(seqs: list[str], names=[]) -> pd.DataFrame:
    if not names:
        names = range(0, len(seqs))

    df = calculator.amps_analysis(seqs)

    macrel = predictor.MacrelPredictor()
    predictions = macrel.calculate_and_predict_seqs(seqs)

    df["Macrel"] = predictions

    return df