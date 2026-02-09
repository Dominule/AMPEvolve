"""
This module goes through the big fat AMPSphere predictions
and critically extracts all potential amps.
"""

import pandas as pd

FILE_PATH = "../../inputs/AMPSphere_all_predicted.csv" # file available on zenodo: http://zenodo.org/records/15471481
FILTERED_PATH = "../../outputs/AMPSphere_filtered.csv"

def read_ampsphere():
    df = pd.read_csv(FILE_PATH)
    mask = df.apply(lambda row: (row["HIGH_PROBA"]
                          and row["PRED_01"] >= 0.87
                          and row["PRED_02"] >= 0.97
                          and row["PRED_03"] >= 0.87), axis=1)
    filtered = df[mask]
    filtered.sort_values("PRED_02")
    filtered.to_csv(FILTERED_PATH)

