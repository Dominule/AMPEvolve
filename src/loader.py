"""
Loads from "<root>/data" directory files,
returns list or pd.DataFrame.
"""

import pandas as pd

def load_csv(file_path: str, column_name: str) -> list[str]:
    """
    Returns list from specified column.
    """
    df = pd.read_csv(file_path)
    return df[column_name]


def load_fasta(file_path: str) -> list[str]:
    sequences = []
    with open(file_path, 'r') as f:
        current_seq = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_seq:
                    sequences.append("".join(current_seq))
                current_seq = []
            else:
                current_seq.append(line)
        if current_seq:
            sequences.append("".join(current_seq))
    return sequences

