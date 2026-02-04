"""
Calculates features for peptide sequences using Macrel library.
"""

import numpy as np
from macrel.macrel_features import compute_all


def get_peptide_features(sequence: str) -> np.ndarray:
    """
    Extracts 22 features using the compute_all function from macrel.
    """
    features = compute_all(sequence)

    # Convert to float32 numpy array for ONNX compatibility
    return np.array(features, dtype=np.float32).reshape(1, -1)
