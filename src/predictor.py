"""
Loads macrel onnx model and predicts peptide properties.
"""

import gzip
from abc import abstractmethod, ABC

import onnxruntime as ort
import numpy as np
from pathlib import Path

import calculator

MODEL_PATH = Path(__file__).parents[1] / "models/macrel.onnx.gz"


class Predictor(ABC):
    @abstractmethod
    def calculate_and_predict_seqs(self, sequences: list[str]) -> list[float]:
        pass


class MacrelPredictor(Predictor):
    def __init__(self, model_path=MODEL_PATH):
        with gzip.open(model_path, 'rb') as f:
            model_bytes = f.read()

        self._model = ort.InferenceSession(model_bytes, providers=['CPUExecutionProvider'])
        self._input_name = self._model.get_inputs()[0].name

    def predict_seq(self, features: np.ndarray) -> np.ndarray:
        """
        Predicts peptide properties using the macrel ONNX model.
        Args:
            features (np.ndarray): Array of shape (n_samples, n_features) containing peptide features.
        Returns:
            np.ndarray: Predicted properties.
        """
        inputs = {self._input_name: features}
        outputs = self._model.run(None, inputs)
        return outputs[1]

    def predict_seqs(self, features_list: list[np.ndarray]) -> list:
        """
        Predicts peptide properties for a list of feature arrays.
        Args:
            list of np.ndarrays (features)
        Returns:
            list of floats (proba of AMP)
        """
        all_features = np.vstack(features_list)
        raw_results = self.predict_seq(all_features)
        return [round(p["AMP"], 2) for p in raw_results]

    def calculate_and_predict_seq(self, sequence: str) -> float:
        """
        Calculates descriptors and predicts from aa sequence.
        Args:
            str (seq)
        Returns:
            float (proba of AMP)
        """
        features = calculator.macrel_descriptors_from_seq(sequence)
        return round(self.predict_seq(features)[0]["AMP"], 3)

    def calculate_and_predict_seqs(self, sequences: list[str]) -> list[float]:
        """
        Calculates descriptors and predicts from aa sequences.
        Args:
            list of str (seqs)
        Returns:
            list of floats (proba of AMP)
        """
        features_list = [calculator.macrel_descriptors_from_seq(seq) for seq in sequences]
        return self.predict_seqs(features_list)
