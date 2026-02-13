import random

import json
from pydantic import BaseModel, Field
from pydantic_core import to_jsonable_python
from tqdm import tqdm

from generator import aabet_without_C
from predictor import MacrelPredictor, Predictor


class HillClimbingResult(BaseModel):
    sequence: str
    score: float = 0.0
    improvement: float = Field(default=0.0, ge=0.0)


class HillClimbingResults(BaseModel):
    results: list[HillClimbingResult]

    def __len__(self):
        return len(self.results)

    def __getitem__(self, index):
        return self.results[index]

    def __iter__(self):
        return iter(self.results)


class HillClimber(BaseModel):
    alphabet: list[str] = aabet_without_C
    change_multiple: bool = False
    scorer: Predictor = MacrelPredictor()
    epochs: int = Field(default=100, gt=0)

    class Config:
        arbitrary_types_allowed = True

    def do_one_step(self, original_seq: str) -> HillClimbingResult:
        original_score = self.scorer.calculate_and_predict_seqs([original_seq])[0]
        best_score = original_score
        best_seq = original_seq
        for position in range(len(original_seq)):
            for letter in self.alphabet:
                if self.change_multiple:
                    new_seq = best_seq[:position] + letter + best_seq[position + 1:]
                else:
                    new_seq = original_seq[:position] + letter + original_seq[position + 1:]
                new_score = self.scorer.calculate_and_predict_seqs([new_seq])[0]
                if new_score > best_score:
                    best_score = new_score
                    best_seq = new_seq
        return HillClimbingResult(sequence=best_seq, score=best_score, improvement=best_score - original_score)

    def fit(self, sequence: str, verbose=False) -> HillClimbingResults:
        best_sequence = sequence
        step_results = []
        for epoch in range(self.epochs):
            result = self.do_one_step(best_sequence)
            best_sequence = result.sequence
            if verbose:
                print(f"Epoch {epoch + 1}: Best sequence score: {result.score:.6f} ({best_sequence})")
            if result.improvement == 0:
                if verbose: print(f"Converged at epoch {epoch + 1}")
                break
            else:
                step_results.append(result)
        return HillClimbingResults(results=step_results)


if __name__ == '__main__':
    hill_climber = HillClimber()
    results = [hill_climber.fit("".join(random.choices(aabet_without_C, k=random.randint(18, 25)))) for _ in
               tqdm(list(range(250)))]
    with open("hill_climber_results.json", "w") as f:
        json.dump(results, f, default=to_jsonable_python, indent=2)

    print(f"Saved {len(results)} results to hill_climber_results.json")

# Epoch 1: Best sequence score: 0.38 (FHGLWQIGVQTYPPKSAYIA)
# Epoch 2: Best sequence score: 0.59 (FHGLWRIGVQTYPPKSAYIA)
# Epoch 3: Best sequence score: 0.78 (FHGLWRIGVKTYPPKSAYIA)
# Epoch 4: Best sequence score: 0.88 (FLGLWRIGVKTYPPKSAYIA)
# Epoch 5: Best sequence score: 0.95 (FLGLWRIGVKIYPPKSAYIA)
# Epoch 6: Best sequence score: 0.96 (FVGLWRIGVKIYPPKSAYIA)
# Epoch 7: Best sequence score: 0.98 (FVGLWRIGVKIYPPKSAGIA)
# Epoch 8: Best sequence score: 0.99 (FVGLMRIGVKIYPPKSAGIA)
# Epoch 9: Best sequence score: 1.00 (FVGLMRIFVKIYPPKSAGIA)
