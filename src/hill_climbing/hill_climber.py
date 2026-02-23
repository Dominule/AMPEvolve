import json
import random
from pathlib import Path

import argparse
from joblib import Parallel, delayed
from pydantic import BaseModel, Field
from pydantic_core import to_jsonable_python

from calculator import AMPKillerPredictor
from generator import aabet_without_C
from hill_climbing.json_to_fasta import json_to_fasta
from predictor import MacrelPredictor, Predictor


class HillClimbingResult(BaseModel):
    sequence: str
    score: float = 0.0
    improvement: float = 0.0


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

    def optimize_sequence(self, sequence: str, verbose=False) -> HillClimbingResults:
        best_sequence = sequence
        step_results = [HillClimbingResult(sequence=sequence,
                                           score=self.scorer.calculate_and_predict_seqs([sequence])[0],
                                           improvement=0) ]

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

    def optimize_sequence_just_string(self, sequence: str, verbose=False) -> str:
        return self.optimize_sequence(sequence, verbose).results[-1].sequence


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Hill Climbing Optimization for Peptide Sequences")
    parser.add_argument('--input_seqs', type=str, help='A fasta file with starting sequences')
    parser.add_argument('--same_as', type=str,
                        help='JSON file containing list of HillClimbingResults to use first sequence from each')
    parser.add_argument('--num_sequences', type=int, default=5,
                        help='Number of sequences to generate (ignored if --input_seqs or --same_as provided)')
    parser.add_argument('--min_length', type=int, default=18,
                        help='Minimum length of generated sequences (ignored if --input_seqs or --same_as provided)')
    parser.add_argument('--max_length', type=int, default=25,
                        help='Maximum length of generated sequences (ignored if --input_seqs or --same_as provided)')
    parser.add_argument('--output', type=str, default="data", help='Output folder path for results')

    args = parser.parse_args()

    # Determine starting sequences
    starting_sequences = []

    if args.input_seqs:
        # Load sequences from fasta file
        with open(args.input_seqs, 'r') as f:
            starting_sequences = [line.strip() for line in f if line.strip() if line[0]!=">"]
    elif args.same_as:
        # Load sequences from JSON file of HillClimbingResults
        with open(args.same_as, 'r') as f:
            results_list = json.load(f)
            for result_data in results_list:
                if isinstance(result_data, dict) and 'results' in result_data and result_data['results']:
                    starting_sequences.append(result_data['results'][0]['sequence'])
    else:
        # Generate random sequences
        starting_sequences = [
            "".join(random.choices(aabet_without_C, k=random.randint(args.min_length, args.max_length)))
            for _ in range(args.num_sequences)
        ]

    hill_climber = HillClimber(scorer=AMPKillerPredictor())

    results = Parallel(n_jobs=-1, backend="threading", verbose=10)(
        delayed(hill_climber.optimize_sequence)(seq)
        for seq in starting_sequences
    )

    output_file_name = f"hill_climber_results_{hill_climber.scorer.__class__.__name__}.json"
    outputs_directory = Path(args.output or ".")
    outputs_directory.mkdir(exist_ok=True, parents=True)
    results_file_path = outputs_directory / output_file_name

    with open(results_file_path, "w") as f:
        json.dump(results, f, default=to_jsonable_python, indent=2)
    json_to_fasta(results_file_path, results_file_path.with_suffix(".fasta"))
    print(f"Saved {len(results)} results to {results_file_path.resolve()}")
