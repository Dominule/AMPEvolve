import json
import argparse
import csv
from pathlib import Path

from pydantic_core import to_jsonable_python

from hill_climbing.hill_climber import HillClimbingResult, HillClimbingResults


def tsv_to_json(tsv_path, fasta_path, output_json_path):
    # Read fasta for name to seq dictionary
    name_to_seq = {}
    with open(fasta_path, 'r') as f:
        text = f.read()
        for lines in [line.split('\n') for line in text.split('>')[1:]]:
            name, seq = lines[0], lines[1]
            name_to_seq[name.strip()] = seq.strip()

    results = []
    current_run = []
    current_group = str(0)
    # Read TSV and update scores
    with open(tsv_path, 'r', newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if not row:
                continue
            # first column: name of the sequence
            # second column: ignored
            # third column: Active or Not Active
            # forth column: confidence
            name = row[0]
            group, run = name.split('_')[1], name.split('_')[2]
            if group != current_group:
                results.append(current_run)
                current_run = []
                current_group = group
            status = row[2]
            try:
                confidence = float(row[3])
            except (IndexError, ValueError):
                confidence = 0.0
            if status == "Active":
                score = confidence
            else:
                score = 1-confidence
            current_run.append(HillClimbingResult(sequence=name_to_seq[name], score=score))
    results.append(current_run)
    # Save the updated JSON
    with open(output_json_path, 'w') as f:
        json.dump([HillClimbingResults(results=x) for x in results], f, indent=2, default=to_jsonable_python)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert TSV back to HillClimbingResults JSON")
    parser.add_argument("input_tsv", help="Path to the input TSV file")
    parser.add_argument("original_fasta", help="Path to the original fasta file (to match sequences)")
    parser.add_argument("output_json", help="Path to the output JSON file")

    args = parser.parse_args()
    tsv_to_json(args.input_tsv, args.original_fasta, args.output_json)
