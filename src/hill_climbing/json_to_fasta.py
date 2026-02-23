import argparse
import json


def json_to_fasta(json_path, fasta_path):
    with open(json_path, 'r') as f:
        data = json.load(f)
    
    with open(fasta_path, 'w') as f:
        for i, results_group in enumerate(data):
            # The input is a list of HillClimbingResults, each has a 'results' field
            # which is a list of HillClimbingResult
            results_list = results_group.get('results', [])
            for j, result in enumerate(results_list):
                sequence = result.get('sequence')
                if sequence:
                    # Header format: group_{i}_step_{j}
                    f.write(f">group_{i}_step_{j}\n")
                    f.write(f"{sequence}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert HillClimbingResults JSON to FASTA")
    parser.add_argument("input_json", help="Path to the input JSON file")
    parser.add_argument("output_fasta", help="Path to the output FASTA file")
    
    args = parser.parse_args()
    json_to_fasta(args.input_json, args.output_fasta)
