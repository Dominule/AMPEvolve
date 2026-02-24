import json
import argparse
import matplotlib.pyplot as plt
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description='Plot Hill Climbing results from JSON files.')
    parser.add_argument('files', nargs='+', help='JSON files containing HillClimbingResults')
    parser.add_argument('--output', help='Output filename (default: hill_climbing_plot.png)')
    args = parser.parse_args()

    if not args.files:
        print("No files provided.")
        return

    first_file = Path(args.files[0])
    output_dir = first_file.parent
    output_name = args.output if args.output else "hill_climbing_plot.png"
    output_path = output_dir / output_name

    plt.figure(figsize=(10, 6))

    for file_path in args.files:
        path = Path(file_path)
        with open(path, 'r') as f:
            data = json.load(f)
        
        color = None
        label_added = False
        all_scores_per_epoch = []

        for hc_results in data:
            # Handle both list of results and results key
            if isinstance(hc_results, dict) and 'results' in hc_results:
                results = hc_results['results']
            elif isinstance(hc_results, list):
                results = hc_results
            else:
                continue

            scores = [r.get('score', 0.0) for r in results]
            epochs = list(range(len(scores)))

            if not scores:
                continue

            all_scores_per_epoch.append(scores)

            if not label_added:
                line, = plt.plot(epochs, scores, label=path.stem, alpha=scores[-1])
                color = line.get_color()
                label_added = True
            else:
                plt.plot(epochs, scores, color=color, alpha=scores[-1])

        # Calculate and plot average line
        if all_scores_per_epoch:
            max_length = max(len(s) for s in all_scores_per_epoch)
            avg_scores = []
            for i in range(max_length):
                epoch_scores = [s[i] for s in all_scores_per_epoch if i < len(s)]
                avg_scores.append(sum(epoch_scores) / len(epoch_scores))
            avg_epochs = list(range(len(avg_scores)))
            plt.plot(avg_epochs, avg_scores, color=color, linewidth=2.5, linestyle='--', alpha=1.0,
                     label=f'{path.stem} (avg)')
    plt.ylabel('Score')
    plt.xlabel('Epoch')
    plt.title('Hill Climbing Progress')
    plt.legend()
    plt.grid(True)
    
    plt.savefig(output_path)
    print(f"Plot saved to {output_path}")

if __name__ == '__main__':
    main()
