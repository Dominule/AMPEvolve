from pathlib import Path

from generator import generate_all_neighbors


def fasta_to_neighbors(fasta_path: Path, output_folder: Path):
    with open(fasta_path, 'r') as f:
        text = f.read()
    name = ""
    for line in text.splitlines():
        if line.startswith(">"):
            name = line[1:].strip()
        else:
            neighbors = generate_all_neighbors(line.strip())
            with open(f"{output_folder}/{name}.fasta", 'w') as f:
                f.write(f">{name}\n")
                f.write(f"{line.strip()}\n")
                for i, neighbor in enumerate(neighbors):
                    f.write(f">{name}_{i}\n")
                    f.write(f"{neighbor}\n")


if __name__ == '__main__':
    output_path = Path(r"/outputs/native_neighbors")
    output_path.mkdir(exist_ok=True)
    fasta_to_neighbors(Path(r"/outputs/natives.fasta"),
                       output_path)
