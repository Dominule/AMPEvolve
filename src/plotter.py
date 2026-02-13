from modlamp.plot import helical_wheel

# Generate the plot
# moment=True will draw an arrow showing the direction of the hydrophobic moment
def save_helical_wheel(seq: str, path):
    helical_wheel(seq, moment=True, filename=path)


def save_helical_wheels(seqs: list[str], path, names=[]):
    if not names:
        names = range(0, len(seqs))
    for seq, name in zip(seqs, names):
        if len(seq)>1:
            save_helical_wheel(seq, path + f"seq_{name}")


