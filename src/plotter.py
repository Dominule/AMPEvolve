from modlamp.plot import helical_wheel

# Generate the plot
# moment=True will draw an arrow showing the direction of the hydrophobic moment
def show_helical_wheel(seq: str, path):
    helical_wheel(seq, moment=True, filename=path)
    print("Helical wheel saved as melittin_wheel.png")