# file intended to output the center x-, y-, and z- coordinates of a molecule.

import sys
from collections import defaultdict

pdb_file = sys.argv[1]

# Dictionary to hold coordinates per residue
# Key = (chain, resi, resn), Value = list of (x, y, z) tuples
residue_coords = defaultdict(list)

with open(pdb_file, "r") as f:
    for line in f:
        if line.startswith(("ATOM", "HETATM")):
            atom_name = line[12:16].strip()
            element = line[76:78].strip().upper()

            # Skip hydrogens
            if element == "H" or atom_name.startswith("H"):
                continue

            chain = line[21]
            resi = line[22:26].strip()
            resn = line[17:20].strip()

            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])

            key = (chain, resi, resn)
            residue_coords[key].append((x, y, z))

# Compute and center of each residue and write to a .txt file
with open(output_file, "w") as out:
    for (chain, resi, resn), coords in residue_coords.items():
        n = len(coords)
        if n == 0:
            continue

        x_avg = sum(x for x, _, _ in coords) / n
        y_avg = sum(y for _, y, _ in coords) / n
        z_avg = sum(z for _, _, z in coords) / n

        line = f"{resn} {resi} (chain {chain}): center = ({x_avg:.3f}, {y_avg:.3f}, {z_avg:.3f})\n"
        out.write(line)
        
print(f"{pdb_file} done")
