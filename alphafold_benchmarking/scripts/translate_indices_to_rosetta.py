# takes a .pdb of anchor residues and a rosetta-indiced .pdb file, and prints the rosetta-translated indices of the anchor residues
# example usage: python path/to/translate_indices_to_rosetta.py anchors.pdb chain_A_ligand_0.pdb


import sys, os
from collections import defaultdict

anchor_file = sys.argv[1]
rosetta_file = sys.argv[2]

def get_centers(pdb_file):
    # Dictionary to hold coordinates per residue
    # Key = (chain, resi, resn), Value = list of (x, y, z) tuples
    residue_coords = defaultdict(list)
    center_dict = defaultdict(list)

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

    # Compute center of each residue
    for (chain, resi, resn), coords in residue_coords.items():
        n = len(coords)
        if n == 0:
            continue

        x_avg = sum(x for x, _, _ in coords) / n
        y_avg = sum(y for _, y, _ in coords) / n
        z_avg = sum(z for _, _, z in coords) / n

        res_id = (resi, resn)
        center_dict[res_id].append((round(x_avg, 3), round(y_avg, 3), round(z_avg, 3)))

    return center_dict


anchor_centers = get_centers(anchor_file)
rosetta_centers = get_centers(rosetta_file)

rosetta_idxs = []

for (a_resi, a_resn), a_coords in anchor_centers.items():
    for (r_resi, r_resn), r_coords in rosetta_centers.items():
        if (a_resn == r_resn) and (a_coords == r_coords):
            rosetta_idxs.append(r_resi)

print(rosetta_idxs)













      






