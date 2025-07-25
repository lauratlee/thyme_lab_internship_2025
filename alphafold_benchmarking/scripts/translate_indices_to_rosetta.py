# takes a .pdb of anchor residues and a rosetta-indiced .pdb file, and prints the rosetta-translated indices of the anchor residues
# example usage: python path/to/translate_indices_to_rosetta.py anchors.pdb chain_A_ligand_0.pdb


import sys, os
from collections import defaultdict
import math

anchor_file = sys.argv[1]
rosetta_file = sys.argv[2]

THRESHOLD = 0.5

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
        center_dict[res_id].append((x_avg, y_avg, z_avg))

    return center_dict

def euclidean_distance(coord1, coord2):
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(coord1, coord2)))

anchor_centers = get_centers(anchor_file)
rosetta_centers = get_centers(rosetta_file)

rosetta_idxs = []

for (a_resi, a_resn), a_coords_list in anchor_centers.items():
    for a_coord in a_coords_list:
        for (r_resi, r_resn), r_coords_list in rosetta_centers.items():
            if a_resn != r_resn:
                continue
            for r_coord in r_coords_list:
                if euclidean_distance(a_coord, r_coord) <= THRESHOLD:
                    rosetta_idxs.append(r_resi)

print(rosetta_idxs)













      






