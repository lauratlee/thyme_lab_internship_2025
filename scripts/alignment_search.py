# script generated in ChatGPT and adapted.
# iterates through each gpcr directory and creates .csv files of common residues within a specified distance (in angstroms)

import sys
import os
import math
import csv
from glob import glob

def parse_centers(filename):
    residues = []
    with open(filename, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 8:
                continue
            resn = parts[0]
            resi = parts[1]
            chain = parts[3].strip("():")
            coord_str = line.split("center = ")[-1].strip()
            coord_str = coord_str.strip("()")
            x_str, y_str, z_str = coord_str.split(",")
            x, y, z = float(x_str), float(y_str), float(z_str)
            residues.append({
                "resn": resn,
                "resi": resi,
                "chain": chain,
                "coords": (x, y, z)
            })
    return residues

def distance(c1, c2):
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(c1, c2)))

def main():
    # Get threshold from command line argument, default to 2.0 if missing or invalid
    try:
        threshold = float(sys.argv[1])
    except (IndexError, ValueError):
        print("Warning: Invalid or missing threshold argument. Using default threshold = 2.0 Å.")
        threshold = 2.0

    center_files = sorted(glob("*_centers.txt"))
    if len(center_files) < 2:
        print("Need at least two *_centers.txt files in the current directory.")
        return

    for i in range(len(center_files)):
        for j in range(i + 1, len(center_files)):
            file_a = center_files[i]
            file_b = center_files[j]

            residues_a = parse_centers(file_a)
            residues_b = parse_centers(file_b)

            output_csv = f"{os.path.splitext(file_a)[0]}-{os.path.splitext(file_b)[0]}_alignments.csv"

            with open(output_csv, "w", newline='') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow([f"Residue_in_{os.path.splitext(file_a)[0]}", f"Closest_Residue_in_{os.path.splitext(file_b)[0]}"])

                for res_a in residues_a:
                    closest_residue = None
                    min_dist = float('inf')

                    for res_b in residues_b:
                        dist = distance(res_a["coords"], res_b["coords"])
                        if dist < min_dist:
                            min_dist = dist
                            closest_residue = res_b

                    res_a_str = f"{res_a['resn']} {res_a['resi']} (chain {res_a['chain']})"

                    if min_dist <= threshold:
                        res_b_str = f"{closest_residue['resn']} {closest_residue['resi']} (chain {closest_residue['chain']})"
                    else:
                        res_b_str = "NO ALIGNED RESIDUE"

                    writer.writerow([res_a_str, res_b_str])

            print(f"Wrote alignment file: {output_csv}")

if __name__ == "__main__":
    main()
