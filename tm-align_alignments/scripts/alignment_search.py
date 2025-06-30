# script generated in ChatGPT and adapted.
# iterates through each gpcr directory and creates .csv files of common residues within a specified distance (in angstroms)
# example usage: python ../scripts/alignment_search.py Class_A 2.0

import sys
import os
import math
import csv
from glob import glob

class_name = sys.argv[1]

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

def process_subfolder(folder, threshold): 
    center_files = sorted(glob.glob("*_centers.txt"))

    if len(center_files) < 2:
        print(f"[{folder}] Skipping: fewer than 2 *_centers.txt files.")
        return

    for i in range(len(center_files)):
        for j in range(i + 1, len(center_files)):
            file_a = center_files[i]
            file_b = center_files[j]

            residues_a = parse_centers(file_a)
            residues_b = parse_centers(file_b)

            gene_a = file_a[:-len("_centers.txt")]
            gene_b = file_b[:-len("_centers.txt")]

            output_csv = f"{gene_a}-{gene_b}_{threshold:.1f}.csv"
            output_path = os.path.join(folder, output_csv)

            with open(output_path, "w", newline='') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow([
                    f"{gene_a}_pocket_residue",
                    f"closest {gene_b}_pocket_residue within {threshold:.1f} angstroms"
                ])

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

            print(f"[{folder}] Wrote: {output_csv}")    

def main():
    gpcr_class_list = ["Class_A", "Class_B1", "Class_C", "Class_F"]
    try:
        threshold = float(sys.argv[2])
    except (IndexError, ValueError):
        print("Warning: Invalid or missing threshold argument. Using default threshold = 2.0 Å.")
        threshold = 2.0

    for sub in os.listdir("."):
      print(f"Running process_subfolder on {sub}/{class_name}")
      process_subfolder(os.path.join(sub, class_name), threshold)


while True:
  answer = input(f"Run alignments for {class_name}? [y/n]").strip().lower()
  if answer == "y":
    main()
    break
  elif answer == "n":
    print("Exiting program.")
    sys.exit(0)
  else:
    print("Invalid input, please answer y/n")












