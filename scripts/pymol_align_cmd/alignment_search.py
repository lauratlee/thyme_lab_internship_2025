# script generated in ChatGPT and adapted.
# iterates through each gpcr directory and creates .csv files of common residues within a specified distance (in angstroms)
# example usage: python ../scripts/alignment_search.py 2.0

import sys, os, math, csv, re
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

def process_subfolder(threshold):
    folder = os.path.basename(".")
    if not os.path.exists("cmd_align"):
        print(f"[{folder}] Skipping: 'cmd_align' folder not found.")
        os.chdir("..")
        return
    
    os.chdir("cmd_align")

    center_files = sorted(glob("*_centers.txt"))

    if len(center_files) < 2:
        print(f"[{folder}] Skipping: fewer than 2 *_centers.txt files.")
        os.chdir("../..")
        return

    for i in range(len(center_files)):
        for j in range(i + 1, len(center_files)):
            file_a = center_files[i]
            file_b = center_files[j]

            residues_a = parse_centers(file_a)
            residues_b = parse_centers(file_b)

            gene_a = os.path.splitext(file_a)[0].replace("_pocket(align)_centers", "")
            gene_b = os.path.splitext(file_b)[0].replace("_pocket(align)_centers", "")
            output_csv = f"{gene_a}-{gene_b}_{threshold:.1f}.csv"

            with open(output_csv, "w", newline='') as csvfile:
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
    os.chdir("../..")

def main():
    try:
        threshold = float(sys.argv[1])
    except (IndexError, ValueError):
        print("Warning: Invalid or missing threshold argument. Using default threshold = 2.0 Å.")
        threshold = 2.0

    subfolders = [f for f in os.listdir(".") if os.path.isdir(f)]
    if not subfolders:
        print("No subfolders found.")
        return

    for folder in sorted(subfolders):
        os.chdir(folder)
        process_subfolder(threshold)

if __name__ == "__main__":
    main()

