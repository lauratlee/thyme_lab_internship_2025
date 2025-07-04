# takes a zf pocket file, a human centers file, and a threshold and generates a center coordinates file that is then used to search for similar residues
# example usage (run from desired output folder, e.g. /outputs/blast/2.0): python ../../scripts/single_alignment_search.py [zf pocket file] [human center file] [threshold] [output file name]

import sys, os, csv, math
from collections import defaultdict

# --------------------------- GETTING CENTER COORDINATES ----------------------------------- #

zf_center_file = sys.argv[4]
threshold = float(sys.argv[3])
human_center_file = sys.argv[2]

# Dictionary to hold coordinates per residue
# Key = (chain, resi, resn), Value = list of (x, y, z) tuples
residue_coords = defaultdict(list)

file_path = sys.argv[1]
with open(file_path, "r") as f:
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
with open(zf_center_file, "w") as out:
  for (chain, resi, resn), coords in residue_coords.items():
      n = len(coords)
      if n == 0:
          continue

      x_avg = sum(x for x, _, _ in coords) / n
      y_avg = sum(y for _, y, _ in coords) / n
      z_avg = sum(z for _, _, z in coords) / n

      line = f"{resn} {resi} (chain {chain}): center = ({x_avg:.3f}, {y_avg:.3f}, {z_avg:.3f})\n"
      out.write(line)

print(f"Centers calculated.")

# --------------------------- COMPUTING % SIMILARITY ----------------------------------- #

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




center_files = [human_center_file, zf_center_file]
file_a = center_files[0]
file_b = center_files[1]

residues_a = parse_centers(file_a)
residues_b = parse_centers(file_b)

gene_a = os.path.basename(file_a)[:-len("_centers.txt")]
gene_b = input("Enter zf gene name:")

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










