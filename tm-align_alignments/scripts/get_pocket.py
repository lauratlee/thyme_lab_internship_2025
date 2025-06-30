# collects binding pocket residues within 5 angstroms of a reference ligand
# example usage (run from gpcr_pocket_dir): python ../scripts/get_pocket.py Class_A

import os, sys, re
import pymol2

class_map = {
  "Class_A": ("2rh1_chainA.pdb", "cau"),
  "Class_B1": ("4k5y_chainA.pdb", "1q5"),
  "Class_C": ("7m3g_chainA.pdb", "h43"),
  "Class_F": ("4jkv_chainA.pdb", "1ks")
}

gpcr_class = sys.argv[1]

if gpcr_class in class_map:
  ref_file, ligand_name = class_map[gpcr_class]
else:
  print("ERROR: Please provide a valid gpcr class.")
  sys.exit(1)

# helper function to extract pocket residues
def find_pocket(gpcr_dir, reference, ligand):
  for root, _, files in os.walk(gpcr_dir):
    for file in files:
      if not file.endswith("pdb.pdb"):
        continue
      print(f"Found alignment file: {file}")
      name = file.removesuffix(".pdb.pdb")
      print(f"GENE NAME: {name}")

      # initialize pymol2 in headless mode
      with pymol2.PyMOL() as pymol:
        # set the internal gui width
        pymol.cmd.set('internal_gui_width', 600)

        # load in reference
        pymol.cmd.load(f"~/thyme_lab_internship_2025/gpcr_class_reps/{ref_file}", "ref")
        if pymol.cmd.count_atoms("ref") == 0: 
          print("[WARNING] Failed to load ref structure")
          sys.exit(1)

        # load in target
        pymol.cmd.load(f"~/thyme_lab_internship_2025/tm-align_alignments/gpcr_pocket_dir/{gpcr_dir}/{gpcr_class}/{file}", "target")
        if pymol.cmd.count_atoms("target") == 0:
          print(f"[WARNING] Failed to load target structure: {file}")
          sys.exit(1)

        # select reference ligand 
        pymol.cmd.select("ligand", f"resn {ligand_name}")
        if pymol.cmd.count_atoms("ligand") == 0:
          print(f"[WARNING] No atoms found for ligand {ligand_name}")

        # locate and select pocket residues of target
        pocket_sele = "byres (target within 5 of ligand) and target"
        pymol.cmd.select("pocket", pocket_sele)
        print(f"Atom count in pocket: {pymol.cmd.count_atoms('pocket')}")

        # save pocket
        output_dir = f"~/thyme_lab_internship_2025/tm-align_alignments/gpcr_pocket_dir/{gpcr_dir}/{gpcr_class}/"
        os.makedirs(output_dir, exist_ok = True)
        output_file_name = f"{name}_pocket.pdb"
        output_path = os.path.join(output_dir, output_file_name)
        pymol.cmd.save(output_path, "pocket")
        print(f"saved {output_file_name}")
        



def main():
  #walk through gpcr directory and run find_pocket on each gpcr subdir
  for sub in os.listdir("."):
      print(f"Calling find_pocket on: {sub}")
      find_pocket(sub, ref_file, ligand_name)



