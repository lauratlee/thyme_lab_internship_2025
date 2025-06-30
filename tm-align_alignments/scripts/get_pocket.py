# collects binding pocket residues within 5 angstroms of a reference ligand

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


def find_pocket(gpcr_dir, reference, ligand):
  for root, _, files in os.walk(gpcr_dir):
    
  










def main():
  #walk through gpcr directory and run find_pocket on each gpcr subdir
  for sub in os.listdir("."):
      print(f"Calling find_pocket on: {sub}")
      find_pocket(sub, ref_file, ligand_name)



