# given a .pdb file, gets the FASTSA sequence and saves it as a .fasta file
# example usage (to be run from gpcr_pocket_dir/): python ../scripts/extract_fasta.py ACKR3/Class_A/ackr3a_pocket.pdb

import os, sys, pymol2

pdb_file = sys.argv[1]

pocket_name = os.path.splitext(os.path.basename(pdb_file))[0]

out_file = f"{pocket_name}.fasta"
out_path = os.path.join(os.path.dirname(pdb_file), out_file)

# get sequence of pocket
with pymol2.PyMOL() as pymol:
  # set the internal gui width
  pymol.cmd.set('internal_gui_width', 600)

  pymol.cmd.load(pdb_file, "pocket")
  sequence = pymol.cmd.get_fastastr("pocket")


with open(out_path, "w") as f:
  f.write(sequence)

print(f"Saved to {out_path}")
  
