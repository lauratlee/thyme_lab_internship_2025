# makes a subfolder in each gpcr directory for all outputs and consequent data from using the "align" pymol command

import os

for gpcr_dir in os.listdir("."):
  os.chdir(gpcr_dir)
  os.makedirs("cmd_align")
  print(f"{gpcr_dir} align folder created")
  os.chdir("..")
