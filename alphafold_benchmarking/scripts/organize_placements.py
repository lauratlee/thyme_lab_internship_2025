#organizes residue folder placements into groups of 1000 by making folders named in multiples (0, 1, 2, ...n(x1000)) and moving files into these new folders
#run from a system directory, e.g. system_dir_h_bonds/

import os

for dir in os.listdir(os.getcwd()):
  #filter for residue placement folders
  if "res_" in dir:
    os.chdir(f"{dir}/")

    #iterate thru files in residue folder to count total number of placements
    count = 0
    for r,d,f in os.walk(os.getcwd()):
      for file in f:
        if ".pdb" in file:
          count += 1

    print(f"{count} placements in {dir}")
  os.chdir("..")












