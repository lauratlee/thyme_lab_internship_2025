#organizes residue folder placements into groups of 1000 by making folders named in multiples (0, 1, 2, ...n(x1000)) and moving files into these new folders
#run from a system directory, e.g. system_dir_h_bonds/

import os

for r1,d1,f1 in os.walk(os.getcwd()):
  #filter for residue placement folders
  if "res_" in d1:
    os.chdir(d1)

    #iterate thru files in residue folder to count total number of placements
    count = 0
    for r2,d2,f2 in os.walk(os.getcwd()):
      for file in f2:
        if ".pdb" in file:
          count += 1

    print(f"{count} placements in {d1}")
  os.chdir("..")












