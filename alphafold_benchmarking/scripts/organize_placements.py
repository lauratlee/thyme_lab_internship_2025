#organizes residue folder placements into groups of 1000 by making folders named in multiples (0, 1, 2, ...n(x1000)) and moving files into these new folders
#run from a system directory, e.g. system_dir_h_bonds/

import os, math

for r1,d1,f1 in os.walk(os.getcwd()):
  for dir1 in d1:
    #filter for residue placement folders
    if "res_" in dir1:
      os.chdir(dir1)
  
      #iterate thru files in residue folder to count total number of placements
      count = 0
      for r2,d2,f2 in os.walk(os.getcwd()):
        for file2 in f2:
          if ".pdb" in file2:
            count += 1
  
      print(f"{count} placements in {dir1}")

      #calculate total # of groups of 1000
      groups = math.ceil(count / 1000)

      #create directories to sort placements into
      for i in range(0, groups):
        os.system(f"mkdir {i}")

      #iterate through files again to sort them into their respective folders
      for r3,d3,f3 in os.walk(os.getcwd()):
        for file3 in f3:
          if ".pdb" not in file3:
            continue

          #split a file name like "chain_A_ligand_0.pdb" to extract "0.pdb", then split again to extract the placement index (in this case 0)
          index = file3.split("_")[-1].split(".")[0]
          group_id = index // 1000

          os.system(f"mv {file3} {group_id}/{file_3}")
    


      #let user know that residue is done
      print(f"{dir1} DONE")

    os.chdir("..")

  #when entire system is done, let user know
  print("all residues in system done")












