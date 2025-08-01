#iterates through residue folders within a system and deletes any placements with a ddG >= 0
#run from a system folder, e.g. 9HZ0/

import os


#helper function for cleaning
def clean_ddgs(residue_folder):
  #add group folders to look at to a list
  group_list = []

  
  for folder in os.listdir(os.getcwd()):
    if os.path.isdir(folder):
      try:
        int(folder)
        group_list.append(folder)
      except ValueError:
        continue

  #print group list
  print(group_list)

  #iterate through groups for cleaning
  for group in group_list:
    print(group)
    os.chdir(group)
    for file in os.listdir(os.getcwd()):
      #confirm that file is a placement
      if ".pdb" not in file:
        continue

      #open pdb and parse
      with open(file, 'r') as placement:
        for line in placement:
          if "Scoring: Post-HighResDock system ddG:" in line:
            parts = line.strip().split()
            ddg_value = float(parts[-1])
            break

      if ddg_value is not None and ddg_value >= 0:
        os.remove(file)

    #exit group and move on to the next
    os.chdir("..")
        



#add residue folders to look at to a list
residue_list = []

for dir in os.listdir(os.getcwd()):
  if "res_" in dir:
    residue_list.append(dir)

for res in residue_list:
  print(res)

  #enter residue folder
  os.chdir(res)

  #clean out files with ddgs >= 0
  clean_ddgs(res)

  #exit residue folder to move on to the next residue
  os.chdir("..")





