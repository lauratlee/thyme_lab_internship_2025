#this is the 2nd step after running residue-specific rmsd calculations within a system (mainly for 9M3N, 9QEK, 9QEL)
#run from system folder, e.g. 9M3N/
#example usage: python /path/to/script/ res_166

import os, sys, csv, ast

system = os.path.basename(os.getcwd())
residue = sys.argv[1]

output_path = f"{system}_{residue}_placements_summary.csv"

rmsd_list = []

for file in os.listdir(os.getcwd()):
  if "res_" in file and "group" in file:
    print(file)
    with open(file, newline="") as rmsd_file:
      reader = csv.reader(rmsd_file)
      header = next(reader, None) #skip header
      rows = list(reader)
      if not rows:
        print(f"WARNING: {file} is empty")
        continue
      for row in rows:
        rmsd_list.append([row[0], row[1], float(row[2])])

#sort best_rmsds by rmsd value
rmsd_list.sort(key=lambda x: x[2])
best_system_rmsd = rmsd_list[0]

#write to system summary file
with open(f"{system}_{residue}_placements_summary.csv", "w") as system_file:
  system_file.write("residue,file,rmsd\n")
  system_file.write(f"{best_system_rmsd[0]},{best_system_rmsd[1]},{best_system_rmsd[2]}\n")
  
          
