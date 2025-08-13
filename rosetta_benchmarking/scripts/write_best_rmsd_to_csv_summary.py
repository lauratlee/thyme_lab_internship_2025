#this is the 2nd step after running residue-specific rmsd calculations within a system (mainly for 9M3N, 9QEK, 9QEL)
#run from system folder, e.g. 9M3N/

import os, sys, csv, ast

system = os.path.basename(os.getcwd())
best_rmsds = []

for file in os.listdir(os.getcwd()):
  if "res_" in file and "placements_summary.csv" in file:
    with open(file, newline="") as rmsd_file:
      reader = csv.reader(rmsd_file)
      header = next(reader, None)
      rows = list(reader)

      for row in rows:
        best_rmsds.append([row[0], row[1], float(row[2])]


if best_rmsds:
  best_rmsds.sort(key=lambda x: x[2])
  best_system_rmsd = best_rmsds[0]
  print(best_system_rmsd)


with open(f"{system}_placements_summary.csv", "w") as system_file:
  system_file.write("residue, file, rmsd\n")
  system_file.write(f"{best_system_rmsd[0]},{best_system_rmsd[1]},{best_system_rmsd[2]}\n")


          
