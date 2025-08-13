#this is the 2nd step after running residue-specific rmsd calculations within a system (mainly for 9M3N, 9QEK, 9QEL)
#run from system folder, e.g. 9M3N/

import os, sys, csv, ast

system = os.path.basename(os.getcwd())
best_rmsds = []

for file in os.listdir(os.getcwd()):
  if "res_" in file and "placements_summary" in file:
    with open(file, newline="") as rmsd_file:
      reader = list(csv.reader(rmsd_file))
      if len(reader) >= 2:
        #see if the second to last row has the best rmsd for that residue
        second_last_row = reader[-2]
        second_last_row_text = "".join(second_last_row)
        print(second_last_row_text)
        
        if "BEST RMSD 1 ENTRY:" in second_last_row_text:
          #get entry as a string
          entry_str = second_last_row_text.split(":", 1)[1].strip()
          
          #convert entry to a list
          entry_list = ast.literal_eval(entry_str)

          best_rmsds.append(entry_list)
        else:
          continue
      else:
        continue

#sort best_rmsds by rmsd value
best_rmsds.sort(key=lambda x: x[2])
best_system_rmsd = best_rmsds[0]

#write to system summary file
with open(f"{system}_placements_summary.csv", "w") as system_file:
  system_file.write("residue,file,rmsd\n")
  system_file.write(f"{best_system_rmsd[0]},{best_system_rmsd[1]},{best_system_rmsd[2]}\n")
  
          
