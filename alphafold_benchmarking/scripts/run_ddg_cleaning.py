#runs the bsub jobs for remove_large_ddgs.py
#run from a system_dir folder, e.g. system_dir_h_bonds/

import os

for r,d,f in os.walk(os.getcwd()):
  for dir in d:
    if dir != "9HZ0" and dir != "9I0T":
      os.chdir(dir)
      os.system(f"bsub {dir}_clean_ddgs.slurm")
      os.chdir("..")
