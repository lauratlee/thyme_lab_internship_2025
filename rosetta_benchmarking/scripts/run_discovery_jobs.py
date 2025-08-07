# run from specific system folder, e.g. 9HZ0/

import os, sys

base = os.path.basename(os.getcwd())

for dir in os.listdir("."):
  if dir.startswith("res_"):
    num = dir.split("_")[1]
    os.chdir(dir)
    job_name = f"{base}_{num}.slurm"
    print(job_name)
    if not os.path.exists(job_name):
      print("[WARNING] file does not exist")
      sys.exit(1)
    else:
      os.system(f"bsub {job_name}")
    os.chdir("..")
  
