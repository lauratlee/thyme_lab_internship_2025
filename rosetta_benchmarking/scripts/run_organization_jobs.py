#run from system library, e.g. system_dir_close_res/

import os

for system in os.listdir(os.getcwd()):
  if os.path.isdir(system) and "9" in system:
    print(system)
    os.chdir(system)
    os.system(f"bsub {system}_organization_close_res.slurm")
    os.chdir("..")
