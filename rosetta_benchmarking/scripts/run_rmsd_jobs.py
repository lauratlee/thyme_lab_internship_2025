#run from system_dir_h_bonds or system_dir_close_res

import os

for file in os.listdir(os.getcwd()):
  if "_rmsd.slurm" in file or "_rmsd_close_res.slurm" in file or "ddg_rmsd.slurm" in file or "ddg_rmsd_close_res.slurm" in file:
    os.system(f"bsub {file}")

