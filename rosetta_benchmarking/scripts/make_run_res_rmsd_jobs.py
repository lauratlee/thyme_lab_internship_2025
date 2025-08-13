#writes residue_specific rmsd job files for a given system

import os, sys

target_system = sys.argv[1]

residue_list = []

for dire in os.listdir(target_system):
  if os.path.isdir(os.path.join(target_system, dire)) and dire.startswith("res_"):
    residue_list.append(dire)

for residue in residue_list:
  with open(f"{target_system}_{residue}_rmsd.slurm", "w") as res_job_file:
    res_job_file.write("#!/bin/bash \n")
    res_job_file.write("#BSUB -n 1 \n")
    res_job_file.write("#BSUB -R 'span[hosts=1]' \n")
    res_job_file.write("#BSUB -W 1600 \n")
    res_job_file.write("#BSUB -R 'rusage[mem=100000]' \n")
    res_job_file.write(f"#BSUB -e {target_system}/{residue}_rmsd_err_log.err \n")
    res_job_file.write(f"python ../scripts/get_placement_rmsd_by_res.py {target_system} {residue} > {target_system}/{target_system}_{residue}_rmsd_out.txt")


  print(f"{target_system}_{residue}_rmsd.slurm")
  os.sys(f"bsub {target_system}_{residue}_rmsd.slurm")








