#writes group_specific rmsd job files for a given system's residue
#run from a system directory, e.g. system_dir_close_res
#example usage: python /path/to/script/ 9M3N res_44

import os, sys

system = sys.argv[1]
residue = sys.argv[2]

group_list = []

for dire in os.listdir(f"{system}/{residue}"):
  if dire.isdigit():
    group_list.append(dire)

for group in group_list:
  with open(f"{system}/{residue}_{group}_rmsd.slurm", "w") as group_job_file:
    group_job_file.write("#!/bin/bash \n")
    group_job_file.write("#BSUB -n 1 \n")
    group_job_file.write("#BSUB -R 'span[hosts=1]' \n")
    group_job_file.write("#BSUB -W 1600 \n")
    group_job_file.write("#BSUB -R 'rusage[mem=100000]' \n")
    group_job_file.write(f"#BSUB -e {system}/{residue}/{residue}_{group}_rmsd_err_log.err \n")
    group_job_file.write(f"python ../scripts/get_placement_rmsd_by_group.py {system} {residue} {group} > {system}/{residue}/{residue}_{group}_rmsd_out.txt\n")


  print(f"{system}/{residue}_{group}_rmsd.slurm")
  os.system(f"bsub {system}/{residue}_{group}_rmsd.slurm")
