# quick prep script to copy over ligand.params files from the h_bonds dir to the close_res dir
# run from close_res dir

import os

for dir in os.listdir(os.getcwd):
  name = os.path.basename(dir)
  os.system(f"cp ../system_dir_h_bonds/{dir}/test_params/ligand.params {dir}/test_params/ligand.params"
  
