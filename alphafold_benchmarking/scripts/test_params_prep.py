# creates the empty patches.txt and exclude_pdb_component_list.txt files in all test_params directories
# run from system_dir/

import os

for _, system, _ in os.getcwd():
  os.chdir(f"{system}/test_params/")
  os.system("touch patches.txt")
  os.system("touch exclude_pdb_component_list.txt")
  os.chdir("../..")


print("test_params_prep done")
