# creates the empty patches.txt and exclude_pdb_component_list.txt files in all test_params directories
# run from system_dir/

import os

for _, dirs, _ in os.walk(os.getcwd()):
  for dir in dirs:
    os.chdir(f"{dir}/test_params/")
    os.system("touch patches.txt")
    os.system("touch exclude_pdb_component_list.txt")
    os.chdir("../..")


print("test_params_prep done")
