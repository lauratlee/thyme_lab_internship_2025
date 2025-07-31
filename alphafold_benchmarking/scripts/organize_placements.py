#organizes placement files into groups of 1000
#run from a system folder, e.g. 9HZ0/

import os
import math

root_dir = os.getcwd()

for dir1 in os.listdir(root_dir):
    #filter for residue folders only
    if "res_" in dir1 and os.path.isdir(os.path.join(root_dir, dir1)):
        res_path = os.path.join(root_dir, dir1)

        #get .pdb placement files
        pdb_files = [f for f in os.listdir(res_path) if f.endswith(".pdb")]

        count = len(pdb_files)
        print(f"{count} placements in {dir1}")

        #calculate number of 1000-size groups
        groups = math.ceil(count / 1000)

        # create group folders
        for i in range(groups):
            os.makedirs(os.path.join(res_path, str(i)), exist_ok=True)

        # move files into correct folders
        for file in pdb_files:
            try:
                index = int(file.split("_")[-1].split(".")[0])
                group_id = index // 1000
                src = os.path.join(res_path, file)
                dst = os.path.join(res_path, str(group_id), file)
                os.rename(src, dst)
            except Exception as e:
                print(f"Error processing {file} in {dir1}: {e}")

        print(f"{dir1} DONE")

print("All residues in system done")
