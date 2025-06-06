# quick script to remove temporary alignment debug files

import os

for folder in os.listdir("."):
  os.chdir(f"{os.path.basename(folder)}")
  for file in os.listdir("."):
    if "_aligned.pdb" in file:
      os.remove(file)
  os.chdir("..")
