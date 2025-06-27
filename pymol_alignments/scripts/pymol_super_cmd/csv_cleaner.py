# quick script to remove csv files

import os

for folder in os.listdir("."):
  os.chdir(f"{os.path.basename(folder)}")
  for file in os.listdir("."):
    if file.endswith(".csv"):
      os.remove(file)
  os.chdir("..")
