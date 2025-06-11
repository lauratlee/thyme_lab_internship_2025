# debug cleaning file

import os, sys

folder = sys.argv[1]

for gpcr_dir in os.listdir("."):
  os.chdir(gpcr_dir)
  os.chdir(folder)
  for file in os.listdir("."):
    if file.endswith("centers.txt"):
      os.remove(file)
      print(f"{file} removed")
  os.chdir("../..")
