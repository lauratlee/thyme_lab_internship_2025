# cleaning out badly named txt files

import os

for gpcr in os.listdir("."):
  os.chdir(gpcr)
  os.chdir("cmd_align")
  for filename in os.listdir("."):
    if filename.count(".txt") == 2:
        print(f"Deleting file: {filename}")
        os.remove(filename)
  os.chdir("../..")
