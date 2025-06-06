#makes directories to hold temporary debug files

import os

for gpcr in os.listdir("."):
  os.chdir(gpcr)
  os.makedirs("debug_files")
  os.chdir("cmd_align")
  os.makedirs("debug_files")
  os.chdir("../..")
