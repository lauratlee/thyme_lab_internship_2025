#makes directories for gpcr class-specific outputs

import os

for gpcr in os.listdir("."):
  os.chdir(gpcr)
  os.makedirs("Class_A")
  os.makedirs("Class_B1")
  os.makedirs("Class_C")
  os.makedirs("Class_F")
  os.chdir("..")
