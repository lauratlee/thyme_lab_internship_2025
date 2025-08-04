#parses system library and returns number of rmsds 0-2A, 2-5A, and >5A from reference for the categories of All Placements, Top 10 DDG, Top 1 DDG
#run from system library, e.g. system_dir_h_bonds

import os

#helper function to update counts of ranges [0,2), [2,5], (5+)
def sort_rmsd(count0, count2, count5, num):
  if 0 <= num < 2: count0 += 1
  elif 2 <= num <= 5: count2 += 1
  elif num > 5: count5 += 1

  return(count0, count2, count5)

#initialize range counts by all placements, top 10 ddg, and top 1 ddg
count0_all, count2_all, count5_all = 0,0,0
count0_10, count2_10, count5_10 = 0,0,0
count0_1, count2_1, count5_1 = 0,0,0


#initialize system count
lib_length = 0

#note system library (h_bonds or close_res)
sys_library = os.path.basename(os.getcwd())

#iterate thru all systems in library
for system in os.listdir(os.getcwd()):
  if os.path.isdir(system):
    lib_length += 1
    os.chdir(system)

    '''#get rmsd from all placements
    with open(f"{system}_placements_summary.csv", 'r') as sum_all:
      entry_all = sum_all.readline()
      rmsd_all = entry_all.strip().split(',')[2]

      data_all = sort_rmsd(count0_all, count2_all, count5_all, rmsd_all)'''

    #get rmsd from top 10 ddg
    with open(f"{system}_placements_summary_ddg_10.csv", 'r') as sum_10:
      entry_10 = sum_10.readline()
      rmsd_10 = entry_10.strip().split(',')[3]

      data_10 = sort_rmsd(count0_10, count2_10, count5_10, rmsd_10)

    #get rmsd from top 10 ddg
      with open(f"{system}_placements_summary_ddg_1.csv", 'r') as sum_1:
        entry_1 = sum_1.readline()
        rmsd_1 = entry_1.strip().split(',')[3]

      data_1 = sort_rmsd(count0_1, count2_1, count5_1, rmsd_1)


data_all = (5,5,5)



print(f"""LIBRARY: {sys_library}\n
      # OF SYSTEMS: {lib_length}\n
      ALL PLACEMENTS: {data_all}\n
      TOP 10 DDG: {data_10}\n
      TOP 1 DDG: {data_1}\n
      """)
  











    


