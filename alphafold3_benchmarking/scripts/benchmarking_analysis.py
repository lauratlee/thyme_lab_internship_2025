import csv

#note systems to skip; default to empty list if no input
systems_to_skip = sys.argv[1].split(',') if len(sys.argv) > 1 else []

summary_csv = "../system_dir/best_placements_summary.csv"

#helper function to update counts of ranges [0,2), [2,5], (5+)
def sort_rmsd(count0, count2, count5, num):
  if 0 <= num < 2: count0 += 1
  elif 2 <= num <= 5: count2 += 1
  elif num > 5: count5 += 1

  return(count0, count2, count5)

#initialize counts
system_counter = 0
count0_all, count2_all, count5_all = 0,0,0
count0_10, count2_10, count5_10 = 0,0,0
count0_1, count2_1, count5_1 = 0,0,0

with open(summary_csv, newline='') as summary:
  reader = csv.reader(summary)
  header = next(reader)

  for row in reader:
    #separate the system name and the system's data
    system = row[0]
    values = row[1:]

    #skip systems that do not have data or were not included in rosetta analysis
    if values == ["X","X","X"] or system in systems_to_skip:
      continue

    system_counter += 1

    rmsd_all = float(values[0])
    rmsd_10 = float(values[1])
    rmsd_1 = float(values[2])

    print(rmsd_all, rmsd_10, rmsd_1)
    
    count0_all, count2_all, count5_all = sort_rmsd(count0_all, count2_all, count5_all, rmsd_all)
    count0_10, count2_10, count5_10 = sort_rmsd(count0_10, count2_10, count5_10, rmsd_10)
    count0_1, count2_1, count5_1 = sort_rmsd(count0_1, count2_1, count5_1, rmsd_1)


data_all = (count0_all, count2_all, count5_all)
data_10 = (count0_10, count2_10, count5_10)
data_1 = (count0_1, count2_1, count5_1)



print(f"""# OF SYSTEMS: {system_counter}\n
      ALL PLACEMENTS: {data_all}\n
      TOP 10 DDG: {data_10}\n
      TOP 1 DDG: {data_1}\n
      """)
    
