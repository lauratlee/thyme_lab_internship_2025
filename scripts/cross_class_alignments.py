# parses through all alignment data for a given threshold distance and outputs a csv with the following columns:
  # human gene, zebrafish gene, match, mismatch, failed alignments, % similarity, best class, method
# where the highest match-mismatch ratio (and thereby highest % similarity) is recorded, along with the class and method (super/align) that produced those results

# example usage (run from thyme_lab_internship_2025/)
# python scripts/cross_class_alignments.py 2.0

# a temporary dictionary is created where the key is a tuple of (human gene, zebrafish gene), and the value is a list of tuples with the corresponding alignment data.
  # value tuple format is (match, mismatch, % similarity, class, method). The classes are A, B1, C, or F, or "A (4s0V)" (referring to the initial alignments done with 4s0v, a class A GPCR, as the reference)

import os, sys, csv, re, pprint

threshold = sys.argv[1]
data_dict = {}

# takes a .csv file along with the file data's method (super/align), and adds the file's data to the data_dict dictionary
def parse_summary(file, gpcr_class, method):
  with open(file, newline = '') as f:
    reader_f = csv.reader(f)
    next(reader_f)
    for row in reader_f:
      key = (row[0], row[1])
      value = tuple(row[2:] + [gpcr_class, method])
      if key in data_dict:
        data_dict[key].append(value)
      else:
        data_dict[key] = [value]

# extracts the gpcr class name from a given file name
def get_class(file):
  match = re.search(r'\[(.*?)\]', file)
  if match:
    return match.group(1)
  else:
    return "Class_A (4S0V)"

# parse align_outputs/
os.chdir("align_outputs")
method = "align"
for summary in os.listdir("."):
  if summary.endswith(f"{threshold}.csv"):
    gpcr_class = get_class(summary)
    parse_summary(summary, gpcr_class, method)
os.chdir("..")

# parse super_outputs/
os.chdir("super_outputs")
method = "super"
for summary in os.listdir("."):
  if summary.endswith(f"{threshold}.csv"):
    gpcr_class = get_class(summary)
    parse_summary(summary, gpcr_class, method)
os.chdir("..")

pprint.pprint(data_dict)


    
  
      
  
  




