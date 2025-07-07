# given a .csv of residue alignment results, prints the # match, # mismatch, # failed, and # similarity.
# example usage (to be run from folder where .csv file is): python [path to script] [alignment .csv file] 

import csv, sys

csv_file = sys.argv[1]

with open(csv_file, newline='') as f:
  reader = csv.reader(f)
  next(reader)  # skip header
  data = list(reader)

matched = 0
mismatched = 0
failed = 0

for row in data:
  # debug check to make sure each row has enough columns
  if len(row) < 2:
    print(f"[WARNING] Row does not have enough columns: {row}")
    sys.exit(1)
  
  ref_id = parse_residue_identity(row[0])
  aligned_id = parse_residue_identity(row[1])

  if aligned_id is None:
    failed += 1
  elif ref_id == aligned_id:
    matched += 1
  else:
    mismatched += 1

total = matched + mismatched
similarity = (matched / total * 100) if total > 0 else 0.0

print(f"""
        MATCHED: {matched}\n
        MISMATCHED: {mismatched}\n
        FAILED: {failed}\n
        % SIMILARITY: {similarity}
        """)





