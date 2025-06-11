# initially generated using ChatGPT and then adjusted accordingly.
# parses alignment data and outputs a single csv file containing summary data for all relevant comparisons (human to zebrafish)
# example usage: python ../scripts/alignment_analysis.py Class_A 2.0
# this script should be run from gpcr_pocket_dir/

import csv, os, sys, re

def parse_residue_identity(label):
    """Extract residue identity (e.g., 'GLY') from 'GLY 53 (chain A)'."""
    if label == "NO ALIGNED RESIDUE":
        return None
    parts = label.strip().split()
    if len(parts) == 0:
        return None
    return parts[0]  # return just the residue identity

def load_alignments(csv_path):
    with open(csv_path, newline='') as f:
        reader = csv.reader(f)
        next(reader)  # skip header
        return list(reader)

def count_matches_by_identity(data):
    matched = 0
    mismatched = 0
    failed = 0

    for row in data:
        ref_id = parse_residue_identity(row[0])
        aligned_id = parse_residue_identity(row[1])

        if aligned_id is None:
            failed += 1
        elif ref_id == aligned_id:
            matched += 1
        else:
            mismatched += 1

    return matched, mismatched, failed

def calculate_percent_similarity(matched, mismatched):
    total = matched + mismatched
    return (matched / total * 100) if total > 0 else 0.0

def extract_genes_from_filename(filename):
    # remove extension
    base = os.path.splitext(filename)[0]  
    
    if class_name == "none":
        name_part = re.sub(r'_\d+(\.\d+)?$', '', base)
    else:
        # Step 1: Remove prefix
        name_part = base.replace("[align]", "")        # "ACKR3-ackr3a_2.0"
        
        # Step 2: Remove the trailing threshold (e.g., _2.0)
        name_part = re.sub(r'_\d+(\.\d+)?$', '', name_part)  # "ACKR3-ackr3a"
    
    # Step 3: Extract gene names using regex
    pattern = r'[a-zA-Z0-9]+(?:-[0-9]+)?'
    genes = re.findall(pattern, name_part)
    
    if len(genes) >= 2:
        return genes[0], genes[1]
    return "Unknown", "Unknown"

def alignment_summary(csv_path):
    data = load_alignments(csv_path)
    matched, mismatched, failed = count_matches_by_identity(data)
    similarity = calculate_percent_similarity(matched, mismatched)
    A, B = extract_genes_from_filename(csv_path)
    return [A, B, matched, mismatched, failed, round(similarity, 2)]

class_name = sys.argv[1]
threshold = sys.argv[2]

results = []

for gpcr in os.listdir("."):
    os.chdir(gpcr)
    
    if class_name == "none":
        os.chdir("cmd_align")
        for f in os.listdir("."):
            if f.endswith(".csv") and f != f.lower() and threshold in f:
                csv_path = f
                results.append(alignment_summary(csv_path))
    else:
        os.chdir(class_name)
        for f in os.listdir("."):
            if f.endswith(".csv") and f != f.lower() and threshold in f and "[align]" in f:
                csv_path = f
                results.append(alignment_summary(csv_path))
    
    os.chdir("../..")
    print(f"{gpcr} done")

#write output file
if class_name == "none":
    output_path = f"../align_outputs/human_zebrafish_alignment_summary_{threshold}.csv"
else:
    output_path = f"../align_outputs/[{class_name}]human_zebrafish_alignment_summary_{threshold}.csv"

with open(output_path, "w", newline='') as out_f:
    writer = csv.writer(out_f)
    writer.writerow(["human gene", "zebrafish gene", "matched", "mismatched", "failed alignments", "% similarity"])
    writer.writerows(results)

print("alignment summary done")

