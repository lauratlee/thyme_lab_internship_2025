# initially generated using ChatGPT and then adjusted accordingly.
# parses alignment data and outputs a single csv file containing summary data for all relevant comparisons (human to zebrafish)
# example usage: python ../scripts/alignment_analysis.py Class_A 2.0
# this script should be run from gpcr_pocket_dir/

import csv, os, sys, re

#dictionary for pairwise comparison of residues for similarity
similar_residues_3letter = {
    # Hydrophobic (aliphatic + sulfur-containing)
    'ALA': {'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'SER', 'THR', 'GLY', 'CYS'},
    'VAL': {'ALA', 'VAL', 'LEU', 'ILE', 'MET'},
    'LEU': {'ALA', 'VAL', 'LEU', 'ILE', 'MET'},
    'ILE': {'ALA', 'VAL', 'LEU', 'ILE', 'MET'},
    'MET': {'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'CYS'},
    # Aromatic
    'PHE': {'PHE', 'TYR', 'TRP'},
    'TYR': {'PHE', 'TYR', 'TRP'},
    'TRP': {'PHE', 'TYR', 'TRP'},
    # Small / flexible / slightly polar
    'GLY': {'GLY', 'ALA', 'SER', 'THR'},
    'SER': {'SER', 'THR', 'ALA', 'GLY', 'ASN', 'GLN', 'CYS'},
    'THR': {'THR', 'SER', 'ALA', 'GLY', 'ASN', 'GLN'},
    # Polar uncharged (amide or hydroxyl)
    'ASN': {'ASN', 'GLN', 'SER', 'THR'},
    'GLN': {'GLN', 'ASN', 'SER', 'THR'},
    # Positively charged (basic)
    'LYS': {'LYS', 'ARG', 'HIS'},
    'ARG': {'ARG', 'LYS', 'HIS'},
    'HIS': {'HIS', 'LYS', 'ARG'},
    # Negatively charged (acidic)
    'ASP': {'ASP', 'GLU'},
    'GLU': {'GLU', 'ASP'},
    # Sulfur-containing (and overlaps with small/polar)
    'CYS': {'CYS', 'SER', 'ALA', 'MET'},
    # Special case: Proline (rigid, unique)
    'PRO': {'PRO'},
}


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
    similar = 0

    for row in data:
        # debug check to make sure each row has enough columns
        if len(row) < 2:
            print(f"[WARNING] Row does not have enough columns: {row}")
            sys.exit(1)
        
        ref_id = parse_residue_identity(row[0])
        aligned_id = parse_residue_identity(row[1])

        #debug print
        print(aligned_id,ref_id)

        if aligned_id is None:
            failed += 1
        elif ref_id == aligned_id:
            matched += 1
            similar += 1
        elif aligned_id in similar_residues_3letter.get(ref_id, set()):
            similar += 1
        else:
            mismatched += 1

    return matched, similar, mismatched, failed

def calculate_percent_identity(matched, mismatched):
    total = matched + mismatched
    return (matched / total * 100) if total > 0 else 0.0

def calculate_percent_similarity(matched, similar, mismatched):
    total = matched + mismatched
    return (similar / total * 100) if total > 0 else 0.0

def extract_genes_from_filename(filename, class_name):
    # remove extension
    base = os.path.splitext(os.path.basename(filename))[0]  
    
    # Step 1: Remove the trailing threshold (e.g., _2.0)
    name_part = re.sub(r'_\d+(\.\d+)?$', '', base)  # "ACKR3-ackr3a"
    
    # Step 2: Extract gene names using regex
    genes = name_part.split("-", 1)  # ['HCAR1', 'hcar1-2']

    if len(genes) >= 2:
        return genes[0], genes[1]
    else:
      print("[WARNING] Could not find 2 gene names. Exiting program.")
      sys.exit(1)

def alignment_summary(csv_path):
    data = load_alignments(csv_path)
    A, B = extract_genes_from_filename(csv_path, class_name)
    print(f"human gene: {A}\nzebrafish gene: {B}")
    matched, similar, mismatched, failed = count_matches_by_identity(data)
    identity = calculate_percent_identity(matched, mismatched)
    similarity = calculate_percent_similarity(matched, similar, mismatched)
    return [A, B, matched, similar, mismatched, failed, round(identity, 2), round(similarity, 2)]

class_name = sys.argv[1]
threshold = sys.argv[2]

results = []

for gpcr in os.listdir("."):
    for file in os.listdir(os.path.join(gpcr, class_name)):
        if file.endswith(f"_{threshold}.csv") and file[0].isupper():
            csv_path = os.path.join(gpcr, class_name, file)
            results.append(alignment_summary(csv_path))
    print(f"{gpcr} done")
    
    
#write output file
output_path = f"../outputs/[{class_name}]human_zebrafish_alignment_summary_{threshold}_with_similarity.csv"

with open(output_path, "w", newline='') as out_f:
    writer = csv.writer(out_f)
    writer.writerow(["human gene", "zebrafish gene", "matched", "similarity", "mismatched", "failed alignments", "% identity", "% similarity"])
    writer.writerows(results)

print("alignment summary done")
