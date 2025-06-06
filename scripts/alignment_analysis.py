import csv, sys, os

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

def analyze_alignment(csv_path):
    data = load_alignments(csv_path)
    matched, mismatched, failed = count_matches_by_identity(data)
    similarity = calculate_percent_similarity(matched, mismatched)

    print(f"Results for {csv_path}:\n")
    print(f"Matched residues     : {matched}")
    print(f"Mismatched residues  : {mismatched}")
    print(f"Failed alignments    : {failed}")
    print(f"Percent similarity   : {similarity:.2f}%")

csv_path = sys.argv[1]

analyze_alignment(csv_path)
