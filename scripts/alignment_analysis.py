import csv, sys, os

file_name = sys.argv[1]

def parse_residue_label(label):
    """Extract residue name and number from a label like 'GLY 53 (chain A)'"""
    if label == "NO ALIGNED RESIDUE":
        return None
    parts = label.split()
    if len(parts) < 2:
        return None
    resn = parts[0]
    resi = parts[1]
    return (resn, resi)

def load_alignments(csv_path):
    with open(csv_path, newline='') as f:
        reader = csv.reader(f)
        next(reader)  # skip header
        return list(reader)

def count_matches(data):
    matched = 0
    mismatched = 0
    failed = 0

    for row in data:
        ref = parse_residue_label(row[0])
        aligned = parse_residue_label(row[1])

        if aligned is None:
            failed += 1
        elif ref == aligned:
            matched += 1
        else:
            mismatched += 1

    return matched, mismatched, failed

def calculate_percent_similarity(matched, mismatched):
    total = matched + mismatched
    if total == 0:
        return 0.0
    return (matched / total) * 100

def analyze_alignment(csv_path):
    data = load_alignments(csv_path)
    matched, mismatched, failed = count_matches(data)
    similarity = calculate_percent_similarity(matched, mismatched)

    print(f"Results for {csv_path}:\n")
    print(f"Matched residues     : {matched}")
    print(f"Mismatched residues  : {mismatched}")
    print(f"Failed alignments    : {failed}")
    print(f"Percent similarity   : {similarity:.2f}%")

# Example usage
analyze_alignment(file_name)
