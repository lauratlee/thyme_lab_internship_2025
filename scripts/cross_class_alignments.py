# parses through all alignment data for a given threshold distance and outputs a csv with the following columns:
  # human gene, zebrafish gene, match, mismatch, failed alignments, % similarity, best class, method
# where the highest match-mismatch ratio (and thereby highest % similarity) is recorded, along with the class and method (super/align) that produced those results

# example usage (run from thyme_lab_internship_2025/)
# python scripts/cross_class_alignments.py 2.0

