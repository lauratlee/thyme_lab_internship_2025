# parses through all alignment data for a given threshold distance and outputs a csv with the following columns:
  # human gene, zebrafish gene, match, mismatch, failed alignments, % similarity, best class, method
# where the highest match-mismatch ratio (and thereby highest % similarity) is recorded, along with the class and method (super/align) that produced those results

# example usage (run from thyme_lab_internship_2025/)
# python scripts/cross_class_alignments.py 2.0

# a temporary dictionary is created where the key is a tuple of (human gene, zebrafish gene), and the value is a list of tuples with the corresponding alignment data.
  # value tuple format is (match, mismatch, % similarity, class, method). The classes are A, B1, C, or F, or "A (4s0V)" (referring to the initial alignments done with 4s0v, a class A GPCR, as the reference)

import os, sys




