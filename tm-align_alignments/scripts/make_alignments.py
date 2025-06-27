# uses tm-align to create alignments against a given reference file.
# example usage: (to be run from gpcr_dirs/) python ../scripts/make_alignments.py ../../gpcr_class_reps/2rh1_chainA.pdb

import sys, os, re


# ensure that given reference file is valid
ref_map = {
    "2rh1_chainA.pdb": ("Class_A", "cau"),
    "4k5y_chainA.pdb": ("Class_B1", "1q5"),
    "7m3g_chainA.pdb": ("Class_C", "h43"),
    "4jkv_chainA.pdb": ("Class_F", "1ks")
}

ref_file = os.path.basename(sys.argv[1])

if ref_file in ref_map:
    gpcr_class, ligand_name = ref_map[ref_file]
    print(f"""
    FILE: {ref_file}
    CLASS: {gpcr_class}
    LIGAND: {ligand_name}
    """)
else:
    print("ERROR: Please provide a reference file that is within the gpcr_class_reps directory.")
    sys.exit(1)

# helper function that takes in a gpcr directory and runs necessary tm-align commands
def tm_align_runner(gpcr_dir):
    os.chdir(gpcr_dir)
    for _, _, genes in os.walk("."):
        for gene in genes:
            if not gene.endswith(".pdb"):
                continue
            print(f"Found gene file: {gene}")
            #make note of gene name for exported pdb
            if "(" not in gene: name = os.path.basename(gpcr_dir)
            else: 
                match = re.search(r'\(([^)]+)\)', gene)
                if match: name = match.group(1)

            print(f"NAME: {name}")

            os.system(f"~/TMalign '{gene}' {sys.argv[1]} -o ~/thyme_lab_internship_2025/tm-align_alignments/gpcr_pocket_dir/{gpcr_dir}/{gpcr_class}/{name}")

            print(f"saved ... {gene} complete")

    os.chdir("..")


def main():
    # walk through gpcr directory and run pymol_runner on each gpcr subdir
    for sub in os.listdir("."):
        print(f"Calling tm_align_runner on: {sub}")
        tm_align_runner(sub)

while True:
    answer = input("Continue with alignments? [y/n]").strip().lower()
    if answer == "y":
        main()
        break
    elif answer == "n":
        print("Exiting program.")
        sys.exit(0)
    else:
        print("Invalid input, please answer y/n")

