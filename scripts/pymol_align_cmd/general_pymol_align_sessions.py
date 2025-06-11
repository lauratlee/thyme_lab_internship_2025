#note, this script was initially written by ChatGPT to write most of the pymol logic, and further modified by me (Laura) to enhance use of inputs

#python (path to script) (REFERENCE FILE)
#this is different than the 4S0V-only script, which was run from the scripts folder. this should be run from the gpcr_dir folder

#breakdown of steps:
#load in reference (with ligand)
#enter directory of choice and go through files one by one
#temporarily store the gene name of the selected file. If there is a parens in the file name, the gene name is between the parens. Otherwise the gene name is the dir name
#for each file, run the "align" command to align both receptors by structure
#select all residues of ONLY the target gpcr (not the reference) within 5 angstroms of ligand; this is the binding pocket
#save the pocket as a .pdb file titled "[gene]_pocket.pdb".

#this uses pymol2, which may be better

import os, sys, re
import pymol2


#-----------------UPDATE------------------------

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
    CLASS: {gpcr_class}
    LIGAND: {ligand_name}
    """)
else:
    print("ERROR: Please provide a reference file that is within the gpcr_class_reps directory.")
    sys.exit(1)
    

#helper function that takes in a gpcr directory and runs necessary pymol commands
def pymol_runner(gpcr_dir): 
    print(f"Entered pymol_runner for {gpcr_dir}")
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

            #initialize pymol2 in headless mode
            with pymol2.PyMOL() as pymol:
                #set the internal gui width
                pymol.cmd.set('internal_gui_width', 600)

                #load in reference
                pymol.cmd.load(os.path.join("../../gpcr_class_reps", ref_file), "ref")
                if pymol.cmd.count_atoms("ref") != 0:
                    print("ref loaded")
                else:
                    print("[WARNING] Failed to load ref structure")
               
                
                #load in target gene
                pymol.cmd.load(gene, "target")
                target_atoms = pymol.cmd.count_atoms("target")
                if target_atoms == 0:
                    print(f"[WARNING] Failed to load target structure: {gene}")

                #use the "align" command to align by sequence
                pymol.cmd.align("target", "ref")
                print("structures aligned.")
                
                
                # ----- TEMPORARY DEBUG -------
                
                # Create a combined selection of the aligned reference and target structures
                combined_name = f"[align]{name}_aligned"
                output_dir = os.path.join("../../gpcr_pocket_dir", os.path.basename(gpcr_dir), gpcr_class)
                os.makedirs(output_dir, exist_ok=True)
                aligned_output_path = os.path.join(output_dir, f"{combined_name}.pdb")

                # Save both aligned structures (ref and target) into a single file
                pymol.cmd.save(aligned_output_path, "ref target")
                print(f"Saved aligned structure to: {aligned_output_path}")
                
                # ------ END OF DEBUG ----------
                

                #select reference ligand
                pymol.cmd.select("ligand", f"resn {ligand_name}")
                ligand_atoms = pymol.cmd.count_atoms("ligand")
                if ligand_atoms == 0:
                    print(f"[WARNING] No atoms found for ligand in {gene}")

                #locate and select pocket residues of target
                pocket_sele = "byres (target within 5 of ligand) and target"
                pymol.cmd.select("pocket", pocket_sele)
                print(f"Atom count in pocket: {pymol.cmd.count_atoms('pocket')}")

                #save pocket
                output_dir = os.path.join("../../gpcr_pocket_dir", os.path.basename(gpcr_dir), gpcr_class)
                os.makedirs(output_dir, exist_ok=True)
                output_path = os.path.join(output_dir, f"[align]{name}_pocket.pdb")
                pymol.cmd.save(output_path, "pocket")
                print(f"saved ... {gpcr_dir} complete")
    os.chdir("..")
                



def main():
    #walk through gpcr directory and run pymol_runner on each gpcr subdir
    for sub in os.listdir("."):
        print(f"Calling pymol_runner on: {sub}")
        pymol_runner(sub)

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
