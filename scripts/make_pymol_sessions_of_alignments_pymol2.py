#note, this script was initially written by ChatGPT to write most of the pymol logic, and further modified by me (Laura) to enhance use of inputs

#usage: [UPDATE]
#python make_pymol_sessions_of_placements_pymol2.py /path/to/placements/directory/ ligandresidueindex list+string+of+residue+indices
#python make_pymol_sessions_of_placements_pymol2.py /scratch/abgvg9/discovery_results/top_1000_placement/agonist_12M_passing_placements/0 282 227+86+253+63+257

#breakdown of steps:
#load in 4S0V (which includes ligand suv)
#enter directory of choice and go through files one by one
#temporarily store the gene name of the selected file. If there is a parens in the file name, the gene name is between the parens. Otherwise the gene name is the dir name
#for each file, run the "super" command to align both receptors by structure
#select all residues of ONLY the target gpcr (not the reference) within 5 angstroms of ligand suv; this is the binding pocket
#save the pocket as a .pdb file titled "[gene]_pocket.pdb".

#this uses pymol2, which may be better

print("script started.")

import os, sys, re
import pymol2

#define the directory containing your files as a command line argument
#directory = 'path/to/your/files'
directory = sys.argv[1]
os.chdir(directory)

#have the directory end with a backslash if it doesn't from the input
if directory.endswith("/") == False:
    directory = directory + "/"



#helper function that takes in a gpcr directory and runs necessary pymol commands
def pymol_runner(gpcr_dir): 
    os.chdir(gpcr_dir)
    for _, _, genes in os.walk(gpcr_dir):
        for gene in genes:
            #make note of gene name for exported pdb
            if "(" not in gene: name = gpcr_dir
            else: 
                match = re.search(r'\(([^)]+)\)', gene)
                if match: name = match.group(1)

            print(f"NAME: {name}")

            #initialize pymol2 in headless mode
            with pymol2.PyMOL() as pymol:
                #set the internal gui width
                pymol.cmd.set('internal_gui_width', 600)

                #load in reference (4S0V)
                pymol.cmd.fetch("4S0V", "ref")
                print("loaded 4S0V")
                
                #load in target gene
                pymol.cmd.load(gene, "target")
                target_atoms = pymol.cmd.count_atoms("target")
                if target_atoms == 0:
                    print(f"[WARNING] Failed to load target structure: {gene}")

                #use the "super" command to align by structure
                pymol.cmd.super("ref", "target")
                print("structures aligned.")

                #select reference ligand
                pymol.cmd.select("ligand", "resn suv")
                ligand_atoms = pymol.cmd.count_atoms("ligand")
                if ligand_atoms == 0:
                    print(f"[WARNING] No atoms found for ligand in {gene}")

                #locate and select pocket residues of target
                pocket_sele = "byres (target within 5 of ligand) and target"
                pymol.cmd.select("pocket", pocket_sele)
                print(f"Atom count in pocket: {pymol.cmd.count_atoms('pocket')}")

                #save pocket
                output_dir = os.path.join("../../gpcr_pocket_dir", os.path.basename(gpcr_dir))
                os.makedirs(output_dir, exist_ok=True)
                output_path = os.path.join(output_dir, f"{name}_pocket.pdb")
                pymol.cmd.save(output_path, "pocket")
                print(f"saved ... {gpcr_dir} complete")
    os.chdir("..")
                




#walk through gpcr directory and run pymol_runner on each gpcr subdir
for sub in os.listdir(directory):
    sub_path = os.path.join(directory, sub)
    if os.path.isdir(sub_path):
        pymol_runner(sub_path)
