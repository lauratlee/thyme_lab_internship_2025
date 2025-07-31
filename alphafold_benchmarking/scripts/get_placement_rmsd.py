#the purpose of this script is to derive the rmsd of the rosetta placement rmsd from native
#this will create csv files that note the closest placement to native for the closest of all placements
#Since the placements are not aligned with the dude library (due to being made by rosetta), we need to align them. The ligand atom indices also do not match, so we will use a heuristic method from rdkit to get rmsd.
#run from system_dir_h_bonds/ (or system_dir_close_res/)

#import os,sys
import os,sys
import pymol2
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
import numpy as np
from openbabel import openbabel as ob
from openbabel import pybel
import re
from rdkit.Chem import SanitizeFlags

#begin a pymol session
with pymol2.PyMOL() as pymol:
	cmd = pymol.cmd

	#write csv file for best csvs and write a header line for the system, the placement file, and the rmsd
	best_1 = open("best_placements_1.csv", "w")
	best_1.write("system,placement,rmsd\n")


	#iterate over each system in the library
	for r,d,f in os.walk(os.getcwd()):
		for dire in d:
			print(dire)

			#temporary filter so we only test this on 9HZ0 until it is time to run on the full set
			#comment/delete when done testing!
			if dire != "9HZ0":
				continue

			#skip 9LE4, which had 2 identical ligands and was not ideal for rosetta docking
			if dire == "9LE4":
				continue

			#enter system directory
			os.chdir(dire)

			#open a file to write pairings of the files with rmsd
			#open it in the respective folder in the alphafold section of the repository
			system_file = open(f"{dire}_placements_summary.csv", "w")
			system_file.write("residue,file,rmsd\n")

			#declare placeholder variables to hold the best placement
			best_rmsd_1 = ["X","X","X"]
			
			#get the original placement from the dude library and open it in pymol
			cmd.load(f"{dire}.pdb", "reference")

			#create a dictionary the holds the placement files and the corresponding confidence and rmsd values
			#the file is the key and the value is a 2 entry list of confidence then rmsd
			placements_data = {}


			#iterate over the placements by residue folder by creating a list of folders to look at per system
			residue_list = []
			for folder in os.listdir(os.getcwd()):
				if "res_" in folder:
					residue_list.append(folder)

			for residue in residue_list:
				for file in os.listdir(residue):
					if ".pdb" in file:
						#load placement into pymol
						cmd.load(f"{residue}/{file}", "placement")

						#align placement to reference
						cmd.align("placement", "reference")

						#select the aligned ligand and save as .pdb 
						cmd.select("aligned_lig", "placement and not polymer.protein")
			
						#derive save name for ligand and save
						file_basename = file.split(".")[0]
						cmd.save(f"{residue}/{file_basename}_aligned_lig.pdb", "aligned_lig")

						#clear the aligned ligand and placement from session but keep reference
						cmd.delete("aligned_lig")
						cmd.delete("placement")

						#make a fixed version of the reference from the original so that the element is recognized
						old_ref_file = open("ligand.mol2", "r")
						fixed_ref_file = open("ligand_fixed.mol2", "w")

						#isolate section with atoms into a list 
						atom_section = []
						inside_atom_block = False
			
						for line in old_ref_file.readlines():
							if line.strip() == "@<TRIPOS>ATOM":
								inside_atom_block = True
								continue
							elif line.strip() == "@<TRIPOS>BOND":
								break
							elif inside_atom_block:
								atom_section.append(line.strip())

						for atom_line in atom_section:
							#check for and remove waters
							if "HOH" in atom_line:
								continue

							#get atom name
							atom_name = atom_line.split()[1]

							#derive element
							element = re.match(r"[A-Za-z]+", atom_name.strip()).group(0).capitalize()

							#skip hydrogens
							if element == "H":
								continue

							stripped_line = atom_line.rstrip("\n")

							fixed_line = stripped_line[:76].ljust(76) + element.rjust(2) 

							fixed_ref_file.write(fixed_line + "\n")

						#close streams
						old_ref_file.close()
						fixed_ref_file.close()

						ref_ligand = Chem.MolFromMol2File("ligand_fixed.mol2", removeHs=True, sanitize=False)

						placement_ligand = Chem.MolFromPDBFile(f"{residue}/{file_basename}_aligned_lig.pdb", removeHs=True, sanitize=False)

						try:
							Chem.SanitizeMol(placement_ligand, sanitizeOps=SanitizeFlags.SANITIZE_ALL ^ SanitizeFlags.SANITIZE_PROPERTIES)
						except Exception as e:
							print("Sanitization failed:", e)


						ref_smiles = Chem.MolToSmiles(ref_ligand)
						pla_smiles = Chem.MolToSmiles(placement_ligand)

						rmsd = "X"

						#use the get best RMS function to derive the rmsd
						if ref_ligand and placement_ligand:
							try:
								rmsd = rdMolAlign.GetBestRMS(ref_ligand, placement_ligand)
								print(f"{residue}/{file_basename}_aligned_lig.pdb", rmsd)
							except RuntimeError as e:
								print("Alignment failed:", e)


						#if the rmsd is X, don't add it
						if str(rmsd) == "X":
							continue

						#store the rmsd in the dictionary by the residue and file name
						placements_data[(residue, file)] = ["X",rmsd]

						#we are now done with the placement, and can move to the next

			#done with all placements for the system, correlate rmsd and confidence and update dictionaries and finish the system-specific csv

			placements_list = []

			for entry in placements_data.keys():
				placement_residue = entry[0]
				placement_file = entry[1]

				#add placements to a list
				placements_list.append([placement_residue, placement_file, float(placements_data[entry][1])])


			#sort placements_list by rmsd in ascending order
			sorted_list = sorted(placements_list, key=lambda x: x[2])

			#get the best (lowest) rmsd entry
			best_rmsd_1_entry = sorted_list[0]
			best_rmsd_1 = [best_rmsd_1_entry[0], best_rmsd_1_entry[1], best_rmsd_1_entry[2]]

			#write result to system file
			system_file.write(f"{best_rmsd_1[0]},{best_rmsd_1[1]},{best_rmsd_1[2]:.3f}\n")

			#exit system directory so we can write the best rmsd to the general best rmsds file
			os.chdir("..")

			best_1.write(f"{dire},{best_rmsd_1[0]},{best_rmsd_1[1]},{best_rmsd_1[2]:.3f}\n")

			#clear reference from pymol session
			cmd.delete("reference")

			#print message to let us know that analysis of current system is done
			print(f"{dire} DONE")
			
			


							


		
