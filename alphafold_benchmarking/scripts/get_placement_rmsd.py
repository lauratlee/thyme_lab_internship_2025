#the purpose of this script is to derive the rmsd of the rosetta placement rmsd from native
#this will create csv files that note the closest placement to native for the closest of all placements
#Since the placements are not aligned with the dude library (due to being made by rosetta), we need to align them. The ligand atom indices also do not match, so we will use a heuristic method from rdkit to get rmsd.
#run from system_dir_h_bonds/ (or system_dir_close_res/)

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

	#iterate through systems
	for system in os.walk(os.getcwd()):
		#temporary filter to just test on 9HZ0. delete second condition when doing all systems
		if os.path.isdir(system) and system == "9HZ0":
			print(system)

			#enter system directory
			os.chdir(system)

			#prepare original ligand for rdkit using openbabel
			os.system("obabel ligand.mol2 -O ligand.sdf --gen3d")	

			#open a system-specific file to write pairings of the files with rmsd 
			system_file = open(f"{system}_placements_summary.csv", 'w')
			#write header
			system_file.write("residue,file,rmsd\n")

			#declare placeholder variables to hold the best placement
			best_rmsd_1 = ["X","X","X"]

			#locate original reference file 
			orig_file = f"{system}.pdb"

			#load original file in pymol
			cmd.load(f"{orig_file}", "reference")

			#ensure reference was loaded properly
			num_ref_atoms = cmd.count_atoms("reference")
			if num_ref_atoms == 0:
				print("WARNING: no atoms in reference. Exiting.")
				sys.exit(1)
			else:
				print(f"ATOMS IN REFERENCE: {num_ref_atoms}")


			#create a dictionary the holds the placement files (and the residue they were derived from) and the corresponding rmsd values
			#key is a tuple of (residue, file), and the value is the rmsd
			placements_data = {}

			#iterate over the placements by residue folder by creating a list of folders to look at per system
			residue_list = []
			for folder in os.listdir(os.getcwd()):
				if os.path.isdir(folder) and "res_" in folder:
					residue_list.append(folder)


			#iterate over residues for analysis
			for residue in residue_list:
				print(residue)

				#iterate through folders for groups and add to a list
				group_list = []

				for res_folder in os.listdir(residue):
					if os.path.isdir(os.path.join("residue", "res_folder"):
						try:
							int(res_folder) #group folder names are all integers
							group_list.append(res_folder)
						except ValueError:
							continue

				#iterate through each group
				for group in group_list:
					#construct path to group folder
					group_path = os.path.join("residue", "group")

					#iterate through files in each group folder
					for group_file in os.listdir(group_path):
						if ".pdb" in group_file: #confirm file is a placement
							#load placement into pymol
							cmd.load(os.path.join(group_path, group_file), "placement")

							#ensure placement was loaded properly
							num_pla_atoms = cmd.count_atoms("placement")
							if num_pla_atoms == 0:
								print("WARNING: no atoms in placement. Exiting.")
								sys.exit(1)
							else:
								print(f"ATOMS IN PLACEMENT: {num_pla_atoms}")

							#align placement to reference
							cmd.align("placement", "reference")

							#select aligned ligand
							cmd.select("aligned_lig", "placement and not polymer.protein")

							#construct save name for aligned ligand and save
							aligned_lig_basename = group_file.split(".")[0] + "_aligned_lig.pdb"
							
							cmd.save(os.path.join(group_path, aligned_lig_basename), "aligned_lig")

							#clear the aligned ligand and placement from session but keep reference 
							cmd.delete("aligned_lig")
							cmd.delete("placement")


							#convert aligned ligand .pdb into .sdf format for rdkit
							aligned_lig_sdf_basename = group_file.split(".")[0] + "_aligned_lig.sdf"
							os.system(f"obabel {group_path}/{aligned_lig_basename} -O {group_path}/{aligned_lig_sdf_basename} --gen3d")

							#read reference ligand into rdkit without hydrogens
							ref_ligand = Chem.MolFromMolFile("ligand.sdf", removeHs=True)

							#check that reference loaded successfully
							if ref_ligand is None:
								print("WARNING: reference ligand did not read into rdkit. Exiting.")
								sys.exit(1)


							#read placement ligand into rdkit without hydrogens
							pla_ligand = Chem.MolFromMolFile(f"{group_path}/{aligned_lig_sdf_basename}", removeHs=True)

							#check that placement loaded successfully
							if pla_ligand is None:
								print("WARNING: placement ligand did not read into rdkit. Exiting.")
								sys.exit(1)


							

							
							
					
				
				
			
			
			








	

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
			
			#get the original placement and open it in pymol
			cmd.load(f"{dire}.pdb", "reference")
			#cmd.load("test_params/chain_A.pdb", "reference")

			
			if cmd.count_atoms("reference") == 0:
				print("WARNING: no reference loaded")
				sys.exit(1)
			print(f"ATOMS IN REFERENCE: {cmd.count_atoms("reference")}")

			#create a dictionary the holds the placement files and the corresponding confidence and rmsd values
			#the file is the key and the value is a 2 entry list of confidence then rmsd
			placements_data = {}


			#iterate over the placements by residue folder by creating a list of folders to look at per system
			residue_list = []
			for folder in os.listdir(os.getcwd()):
				if "res_" in folder:
					residue_list.append(folder)

			for residue in residue_list:
				print(residue)
				#iterate through group folders
				for folder in os.listdir(residue):
					if os.path.isdir(f"{residue}/{folder}/"):
						for file in os.listdir(f"{residue}/{folder}/"):
							if ".pdb" in file:
								#load placement into pymol
								cmd.load(f"{residue}/{folder}/{file}", "placement")
								if cmd.count_atoms("placement") == 0:
									print("WARNING: no placement loaded")
									sys.exit(1)
								print(f"ATOMS IN PLACEMENT: {cmd.count_atoms("placement")}")
		
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
		
								
		
								ref_ligand = Chem.MolFromMol2File("ligand.mol2", removeHs=True, sanitize=False)
								
								#check for successful load of ref_ligand and exit otherwise
								if ref_ligand is None:
									print("ERROR: Failed to load ligand.mol2 as RDKit molecule!")
									sys.exit(1)
								
		
								placement_ligand = Chem.MolFromPDBFile(f"{residue}/{file_basename}_aligned_lig.pdb", removeHs=True, sanitize=False)
		
								#check for successful load of placement_ligand and exit otherwise
								if placement_ligand is None:
									print("ERROR: Failed to load placement ligand as RDKit molecule!")
									sys.exit(1)
		
								"""try:
									Chem.SanitizeMol(placement_ligand, sanitizeOps=SanitizeFlags.SANITIZE_ALL ^ SanitizeFlags.SANITIZE_PROPERTIES)
								except Exception as e:
									print("Sanitization failed:", e)"""
		
		
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
			
			


							


		
