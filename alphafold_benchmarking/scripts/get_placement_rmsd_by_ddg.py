#the purpose of this script is to derive the rmsd of the rosetta placement rmsd from native
#this will create csv files that note the closest placement to native for the closest of all placements
#Since the placements are not aligned with the dude library (due to being made by rosetta), we need to align them. The ligand atom indices also do not match, so we will use a heuristic method from rdkit to get rmsd.
#run from system_dir_h_bonds/ (or system_dir_close_res/)

import os,sys
import pymol2
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
#import numpy as np
#from openbabel import openbabel as ob
from openbabel import pybel
#import re
from rdkit.Chem import SanitizeFlags
from rdkit.Chem import AllChem
from pymol import CmdException


#helper function to strip all bond orders by setting everything to single bonds
def strip_bond_orders(mol):
	rw_mol = Chem.RWMol(mol)
	for bond in rw_mol.GetBonds():
		bond.SetBondType(Chem.rdchem.BondType.SINGLE)
	Chem.SanitizeMol(rw_mol)
	return rw_mol.GetMol()

#helper function to get ddg value from a placement .pdb file
def get_ddg(file):
  #open pdb and parse
  with open(file, 'r') as placement:
    for line in placement:
      if "Scoring: Post-HighResDock system ddG:" in line:
        parts = line.strip().split()
        ddg_value = float(parts[-1])
        break
  
  return ddg_value


#note system to work on based on user argument
target_system = sys.argv[1]

#begin a pymol session
with pymol2.PyMOL() as pymol:
	cmd = pymol.cmd

	#iterate through systems
	for system in os.listdir(os.getcwd()):
    
		if os.path.isdir(system) and system == target_system:
			print(system)

			#enter system directory
			os.chdir(system)

			# ---------------------------------- PREPARE ORIGINAL LIGAND ---------------------------------- #

			#prepare original ligand for rdkit using openbabel
			os.system("obabel ligand.mol2 -O ligand.sdf")	
			
			'''#read reference ligand into rdkit without hydrogens
			ref_ligand = Chem.MolFromMolFile("ligand.sdf", removeHs=True)
			Chem.SanitizeMol(ref_ligand)'''

			#read reference ligand into rdkit without hydrogens
			ref_ligand = Chem.MolFromMolFile("ligand.sdf", removeHs=True)

			#strip bond orders from reference
			ref_ligand = strip_bond_orders(ref_ligand)

			#check that reference loaded successfully
			if ref_ligand is None:
				print("WARNING: reference ligand did not read into rdkit. Exiting.")
				sys.exit(1)
				
			#get smiles of reference ligand
			ref_smiles = Chem.MolToSmiles(ref_ligand)

			# ---------------------------------------------------------------------------------------------- #

			
			#open a system-specific file to write pairings of the files with rmsd 
			with open(f"{system}_placements_summary.csv", "w") as system_file:
				#write header
				system_file.write("residue,file,rmsd\n")

			#declare placeholder variables to hold the best placement
			best_rmsd_1 = ["X","X","X","X"]

			#locate original reference file 
			orig_file = f"{system}.pdb"

			#load original file in pymol
			cmd.load(f"{orig_file}", "reference")

			#ensure reference was loaded properly
			num_ref_atoms = cmd.count_atoms("reference")
			if num_ref_atoms == 0:
				print("WARNING: no atoms in reference. Exiting.")
				sys.exit(1)
			#else:
				#print(f"ATOMS IN REFERENCE: {num_ref_atoms}")


			#create a dictionary the holds the placement files (and the residue they were derived from) and the corresponding rmsd values
			#key is a tuple of (residue, file), and the value is the rmsd
			placements_data = {}

			#iterate over the placements by residue folder by creating a list of folders to look at per system
			residue_list = []
			for folder in os.listdir(os.getcwd()):
				if os.path.isdir(folder) and "res_" in folder:
					residue_list.append(folder)


      #initialize a list that will hold the files with the top 10 ddgs
      ddg_files = []

			#iterate over residues for analysis
			for residue in residue_list:
				print(residue)

				#iterate through folders for groups and add to a list
				group_list = []

				for res_folder in os.listdir(residue):
					if os.path.isdir(os.path.join(residue, res_folder)):
						try:
							int(res_folder) #group folder names are all integers
							group_list.append(res_folder)
						except ValueError:
							continue

				#iterate through each group
				for group in group_list:
					print(f"residue {residue}, group {group}")
					#construct path to group folder
					group_path = os.path.join(residue, group)

          
          for group_file in os.listdir(group_path):
            if ".pdb" in group_file:
              group_file_path = os.path.join(group_path, group_file)
              ddg = get_ddg(group_file)

              #append a tuple of (file path, ddg) to ddg_files
              ddg_files.append((group_file_path, ddg))


          #keep only the top 10 ddg files
          ddg_files[:] = sorted(ddg_files, key=lambda x: x[1], reverse=True)[:10]
              

					#iterate through each file in ddg_files
					for ddg_file in ddg_files::
						if ".pdb" in ddg_file: #confirm file is a placement
							#load placement into pymol
							cmd.load(os.path.join(group_path, group_file), "placement")

							#ensure placement was loaded properly
							num_pla_atoms = cmd.count_atoms("placement")
							if num_pla_atoms == 0:
								print("WARNING: no atoms in placement. Exiting.")
								sys.exit(1)
							#else:
								#print(f"ATOMS IN PLACEMENT: {num_pla_atoms}")

							#align placement to reference. if alignment fails then skip file
							try:
								cmd.align("placement", "reference")
							except CmdException:
								continue

							#select aligned ligand
							cmd.select("aligned_lig", "placement and not polymer.protein")

							#construct save name for aligned ligand and save
							aligned_lig_basename = group_file.split(".")[0] + "_aligned_lig.mol2"
							
							cmd.save(os.path.join(group_path, aligned_lig_basename), "aligned_lig")

							#clear the aligned ligand and placement from session but keep reference 
							cmd.delete("aligned_lig")
							cmd.delete("placement")

							# ---------------------------------- PREPARE PLACEMENT LIGAND ---------------------------------- #
							
							#convert aligned ligand .mol2 into .sdf format for rdkit
							aligned_lig_sdf_basename = group_file.split(".")[0] + "_aligned_lig.sdf"
							os.system(f"obabel {group_path}/{aligned_lig_basename} -O {group_path}/{aligned_lig_sdf_basename}")


							#read placement ligand into rdkit without hydrogens
							pla_ligand = Chem.MolFromMolFile(f"{group_path}/{aligned_lig_sdf_basename}", removeHs=True)

							#strip bond orders from placement
							pla_ligand = strip_bond_orders(pla_ligand)

							#check that placement loaded successfully
							if pla_ligand is None:
								print("WARNING: placement ligand did not read into rdkit. Exiting.")
								sys.exit(1)

							'''#sanitize placement ligand
							Chem.SanitizeMol(pla_ligand)'''

							#get smiles of placement ligand
							pla_smiles = Chem.MolToSmiles(pla_ligand)

							# ---------------------------------------------------------------------------------------------- #

							rmsd = "X"

							#use the get best RMS function to derive the rmsd
							if ref_ligand and pla_ligand:
								try:
									rmsd = rdMolAlign.GetBestRMS(ref_ligand, pla_ligand)
									print(f"{group_path}/{aligned_lig_sdf_basename}", rmsd)
								except RuntimeError as e:
									print("Alignment failed:", e)

							#if the rmsd is X, don't add it
							if str(rmsd) == "X":
								continue
	
							#store the rmsd in the dictionary by the residue and file name
							placements_data[(residue, group_file)] = ["X",rmsd]
	
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
			with open(f"{system}_placements_summary.csv", "w") as system_file:
				system_file.write(f"{best_rmsd_1[0]},{best_rmsd_1[1]},{best_rmsd_1[2]}\n")


			#clear reference from pymol session
			cmd.delete("reference")


			

							


		
