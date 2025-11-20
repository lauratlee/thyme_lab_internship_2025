#the purpose of this script is to go system by system and determine the best RMSDS for the best, top 10, and all placements per system using the 250 conformer schema
#this script is to be run in the scripts directory like: python ari_get_rmsds_rosetta.py

#import
import os,sys

#create a dictionary that holds the systems and a list for the best rmsds per system
systems_rmsds = {}

#store the workign location, which should be thyme_lab_internship_2025/rosetta_benchmarking/scripts/
starting_script_dir = os.getcwd()

print("starting script directory, which should contain thyme_lab_internship_2025/rosetta_benchmarking/scripts", starting_script_dir)

#move to the system directory
os.chdir("../system_dir/")

#store the system directory
system_dir = os.getcwd()

#iterate over the systems directories and go system by system, collecting the best rmsds
for r,d,f in os.walk(system_dir):
	#go over each directory
	for dire in d:
		#ensure that the root is the system directory
		if r == system_dir:
			print("On system: ", dire)

			#we have our system, now go through and get the original ligand placement data
			#store native ligand placement data as a dictionary
			native_ligand_coords = {}

			#read the native ligand file
			native_lig_file = open(r + "/" + dire + "/ligand.mol2","r")
			#note when we are in the atom block
			in_atom_block = False
			for line in native_lig_file.readlines():
				#determine if we are starting the atom block
				if line.startswith("@<TRIPOS>ATOM"):
					in_atom_block = True
					#continue since this line has no data
					continue

				#if we hit the bond block, we are done collecting atoms and can end
				if line.startswith("@<TRIPOS>BOND"):
					in_atom_block = False

				#if we are in the atom block, get the atom data for the line
				if in_atom_block:
					atom_name = line.split()[1]
					atom_x = float(line.split()[2])
					atom_y = float(line.split()[3])
					atom_z = float(line.split()[4])

			#test print of the native atom data
			print("Native atom data: ", native_ligand_coords)

