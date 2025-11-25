#the purpose of this script is to go system by system and determine the best RMSDS for the best, top 10, and all placements per system using the 250 conformer schema
#this script is to be run in the scripts directory like: python ari_get_rmsds_rosetta.py

#if you only want to look at a single system or group of systems, add them after the call

#import
import os,sys

#collect select systems if there are any
select_systems = []

if len(sys.argv) > 1:
	#add each specific system
	for i in range(len(sys.argv)):
		if i > 0:
			select_systems.append(sys.argv[i])

#print the select systems for confirmation
print("Selected systems:")
print(select_systems)

#create a dictionary that holds the systems and a list for the best rmsds per system
systems_rmsds = {}

#store the workign location, which should be thyme_lab_internship_2025/rosetta_benchmarking/scripts/
starting_script_dir = os.getcwd()

print("starting script directory, which should contain thyme_lab_internship_2025/rosetta_benchmarking/scripts", starting_script_dir)

#move to the system directory
os.chdir("../system_dir/")

#store the system directory
system_dir = os.getcwd()

#make a file to write the results
write_file = open("all_ari_rosetta_placement_results.txt", "w")

#iterate over the systems directories and go system by system, collecting the best rmsds
for r,d,f in os.walk(system_dir):
	#go over each directory
	for dire in d:
		#ensure that the root is the system directory
		if r == system_dir:
			#ensure the directory is in the select systems and we are selecting systems
			if dire not in select_systems and len(select_systems) > 0:
				print(dire + " not in selected systems.")
				continue

			print("On system: ", dire)

			#create a list that holds all placements with the rmsd, ddg, and system file
			#this will get sorted at the end by ddg to derive the top ddg for the best of all, 10, and single lowest ddg
			placements_data = []

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

					#skip hydrogens
					if atom_name.startswith("H"):
						continue

					native_ligand_coords[atom_name] = [atom_x,atom_y,atom_z]

			#derive the center of mass of the ligand
			lig_com = [0,0,0]
			natoms = 0
			for atom in native_ligand_coords.keys():
				lig_com[0] = lig_com[0] + native_ligand_coords[atom][0]
				lig_com[1] = lig_com[1] + native_ligand_coords[atom][1]
				lig_com[2] = lig_com[2] + native_ligand_coords[atom][2]
				natoms = natoms + 1

			if natoms != 0:
				lig_com[0] = lig_com[0] / natoms
				lig_com[1] = lig_com[1] / natoms
				lig_com[2] = lig_com[2] / natoms

			#test print of the native atom data
			print("Native atom data: ", native_ligand_coords)

			#now, read through each system and derive the placement
			for r2,d2,f2 in os.walk(r + "/" + dire):
				for placement_file in f2:
					if "_individual_conf_" in placement_file and placement_file.endswith(".pdb"):
						#we have a placement file; derive the ligand data
						placement_ligand_coords = {}

						#variable to store teh ddg
						ddg = 100
						anchor_res = ""

						placement_read_file = open(r2 + "/" + placement_file, "r")
						for line in placement_read_file.readlines():
							#the ligand lines start with HETATM
							if line.startswith("HETATM"):
								#derive the atom data
								atom_name = line.split()[2]
								atom_x = float(line[30:38])
								atom_y = float(line[38:46])
								atom_z = float(line[46:54])		
														
								if atom_name.startswith("H"):
									continue

								placement_ligand_coords[atom_name] = [atom_x,atom_y,atom_z]

							if line.startswith("Scoring: Post-HighResDock system ddG:"):
								ddg = float(line.strip().split()[4])

							if line.startswith("Table Values "):
								anchor_res = line.split(",")[5]


						#now, derive the placement rmsd
						distance_sum = 0
						natoms = 0

						for atom in placement_ligand_coords.keys():
							if atom in native_ligand_coords.keys():
								#derive the distance 
								distance = (((placement_ligand_coords[atom][0] - native_ligand_coords[atom][0]) ** 2) + ((placement_ligand_coords[atom][1] - native_ligand_coords[atom][1]) ** 2) + ((placement_ligand_coords[atom][2] - native_ligand_coords[atom][2]) ** 2)) ** 0.5

								#add the distance to the distance sum and increment the natoms
								distance_sum = distance_sum + distance
								natoms = natoms + 1

						#derive the rmsd
						rmsd = 100

						pla_com = [0,0,0]
						natoms = 0
						for atom in placement_ligand_coords.keys():
							pla_com[0] = pla_com[0] + placement_ligand_coords[atom][0]
							pla_com[1] = pla_com[1] + placement_ligand_coords[atom][1]
							pla_com[2] = pla_com[2] + placement_ligand_coords[atom][2]
							natoms = natoms + 1

						if natoms != 0:
							pla_com[0] = pla_com[0] / natoms
							pla_com[1] = pla_com[1] / natoms
							pla_com[2] = pla_com[2] / natoms						


						if natoms > 0:
							#rmsd = distance_sum / natoms
							rmsd = (((lig_com[0] - pla_com[0]) ** 2) + ((lig_com[1] - pla_com[1]) ** 2) + ((lig_com[2] - pla_com[2]) ** 2)) ** 0.5

						#send values to the list
						placements_data.append([rmsd,ddg,r2 + "/" + placement_file,anchor_res])

			if len(placements_data) == 0:
				print("NO PLACEMENTS FOR THIS SYSTEM!")
				continue

			#done looking at placements, derive the best rmsd for the ddg cutoffs
			#sort the placements_data by ddg
			placements_data = sorted(placements_data, key=lambda x: x[1])

			#test print of top 10 ddg
			#for i in range(0,10):
			#	print(placements_data[i])

			#determine the best rmsd for all, top 10, and top 1
			top_1_rmsd = placements_data[0]

			#top 10
			#default to keeping the first
			top_10_rmsd = placements_data[0]
			for i in range(0,10):
				#safety check to make sure we don't run out of bounds in case there are fewer than 10 palcements (which seems unlikely)
				if len(placements_data) > i:
					#check if better than currently held
					if placements_data[i][0] < top_10_rmsd[0]:
						top_10_rmsd = placements_data[i]


			#top all
			top_all_rmsd = placements_data[0]
			for i in range(0,len(placements_data)):

				#check if better than currently held
				if placements_data[i][0] < top_all_rmsd[0]:
					top_all_rmsd = placements_data[i]

			#return the top rmsds
			print("top all: ", top_all_rmsd)
			print("top 10: ", top_10_rmsd)
			print("top 1: ", top_1_rmsd)

			write_file.write(dire + "\n")
			write_file.write(f"top all: {top_all_rmsd}\n")
			write_file.write(f"top 10: {top_10_rmsd}\n")
			write_file.write(f"top 1: {top_1_rmsd}\n")

			#open and write a system_specific file to write the data for too
			system_write_file = open(r + "/" + dire + "/" + dire + "_best_rmsds.txt","w")
			system_write_file.write(dire + "\n")
			system_write_file.write(f"top all: {top_all_rmsd}\n")
			system_write_file.write(f"top 10: {top_10_rmsd}\n")
			system_write_file.write(f"top 1: {top_1_rmsd}\n")

			#now, run through each anchor residue index and determine the best per anchor residue, to determine if any anchors are bad and are just placing in a different pocket
			#iterate over the sorted placements data and rehash the data into a dictionary for easier processing

			placements_by_anchor_dict = {}

			for placement in placements_data:
				#skip if somehow no anchor
				if placement[3] == "":
					continue

				#see if the anchor is a key, and if not, create a new key list
				if placement[3] not in placements_by_anchor_dict.keys():
					placements_by_anchor_dict[placement[3]] = []

				#append the placement (should still be in ddg order)
				placements_by_anchor_dict[placement[3]].append(placement)

			#iterate over the anchor dict and then determine the best rmsd for each of the 3 groups and print/write
			for anchor in placements_by_anchor_dict.keys():
				#determine the best rmsd for all, top 10, and top 1
				top_1_rmsd = placements_by_anchor_dict[anchor][0]

				#top 10
				#default to keeping the first
				top_10_rmsd = placements_by_anchor_dict[anchor][0]
				for i in range(0,10):
					#safety check to make sure we don't run out of bounds in case there are fewer than 10 palcements (which seems unlikely)
					if len(placements_by_anchor_dict[anchor]) > i:
						#check if better than currently held
						if placements_by_anchor_dict[anchor][i][0] < top_10_rmsd[0]:
							top_10_rmsd = placements_by_anchor_dict[anchor][i]


				#top all
				top_all_rmsd = placements_by_anchor_dict[anchor][0]
				for i in range(0,len(placements_by_anchor_dict[anchor])):

					#check if better than currently held
					if placements_by_anchor_dict[anchor][i][0] < top_all_rmsd[0]:
						top_all_rmsd = placements_by_anchor_dict[anchor][i]

				#return the top rmsds
				print("anchor res " + str(anchor))
				print("top all: ", top_all_rmsd)
				print("top 10: ", top_10_rmsd)
				print("top 1: ", top_1_rmsd)

				write_file.write(dire + "\n")
				write_file.write("anchor res " + str(anchor) +"\n")
				write_file.write(f"top all: {top_all_rmsd}\n")
				write_file.write(f"top 10: {top_10_rmsd}\n")
				write_file.write(f"top 1: {top_1_rmsd}\n")

				#open and write a system_specific file to write the data for too
				#system_write_file = open(r + "/" + dire + "/" + dire + "_best_rmsds.txt","w")
				system_write_file.write("anchor res " + str(anchor) +"\n")
				system_write_file.write(dire + "\n")
				system_write_file.write(f"top all: {top_all_rmsd}\n")
				system_write_file.write(f"top 10: {top_10_rmsd}\n")
				system_write_file.write(f"top 1: {top_1_rmsd}\n")


