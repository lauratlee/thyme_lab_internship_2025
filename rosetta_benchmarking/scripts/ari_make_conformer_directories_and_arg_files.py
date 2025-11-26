#this script serves to set up directories for Rosetta discovery and generate up to 250 conformers with conformator per ligand against all residue anchors indices that Laura identified
#run this from the thyme_lab_internship_2025/rosetta_benchmarking/scripts/ directory


#imports
import os,sys

#make a dictionary that has all systems of interest and residue indices to use:
"""
systems_dict = {
    "9HZ0": [92, 94, 100],
    "9I1H": [13, 60, 92, 93, 160, 161, 218, 220, 329, 365, 367, 368, 370, 392, 394, 422, 474, 475, 481],
    "9LGN": [112, 115, 176],
    "9LYF": [9, 11, 14],
    "9M3N": [44, 45, 46, 70, 85, 86, 87, 97, 98, 99, 163, 165, 166],
    "9M7M": [21, 24, 28, 39, 44, 47, 48, 49, 50, 52, 55, 113],
    "9N48": [99, 159],
    "9QEK": [42, 91, 93, 96, 147, 148, 150],
    "9I0T": [66, 184],
    "9LNH": [171],
    "9ODR": [61, 63],
    "9LO7": [186, 234, 235, 245, 246, 282],
    "9OG3": [205, 210, 323],
    "9QEL": [18, 22, 23],
    "9M3O": [267],
    "9I4H": [200],
    "9LFU": [89],
    "9LSL": [245],
    "9QAC": [89],
    "9MZX": [144],
    "9I0W": [3, 105, 107]
}
"""
"""
#smaller dict to rerun
systems_dict = {
    "9LSL": [245],
    "9QAC": [89],
    "9MZX": [144],
    "9I0W": [3, 105, 107]
}
"""

#rerun for more accurate residues list
systems_dict = {
    "9HZ0": [26, 42, 44, 91, 100, 145],#
    "9I1H": [60, 92, 93, 218, 220, 329, 365, 392, 475, 481],#
    "9LGN": [112, 115, 176],#good no redo
    "9LYF": [14, 86, 209, 246, 250, 255, 294],#
    "9M3N": [44, 46, 86, 87, 99, 105, 111, 166],#
    "9M7M": [21, 24, 28, 39, 44, 47, 48, 49, 50, 52, 55, 113],#probably fine, no redo
    "9QEK": [42, 91, 93, 96, 147, 148, 150],#probably fine, no redo
    "9I0T": [66, 91, 184, 186],#
    "9ODR": [35, 63, 69, 83, 85],#
    "9LO7": [73, 186, 249, 253, 282, 285],#
    "9OG3": [205, 206, 208, 210, 319, 327, 387],#
    "9QEL": [17, 22, 32, 91, 128],#
    "9M3O": [267],#fine, no redo
    "9I4H": [174, 187 ,200, 267, 282],#
    "9LFU": [27, 42, 65, 86, 88, 141, 152],#
    "9LSL": [223, 234, 241, 245, 248],#
    "9QAC": [34, 38, 46, 89]#

    
}


#temp smaller dictionary to test with and make sure this works
#systems_dict = {"9QEL": [18, 22, 23]}

#iterate over each system
#for each system, make a diverse (up to 250) conformer set for the system ligand
#then separate the conformers file with obabel
#then make a folder for each conformer and make params for each conformer in the folder (1 conformer per folder for quick parallel rosetta runs)
#make new arg files for each conformer as well as a test params directory

#store the workign location, which should be thyme_lab_internship_2025/rosetta_benchmarking/scripts/
starting_script_dir = os.getcwd()

print("starting script directory, which should contain thyme_lab_internship_2025/rosetta_benchmarking/scripts", starting_script_dir)

for system in systems_dict.keys():
	print(system)

	#move to the directory for the system
	os.chdir("../system_dir/" + system)

	#store the system directory
	system_dir = os.getcwd()

	#debug print
	print("now working out of ", os.getcwd())

	#initial cleanup
	os.system("rm -drf individual_conf* confs.mol2")

	#run conformator on the ligand (always named ligand.mol2)
	os.system("/pi/summer.thyme-umw/2024_intern_lab_space/conformator_1.2.1/conformator -i ligand.mol2 -o confs.mol2 --keep3d --hydrogens -n 250 -v 0")

	#split the conformers file into multiple files
	#example files made: individual_conf_1.mol2, individual_conf_2.mol2, individual_conf_3.mol2, ... individual_conf_250.mol2
	os.system("obabel -i mol2 confs.mol2 -o mol2 -O individual_conf_.mol2 -m")

	#for each generated conformer, make a folder based on the conformer and move the file into the folder
	for r,d,f in os.walk(os.getcwd()):
		for file in f:
			if file.startswith("individual_conf_"):
				file_base = file.split(".")[0]

				#make a directory and move the file
				os.system("mkdir " + file_base)

				#move the file to the directory
				os.system("mv " + file + " " + file_base)

				#move into the directory
				os.chdir(file_base)

				#run molfile_to_params.py from our built container
				#singularity exec  /pi/summer.thyme-umw/2024_intern_lab_space/ari_work/rosetta_discovery_benchmark_test_space_september_2025/rosetta_and_conformator.sif python /rosetta/source/scripts/python/public/molfile_to_params.py
				os.system("singularity exec /pi/summer.thyme-umw/2024_intern_lab_space/ari_work/rosetta_discovery_benchmark_test_space_september_2025/rosetta_and_conformator.sif python /rosetta/source/scripts/python/public/molfile_to_params.py " + file + " -n " + file_base + " --keep-names --long-names --clobber --no-pdb")

				#prepare a test_params directory to put the conformer params file in
				os.system("mkdir test_params")

				#make empty needed files in the directory
				os.system("touch test_params/exclude_pdb_component_list.txt test_params/patches.txt")

				#move the params file in the folder
				os.system("mv " + file_base + ".params test_params")

				#write the residue_types file
				restype_file = open("test_params/residue_types.txt", "w")
				restype_file.write("## the atom_type_set and mm-atom_type_set to be used for the subsequent parameter\n")
				restype_file.write("TYPE_SET_MODE full_atom\n")
				restype_file.write("ATOM_TYPE_SET fa_standard\n")
				restype_file.write("ELEMENT_SET default\n")
				restype_file.write("MM_ATOM_TYPE_SET fa_standard\n")
				restype_file.write("ORBITAL_TYPE_SET fa_standard\n")
				restype_file.write("## Params files\n")
				restype_file.write(file_base + ".params\n")
				restype_file.close()

				#now, write the flags file for this Rosetta run

				#derive the anchor residue(s) list from the dictionary
				anchor_res_list = ""
				for resid in systems_dict[system]:
					anchor_res_list = anchor_res_list + str(resid) + ","
				#trim off tailing comma
				if anchor_res_list.endswith(","):
					anchor_res_list = anchor_res_list[:-1]

				#derive the backbone pdb file, the only pdb file in ../test_params/
				backbone_file = ""
				for r2, d2, f2 in os.walk("../test_params/"):
					for file2 in f2:
						if "chain_" in file2 and file2.endswith(".pdb"):
							backbone_file = file2

				#iterate and run over each anchor residue
				for resid in systems_dict[system]:

					os.system("mkdir " + str(resid))
					os.chdir(str(resid))

					argfile = open("args", "w")

					argfile.write("#keep seed constant\n")
					argfile.write("-constant_seed\n")

					argfile.write("#input empty receptor protein\n")
					argfile.write("-s ../../test_params/" + backbone_file + "\n")

					argfile.write("#directory of ligand(s) to attempt to dock\n")
					argfile.write("#POINT TO test_params DIRECTORY\n")
					argfile.write("-params_directory_path ../test_params/\n")

					argfile.write("#ligand motifs library\n")
					argfile.write("#this is the motifs file you will use\n")
					argfile.write("-motif_filename /pi/summer.thyme-umw/2024_intern_lab_space/FINAL_motifs_list_filtered_2_3_2023.motifs\n")

					argfile.write("#index of residue(s) to dock ligands against\n")
					argfile.write("-protein_discovery_locus " + str(resid) + " \n")

					argfile.write("#minimum cutoffs for fa_atr, fa_rep, and combined fa_atr_rep to be under\n")
					argfile.write("#keep constant\n")
					argfile.write("-fa_atr_cutoff = -2\n")
					argfile.write("-fa_rep_cutoff = 150\n")
					argfile.write("#-ddg_cutoff = -9\n")

					argfile.write("#constrain coordinates\n")
					argfile.write("#keep\n")
					argfile.write("-constrain_relax_to_start_coords\n")

					argfile.write("#keep all placements\n")
					argfile.write("-best_pdbs_to_keep 0\n")
					argfile.close()

					#bsub job throttle to make sure we do not exceed our local limit
					#write the length of the bjobs queue to this current location
					os.system("bjobs | wc -l > bjobs_length.txt")
					job_count = 0
					with open("bjobs_length.txt") as f:
						job_count = int(f.read().strip())
					while job_count > 3000:
						#sleep for 1 second to not overburden the system
						os.system("sleep 1")
						os.system("bjobs | wc -l > bjobs_length.txt")
						with open("bjobs_length.txt") as f:
							job_count = int(f.read().strip())
					#remove the length file to avoid clutter
					os.system("rm bjobs_length.txt")

					#now, run Rosetta in a bsub job
					os.system("bsub -q long -M 5000 -R \"rusage[mem=5000]\" -W 8:00 \"/pi/summer.thyme-umw/2024_intern_lab_space/rosetta/source/bin/ligand_discovery_search_protocol.linuxgccrelease @args\"")

					os.chdir("..")


				#at end, move back up
				os.chdir(system_dir)

	#return to the top
	os.chdir(starting_script_dir)