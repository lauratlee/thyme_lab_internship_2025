
# calculates the best RMSD for the outputs of each system and stores it in a .csv file (format SYSTEM, BEST_RMSD)
# run from system_dir/

import os,sys

csv_path = "best_rmsds.csv"

rmsds_dict = {}


def get_rmsd(result_list):
    num_atoms = 0
    native_ligand = {}

    lig_path = "../ligand.mol2"
    
    # reads native ligand file for reference coordinates
    with open(lig_path, 'r') as lig_file:
        
        atom_start_index = 7
        
        for i, line in enumerate(lig_file):
            # reads number of atoms in ligand from file
            if i == 2:
                atom_count_line = line.strip()
                atom_count_line_clean = line.split("    ")
                num_atoms = int(atom_count_line_clean[0])
                

            # adds coordinates of atoms in native ligand to native_ligand dictionary
            if atom_start_index <= i < (atom_start_index + num_atoms):
                atom_id = line.strip()[:8].strip()
                atom_coordinates = line.strip()[8:].strip().split("  ")

                for value in atom_coordinates:
                    if value == "":
                        atom_coordinates.pop(atom_coordinates.index(value))
                
                if "H" not in atom_id:
                    atom_x = float(atom_coordinates[0].strip())
                    atom_y = float(atom_coordinates[1].strip())
                    atom_z = float(atom_coordinates[2].strip(" ").split()[0])                           
                    native_ligand[atom_id] = [atom_x, atom_y, atom_z]
        lig_file.close()

        # add placement files to a list for analysis
        file_list = []
        for _, _, files in os.walk(os.getcwd()):
            for file in files:
                if ".pdb" in file and os.path.getsize(file) > 0:
                    file_list.append(file)

        # iterate through placement files
        for out_file in file_list:
            res_rmsd_list = []
            gen_rmsd_list = []
            placement = {}
            best_res_rmsd = None

            with open(out_file, 'r') as output:
                for i, line in enumerate(output):
                    # filter for ligand placement coordinates only
                    if "HETATM" in line and "lig" in line: 
                        out_atom_id = int(line[6:11].strip())
    
                        split_line = [x for x in line.split(" ") if x.strip() != ""]
    
                        out_atom_x = float(split_line[6])
                        out_atom_y = float(split_line[7])
                        out_atom_z = float(split_line[8].split()[0])
    
                        placement[out_atom_id] = [out_atom_x, out_atom_y, out_atom_z]
    
                        non_hydrogens = len(placement)
    
                        distance_sum = 0
    
                    # iterates through all atoms and calculate distance between output and native placement
                    for atom in placement.keys():
                        a1 = placement[atom]
                        a2 = native_ligand[atom]
    
                        x = (a1[0] - a2[0])**2
                        y = (a1[1] - a2[1])**2
                        z = (a1[2] - a2[2])**2
    
                        dist = x+y+z
    
                        sqrt_dist = dist**0.5
                        
                        distance_sum += sqrt_dist
    

                    # calculate rmsd and select best rmsd for that specific residue, then add that value to a list of best (residue-specific) rmsds
                    rmsd = distance_sum/non_hydrogens
                    res_rmsd_list.append(rmsd)
            output.close()

        best_res_rmsd = float(min(res_rmsd_list))
        result_list.append(best_res_rmsd)

        return(result_list)


for system in os.listdir(os.getcwd()):
    # exclude 9LE4, which had 2 identical ligands and was not optimal for discovery
    if system == "9LE4":
        continue 
    print(system)
    os.chdir(system)

    gen_rmsd_list = []

    # iterate thru all residues for system
    for dir in os.listdir(os.getcwd()):
        if "res_" in dir:
            os.chdir(dir)
            # update gen_rmsd_list with best rmsd value for that residue
            get_rmsd(gen_rmsd_list)
            os.chdir("..")

    # calculate best rmsd across all anchor residues
    best_rmsd = float(min(gen_rmsd_list))

    rmsds_dict[system] = best_rmsd
    os.chdir("..")


with open(csv_path, 'w') as rmsd_file:
    for key in rmsds_dict.keys():
        rmsd_file.write(str(key) + "," + str(rmsds_dict[key]) + "\n")
        






                    

                    
                            
                    

                    

                

            
                 
        







                                                                                       
