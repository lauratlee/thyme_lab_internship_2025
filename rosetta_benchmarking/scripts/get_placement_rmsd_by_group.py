# get_placement_rmsd_best_by_group.py
# Usage: python get_placement_rmsd_best_by_group.py <system> <residue_folder> <group_number>

import os, sys
import pymol2
from rdkit import Chem
from rdkit.Chem import rdMolAlign
from openbabel import pybel
from pymol import CmdException

# Helper function to strip all bond orders
def strip_bond_orders(mol):
    rw_mol = Chem.RWMol(mol)
    for bond in rw_mol.GetBonds():
        bond.SetBondType(Chem.rdchem.BondType.SINGLE)
    Chem.SanitizeMol(rw_mol)
    return rw_mol.GetMol()

# User arguments
target_system = sys.argv[1]
residue_folder = sys.argv[2]
target_group = sys.argv[3]

# Start PyMOL session
with pymol2.PyMOL() as pymol:
    cmd = pymol.cmd

    # Enter system directory
    if not os.path.isdir(target_system):
        print(f"System directory {target_system} not found.")
        sys.exit(1)
    os.chdir(target_system)

    print(f"Processing system: {target_system}")

    # Prepare reference ligand
    os.system("obabel ligand.mol2 -O ligand.sdf")
    ref_ligand = Chem.MolFromMolFile("ligand.sdf", removeHs=True)
    ref_ligand = strip_bond_orders(ref_ligand)
    if ref_ligand is None:
        print("WARNING: reference ligand did not read into RDKit. Exiting.")
        sys.exit(1)

    # Load reference structure in PyMOL
    orig_file = f"{target_system}.pdb"
    cmd.load(orig_file, "reference")
    if cmd.count_atoms("reference") == 0:
        print("WARNING: no atoms in reference. Exiting.")
        sys.exit(1)
    print(f"ATOMS IN REFERENCE: {cmd.count_atoms('reference')}")

    # Locate residue and group folder
    group_path = os.path.join(residue_folder, target_group)
    if not os.path.isdir(group_path):
        print(f"Group folder {group_path} not found.")
        sys.exit(1)

    print(f"Processing residue {residue_folder}, group {target_group}")

    # Track the best RMSD for this group
    best_rmsd = None
    best_file = None

    # Iterate through placement files
    for group_file in os.listdir(group_path):
        if ".pdb" not in group_file:
            continue

        # Load placement into PyMOL
        cmd.load(os.path.join(group_path, group_file), "placement")
        if cmd.count_atoms("placement") == 0:
            cmd.delete("placement")
            continue

        # Align to reference
        try:
            cmd.align("placement", "reference")
        except CmdException:
            cmd.delete("placement")
            continue

        # Save aligned ligand
        cmd.select("aligned_lig", "placement and not polymer.protein")
        aligned_lig_basename = group_file.split(".")[0] + "_aligned_lig.mol2"
        cmd.save(os.path.join(group_path, aligned_lig_basename), "aligned_lig")
        cmd.delete("aligned_lig")
        cmd.delete("placement")

        # Convert to SDF for RDKit
        aligned_lig_sdf_basename = group_file.split(".")[0] + "_aligned_lig.sdf"
        os.system(f"obabel {group_path}/{aligned_lig_basename} -O {group_path}/{aligned_lig_sdf_basename}")

        # Read ligand in RDKit
        pla_ligand = Chem.MolFromMolFile(f"{group_path}/{aligned_lig_sdf_basename}", removeHs=True)
        pla_ligand = strip_bond_orders(pla_ligand)
        if pla_ligand is None:
            continue

        # Compute RMSD
        try:
            rmsd = rdMolAlign.CalcRMS(ref_ligand, pla_ligand)
        except Exception:
            continue

        # Track the best RMSD
        if best_rmsd is None or rmsd < best_rmsd:
            best_rmsd = rmsd
            best_file = group_file

    # Write only the best RMSD to a CSV
    output_csv = f"{target_system}_{residue_folder}_group{target_group}_best_placement.csv"
    if best_rmsd is not None:
        with open(output_csv, "w") as f:
            f.write("residue,file,rmsd\n")
            f.write(f"{residue_folder},{best_file},{best_rmsd}\n")
        print(f"Best RMSD written to {output_csv}")
    else:
        print(f"No valid placements found for group {target_group}")

    # Cleanup
    cmd.delete("reference")
    print(f"{target_system} {residue_folder} group {target_group} DONE")
