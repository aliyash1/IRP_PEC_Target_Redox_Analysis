import json
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser

# Directories
base_dir = Path(__file__).resolve().parent
project_dir = base_dir.parent.parent
monomer_dir = project_dir / "Inputs" / "Structures" / "PEC_monomers"   # PEC_monomers is not actually monomers
# PEC-target cif is a 2-chain complex, but the script only reads structure[0]["A"]
# 'A' is the FER_ECOLI chain , ACKA_ECOLI chain B is not read

output_dir = project_dir / "Results" / "2_manual" / "Binding_CM"

parser = MMCIFParser()


with open("Inputs/final_pec_dataset_fes_only.json") as f:
    pecs = json.load(f)
 
with open("Inputs/cofactor_residues_FeS_only.json") as f:
    cofactor_data = json.load(f)
cofactor_data = {item["protein"]: item for item in cofactor_data}
 

# ==========================================
# get_all_distances.py - Lucas' function
# ==========================================
def get_cofactor_coordinates(structure, chain_id, cofactor_sites):

    chain = structure[0][chain_id]

    cofactor_coordinates = []

    for cofactor in cofactor_sites["cofactors"]:

        name = cofactor["name"]
        label = cofactor["label"]
        residue_nums = cofactor["residues"]

        coordinates = []

        for residue_num in residue_nums:
            residue_id = (" ", residue_num, " ")
            residue = chain[residue_id]
            atom = residue["CA"]
            coord = atom.get_coord()
            coordinates.append(coord)

        cofactor_coordinate = np.mean(coordinates, axis=0)

        cofactor_coordinates.append(
            (name, label, cofactor_coordinate)
        )
    return cofactor_coordinates


# ================================================
# Find PEC monomer structure file in monomer_dir
# ================================================
def find_monomer_structure(pec_name, monomer_dir):
    matches = list(monomer_dir.glob(f"*{pec_name}*.cif"))
    if not matches:
        return None
    return matches[0]

# I manually renamed sp_p0a9r4_fer_ecoli_sp_p0a6a3_acka_ecoli_model.cif to FER_ECOLI_model.cif
 

 
# ===========================================
# Main loop - binding-site centroid per PEC
# ===========================================

centroid_rows = []
 
for entry in pecs:
 
    pec_name = entry["uniProtkbId"]
 
    protein_entry = cofactor_data.get(pec_name)
    if protein_entry is None:
        print(f"skipping {pec_name}: not in cofactor_data")
        continue
 
    structure_path = find_monomer_structure(pec_name, monomer_dir)
    if structure_path is None:
        print(f"skipping {pec_name}: no structure file found in {monomer_dir}")
        continue
 
    try:
        structure = parser.get_structure("protein", structure_path)
    except Exception as e:
        print(f"skipping {pec_name}: could not parse {structure_path} ({e})")
        continue
 
    try:
        # PEC monomer structures -> single chain, same "A" convention as Lucas used
        coords = get_cofactor_coordinates(structure, "A", protein_entry)
    except KeyError as e:
        print(f"skipping {pec_name}: binding-site residue not found in structure ({e})")
        continue
 
    for name, label, com in coords:
        centroid_rows.append({
            "pec_name": pec_name,
            "cofactor_name": name,
            "cofactor_label": label,
            "centroid_x": float(com[0]),
            "centroid_y": float(com[1]),
            "centroid_z": float(com[2]),
        })
 
### Write the output file
output_dir.mkdir(parents=True, exist_ok=True)
output_file = output_dir / "pec_binding_site_centroids.csv"
centroid_df = pd.DataFrame(centroid_rows)
centroid_df.to_csv(output_file, index=False)

print(f"There are {len(centroid_df)} PEC binding-site centroids to {output_file}")