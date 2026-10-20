import json
import pandas as pd
import os
from Bio.PDB import MMCIFParser
import numpy as np
from pathlib import Path

### Settings
parser = MMCIFParser()
num_jobs = 6
output_dir = Path("final_datasets/output_datasets")


with open("JSON_files/final_pec_dataset.json") as f:
    pecs = json.load(f)

with open("final_datasets/working_data/cofactor_residues_for_all_proteins.json") as f:
    cofactor_data = json.load(f)
cofactor_data = {item["protein"]: item for item in cofactor_data}

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

def get_cofactor_distances(cofactor_coordinates_pec, cofactor_coordinates_target, pec_name, target_name):
    distances = []

    for (cofactor_name1, cofactor_label1, coord1) in cofactor_coordinates_pec:
        for (cofactor_name2, cofactor_label2, coord2) in cofactor_coordinates_target:
            
            distance = np.linalg.norm(coord1 - coord2)

            distances.append({
                "pec_name" : pec_name,
                "pec_cofactor_name": cofactor_name1,
                "pec_cofactor_label": cofactor_label1,
                "target_name": target_name,
                "target_cofactor_name": cofactor_name2,
                "target_cofactor_label": cofactor_label2,
                "distance": float(distance)
            })
    return distances

def get_all_distances(pec_name, parser, num_jobs, cofactor_data):
    pec_cofactor_residues = cofactor_data[f"{pec_name}"]
    distances = []

    ### Loop over jobs
    for ii in range(num_jobs):

        test_counter = 1

        # Open AF3 input file
        with open(fr"all_jobs\{pec_name}\data\job_{ii+1}.json") as f:
            input_data = json.load(f)
        input_data = pd.DataFrame(input_data)

        # Visit every fold in job
        for fold in input_data.itertuples(index=False):
            
            print(f"busy with pec: {pec_name}, job: {ii+1}, fold: {test_counter}")
            
            # Get structure as parsed CIF
            fold_name = fold.name.lower()
            job_dir = f"output_job_{ii+1}_{pec_name}"
            fold_dir = fold_name
            fold_path = os.path.join(AF3_outputs_dir, job_dir, fold_dir, f"{fold_name}_model.cif")
            try: 
                structure = parser.get_structure("protein", fold_path)
            except:
                continue
            # get target ID and name
            parts = fold_name.split("_")
            target_name = f"{parts[6]}_{parts[7]}".upper()

            # Get coordinates that bind to the cofacor(s) of target and pec
                # In the pec the residues are always the same, so outside loop
            target_cofactor_residues = cofactor_data[f"{target_name}"]

        
            # Position differs every fold, so pec is inside the loop
            pec_cofactor_coordinates = get_cofactor_coordinates(structure, "A", pec_cofactor_residues)
            target_cofactor_coordinates = get_cofactor_coordinates(structure, "B", target_cofactor_residues)
            
            # Calculate distances
            cofactor_distances = get_cofactor_distances(pec_cofactor_coordinates, target_cofactor_coordinates, pec_name, target_name)
            
            distances.append(cofactor_distances)

            # Test only!!
            test_counter += 1

    return distances

for entry in pecs:

    ### Get info of the PEC
    pec_name = entry["uniProtkbId"]
    pec_id = entry["primaryAccession"]

    AF3_outputs_dir = Path(fr"D:\IRP Lucas\PECs\{pec_name}")
    output_pec_dir = output_dir / pec_name

    ### Get and save all distances of PEC
    distances = get_all_distances(pec_name, parser, num_jobs, cofactor_data)

    output_file = output_pec_dir / f"distances_for_{pec_name}.json"

    with open(output_file, "w") as json_file:
        json.dump(distances, json_file, indent=4)






