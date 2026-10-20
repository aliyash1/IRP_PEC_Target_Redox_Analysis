from Bio.PDB import MMCIFParser
from Bio.PDB.NeighborSearch import NeighborSearch
from pathlib import Path
import json
import requests
from collections import defaultdict 
import os

main_dir = Path("D:/IRP Lucas")
alphafill_dir = main_dir / "alphafill_files"
cofactors_for_all_targets = []
uniprot_link = "https://rest.uniprot.org/uniprotkb/"

def find_residues_for_alphafill_pecs(alphafill_dir):
    ### Gets the cofactor residues from the Alphafill files. 
    ### It takes all residues within 4A of the cofactor

    #loop over all cif files in the alphafill directory
    for protein in alphafill_dir.glob("*.cif"):
        ### Collect information about the protein and cofactor

        file_name = protein.name
        parts = file_name.split("_")
        pec_name = f"{parts[0]}_{parts[1]}"
        cofactor = parts[2]
        cofactor_split = cofactor.split(".")
        cofactor = cofactor_split[0]
        label = 1       # only 1 cofacor put in so label = 1 for all

        ### Find residues

        # load structure
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("prot", protein)

        model = structure[0]

        protein_atoms = []
        ligand_atoms = []

        # separate protein vs ligand
        for chain in model:
            for residue in chain:
                hetflag = residue.id[0]

                for atom in residue:
                    if hetflag == " ":   # standard amino acid
                        protein_atoms.append(atom)
                    else:                 # cofactor / ligand
                        ligand_atoms.append(atom)

        # build spatial index
        ns = NeighborSearch(protein_atoms)

        # find interacting residues
        contact_residues = set()

        cutoff = 4.0  # Angstroms

        for atom in ligand_atoms:
            neighbors = ns.search(atom.coord, cutoff)

            for n in neighbors:
                contact_residues.add(n.get_parent())
        
        residue_numbers = [res.id[1] for res in contact_residues]

        ### Add residues to list

        cofactors_for_all_targets.append({
            "protein": pec_name,
            "cofactors": [{
                "name": cofactor,
                "label": label,
                "residues": residue_numbers
            }]
        })

def find_residues_for_regular_proteins(cofactors_for_all_targets, uniprot_link):

    ### Finds the cofactor residues from UniProt and adds them to the list

    ### PECs

    existing_proteins = {item["protein"] for item in cofactors_for_all_targets}

    #Loop over all pecs in dataset
    with open("JSON_files/final_pec_dataset.json") as f:
        pecs = json.load(f)

    for entry in pecs:

        # Get info of pec
        pec_name = entry["uniProtkbId"]
        pec_id = entry["primaryAccession"]

        # If PEC has a Aplhafill cofactor, this will give an empty result
        if pec_name in existing_proteins:
            continue

        # Collect data from UniProt
        response = requests.get(
            f"{uniprot_link}{pec_id}.json?fields=ft_binding"
        )
        pec_cofactor_data = response.json()

        # Structure data and add it to the list
        cofactor_sites = defaultdict(list)

        for feature in pec_cofactor_data.get("features", []):
            if "binding site" in feature.get("type", "").lower():

                residue = feature["location"]["start"]["value"]

                ligand = feature.get("ligand", {})
                label = int(ligand.get("label") or 1)
                name = ligand.get("name", "unknown")

                cofactor_sites[(name, label)].append(residue)

        cofactors = []
        for (name, label), residues in cofactor_sites.items():
            cofactors.append({
                "name": name,
                "label": label,
                "residues": sorted(set(residues))
            })

        cofactors_for_all_targets.append({
            "protein": pec_name,
            "cofactors": cofactors
        })

    ### Targets

    # Loop over all targets in dataset
    with open("JSON_files/final_target_proteins_dataset.json") as f:
        targets = json.load(f)

    for entry in targets:

        # Get info of target
        target_name = entry["uniProtkbId"]
        target_id = entry["primaryAccession"]

        existing_proteins = {item["protein"] for item in cofactors_for_all_targets}
        if target_name in existing_proteins:
            print(f"Skipping duplicate target: {target_name}")
            continue

        # Get data from UniProt
        response = requests.get(
            f"{uniprot_link}{target_id}.json?fields=ft_binding"
        )
        target_cofactor_data = response.json()

        # Structure data and add it to list
        cofactor_sites = defaultdict(list)

        for feature in target_cofactor_data.get("features", []):
            if "binding site" in feature.get("type", "").lower():

                residue = feature["location"]["start"]["value"]

                ligand = feature.get("ligand", {})
                try:
                    label = int(ligand.get("label") or 1)
                except:
                    label = 1 
                    # Only needed for TKT1 cofactors D-ribose 5-phosphate and D-xylulose 5-phosphate
                    # And TKT1 only has 1 of these so can be fixed by setting label 1 
                name = ligand.get("name", "unknown")

                cofactor_sites[(name, label)].append(residue)

        cofactors = []
        for (name, label), residues in cofactor_sites.items():
            cofactors.append({
                "name": name,
                "label": label,
                "residues": sorted(set(residues))
            })

        cofactors_for_all_targets.append({
            "protein": target_name,
            "cofactors": cofactors
        })


find_residues_for_alphafill_pecs(alphafill_dir)
find_residues_for_regular_proteins(cofactors_for_all_targets, uniprot_link)

#Save
output_dir = Path("final_datasets/working_data")
os.makedirs(output_dir, exist_ok=True)
output_file = output_dir / "cofactor_residues_for_all_proteins.json"
with open(output_file, "w") as json_file:
        json.dump(cofactors_for_all_targets, json_file, indent=4)