# Inputs are AlphaFold protein monomers from 3 folders: AF_PECs, AF_Targets, AF_Prots_both_types
# And json file with cofactor binding-site residue numbers for each protein (cofactor_residues_FeS_only.json)

# The output is a table of binding-site centroids for all proteins
# The centroid calculation is per cofactor per protein (if a protein has multiple cofactors, there will be multiple rows for that protein in the output table).
# Output columns:
# protein, source_folder, cofactor_name, cofactor_label, centroid_x/y/z, and structure_file.


import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser

# Directories
base_dir = Path(__file__).resolve().parent
project_dir = base_dir.parent.parent
af_pec_dir = project_dir / "Inputs" / "Structures" / "Protein_AF_structures" / "AF_PECs"
af_target_dir = project_dir / "Inputs" / "Structures" / "Protein_AF_structures" / "AF_Targets"
af_both_dir = project_dir / "Inputs" / "Structures" / "Protein_AF_structures" / "AF_Prots_both_types"

# label used in the output table for where each structure came from
af_dirs = {
    "AF_PECs": af_pec_dir,
    "AF_Targets": af_target_dir,
    "AF_Prots_both_types": af_both_dir,
}

# The AF structures are AlphaFold Server monomer outputs (single chain "A"), e.g.
#   fold_ygfs_ecoli/fold_ygfs_ecoli_model_0.cif ... model_4.cif  (+ full_data, summary_confidences, msas/, templates/)
# model_0 is the top-ranked model, so that is the one used.

output_dir = project_dir / "Results" / "2_manual" / "Binding_CM"

parser = MMCIFParser()


with open(project_dir / "Inputs" / "json_files"/"cofactor_residues_FeS_only.json") as f:
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


# ================================================================
# List every protein that has a top-ranked AF model in a folder
# ================================================================
# AlphaFold Server file names look like fold_<protein name, lowercase>_model_<0-4>.cif

MODEL_0_PATTERN = re.compile(r"^fold_(.+)_model_0\.cif$", re.IGNORECASE)

def find_structures(af_dir):
    """Return {PROTEIN_NAME: path to its model_0 cif} for everything below af_dir."""
    structures = {}
    if not af_dir.exists():
        print(f"warning: folder not found: {af_dir}")
        return structures

    for path in sorted(af_dir.rglob("*.cif")):
        match = MODEL_0_PATTERN.match(path.name)
        if not match:
            continue
        protein_name = match.group(1).upper()
        if protein_name in structures:
            print(f"warning: {protein_name} has several model_0 files in {af_dir.name}, using {structures[protein_name]}")
            continue
        structures[protein_name] = path
    return structures


# =========================================================================
# Main loop: binding-site centroid per protein, over all three folders
# =========================================================================

centroid_rows = []
seen = {}   # protein -> folders it was found in, to flag duplicates across folders

for source_name, af_dir in af_dirs.items():

    structures = find_structures(af_dir)
    print(f"{source_name}: found {len(structures)} proteins")

    for protein_name, structure_path in structures.items():

        seen.setdefault(protein_name, []).append(source_name)

        protein_entry = cofactor_data.get(protein_name)
        if protein_entry is None:
            print(f"skipping {protein_name} ({source_name}): not in cofactor_data")
            continue

        try:
            structure = parser.get_structure("protein", structure_path)
        except Exception as e:
            print(f"skipping {protein_name} ({source_name}): could not parse {structure_path} ({e})")
            continue

        try:
            # monomer structures -> single chain "A", same convention as Lucas used
            coords = get_cofactor_coordinates(structure, "A", protein_entry)
        except KeyError as e:
            print(f"skipping {protein_name} ({source_name}): binding-site residue not found in structure ({e})")
            continue

        for name, label, com in coords:
            centroid_rows.append({
                "protein": protein_name,
                "source_folder": source_name,
                "cofactor_name": name,
                "cofactor_label": label,
                "centroid_x": float(com[0]),
                "centroid_y": float(com[1]),
                "centroid_z": float(com[2]),
                "structure_file": str(structure_path.relative_to(project_dir)),
            })

duplicates = {p: f for p, f in seen.items() if len(f) > 1}
if duplicates:
    print(f"note: {len(duplicates)} proteins appear in more than one folder (one row per folder in the table):")
    for p, f in duplicates.items():
        print(f"   {p}: {', '.join(f)}")

### Write the output file
output_dir.mkdir(parents=True, exist_ok=True)
output_file = output_dir / "binding_site_centroids_all_proteins.csv"
centroid_df = pd.DataFrame(centroid_rows)
centroid_df.to_csv(output_file, index=False)

print(f"Wrote {len(centroid_df)} binding-site centroids for {centroid_df['protein'].nunique() if len(centroid_df) else 0} proteins to {output_file}")