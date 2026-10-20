# Centroid and Centre of mass of Fe-S clusters
# Input: RCSB ligand CIF files

# 1. Read the element and the x, y, z position of every atom in file
    # each ligand CIF stores kinds of two coordinate sets from the same atoms:
        # model: _chem_comp_atom.model_Cartn_{x,y,z}
        # ideal: _chem_comp_atom.pdbx_model_Cartn_{x,y,z}_ideal

# 2. Calculate geometic centre and centre of mass (for model and ideal coordinats)
# 3. Print the results, saves them to ligand_com.csv.

import csv
from pathlib import Path
import numpy as np
from Bio.PDB.MMCIF2Dict import MMCIF2Dict  # reads a CIF file into a Python dict

# ----------------------------------------------------------------------------
# Initialise
# ----------------------------------------------------------------------------

# Folders
BASE_DIR = Path(__file__).resolve().parent
PROJECT_DIR = BASE_DIR.parent.parent
INPUT_DIR = PROJECT_DIR / "Inputs" / "Structures" / "Cofactor_structures"
# INPUT_DIR = (Path(__file__).resolve().parents[2] / "Inputs" / "Structures" / "Cofactor_structures")


# Dicitonary of atomic masses in amu (g/mol)
ATOMIC_MASS = {
    "H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999,
    "S": 32.06, "FE": 55.845,
}

# Match labels to RCSB names
LIGANDS = {
    "Fe(2+)":                     "FE2",
    "[2Fe-2S] cluster":           "FES",
    "[3Fe-4S] cluster":           "F3S",
    "[4Fe-4S] cluster":           "SF4",
    "hybrid [4Fe-2O-2S] cluster": "T2N",
}

# ----------------------------------------------------------------------------
# Get coordinates & parse
# ----------------------------------------------------------------------------

def ligand_id_of(cif_path: Path) -> str:
    """
    Extracts the  the ligand ID written inside the CIF file
    """
    d = MMCIF2Dict(str(cif_path))
    for tag in ("_chem_comp.id", "_chem_comp_atom.comp_id"):
        if tag in d:
            return as_list(d[tag])[0].upper()
    return cif_path.stem.upper()


def as_list(x):
    """
    The atoms of a ligand should be stored in a list (even if there is one atom).
    Single-atom files return a raw string instead of a list from MMCIF2Dict
    """
    return x if isinstance(x, list) else [x] # object x should be of list type


def read_atoms(cif_path: Path, kind: str):
    """
    For each ligand, return
        elements: list of element symbols
        coordinates: numpy array of shape (n_atoms, 3) with x, y, z in angstroms
            Each column in this array corresponds to a particular atom.
            Based on the kind (model or ideal), use corresponding coordinates (real or idealised)

    Atoms whose coordinates are missing in the file ("?" or ".") are skipped.
    """
    d = MMCIF2Dict(str(cif_path))        # whole CIF file as {tag: values}

    # The two coordinate sets use different column names
    if kind == "model":
        cols = [f"_chem_comp_atom.model_Cartn_{a}" for a in "xyz"]
    else:
        cols = [f"_chem_comp_atom.pdbx_model_Cartn_{a}_ideal" for a in "xyz"]

    # Element symbol of every atom (upper-cased so "Fe" and "FE" both work)
    elements = [e.upper() for e in as_list(d["_chem_comp_atom.type_symbol"])]

    # row 1 = xyz[0] is the list of all x values, xyz[1] all y values, xyz[2] all z
    xyz = [as_list(d[c]) for c in cols]

    els, coords = [], []
    for i, el in enumerate(elements):
        vals = [xyz[0][i], xyz[1][i], xyz[2][i]]   # x, y, z of atom number i
        if any(v in ("?", ".") for v in vals):     # no coordinates means skip
            continue
        els.append(el)
        coords.append([float(v) for v in vals])    # convert to number
    return els, np.array(coords)

# ----------------------------------------------------------------------------
# Calculate CM
# ----------------------------------------------------------------------------

def centre_of_mass(elements, coords):
    """
    Calculate the centre of mass and the geometric centre of a set of atoms:
    centre of mass = sum(mass_i * position_i) / sum(mass_i)
    geometric centre (centroid) = sum(position_i) / number_of_atoms

    Return (com, geom), each a numpy array [x, y, z].
    """
    # Mass of every atom, as a column of numbers in the same order as coords.
    masses = np.array([ATOMIC_MASS[e] for e in elements])

    # Centre of mass = sum(mass * position) / total mass.
    com = (coords * masses[:, None]).sum(axis=0) / masses.sum()
    # masses[:, None] turns the mass list into a column so each atom's x, y and z are
    # all multiplied by that atom's mass.
    # .sum(axis=0) adds up over atoms (one total each for x, y and z).

    # Geometric centre
    geom = coords.mean(axis=0)
    return com, geom


def main():
    cif_files = sorted(INPUT_DIR.glob("*.cif"))
    print(f"Reading {len(cif_files)} CIF files")

    # Get readable label from ligand ID
    label_of = {v: k for k, v in LIGANDS.items()}

    rows = []                                      # one result row per calculation
    for cif in cif_files:
        # Step 1. Find out which ligand file contains by ID
        try:
            lig_id = ligand_id_of(cif)
        except Exception as e:
            print(f"{cif.name}: could not read ({e})")
            continue                               # carry on with the next file
        label = label_of.get(lig_id, lig_id)       # unknown IDs keep their ID as label

        # Step 2. Centre calculations for both coordinate sets
        for kind in ("ideal", "model"):
            try:
                els, coords = read_atoms(cif, kind)
            except KeyError:                    
                continue
            if len(els) == 0:                      
                continue

            com, geom = centre_of_mass(els, coords)

            # Step 3. Store the result
            rows.append({
                "file": cif.name,                  # input file the row came from
                "ligand": lig_id,                  # RCSB ID
                "type": label,                     # readable name
                "coords": kind,                    # "ideal" or "model" coordinates
                "n_atoms": len(els),               # atoms used in the calculation
                "formula": "".join(f"{e.capitalize()}{els.count(e)}"
                                   for e in dict.fromkeys(els)),
                # com_*: centre of mass coordinates
                "com_x": round(float(com[0]), 4),
                "com_y": round(float(com[1]), 4),
                "com_z": round(float(com[2]), 4),
                # geo_*: geometric centre coordinates
                "geo_x": round(float(geom[0]), 4),
                "geo_y": round(float(geom[1]), 4),
                "geo_z": round(float(geom[2]), 4),
            })
            print(f"{cif.name}: {lig_id} ({label}) [{kind}] {rows[-1]['formula']}")
            print(f"   CoM:       {np.round(com, 4)}")
            print(f"   Geometric: {np.round(geom, 4)}")

    # Step 4: Save a CSV file
    if rows:
        out = PROJECT_DIR / "Results" / "2_manual" / "Binding_CM" /"cofactor_cm.csv"

        with open(out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=rows[0].keys())
            w.writeheader()
            w.writerows(rows)

        print(f"\nSaved {len(rows)} rows to {out}")


if __name__ == "__main__":
    main()