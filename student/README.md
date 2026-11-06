# Energy minimisation of Fe-S proteins in GROMACS

This folder turns PDB structures that contain iron-sulfur clusters into GROMACS systems and energy-minimises
them. It handles hundreds of structures in one go. For every cluster, the script `fes_prep.py`:

1. recognises the cluster type and its ligands;
2. picks a literature parameter set;
3. renames the cluster and its Cys/His ligands to match that set;
4. builds the system with AmberTools (tleap), with the protein in ff19SB and solvated in OPC water;
5. converts it to GROMACS with ParmEd.

`run_em.sh` then minimises all systems, and `fes_prep.py report` collects the results in one table.

## Contents

| Path | What it is |
|---|---|
| `fes_prep.py` | preparation (`prep`), results table (`report`), list of parameter sets (`sets`) |
| `config.ini` | settings: force field, water, which parameter set per cluster type |
| `overrides.csv` | per-cluster exceptions, e.g. a reduced cluster in one structure |
| `run_em.sh` | minimises every prepared system (local machine) |
| `em_array.slurm` | the same as a Slurm array job (cluster) |
| `minim.mdp` | GROMACS minimisation settings |
| `parameters/` | 13 parameter sets, one folder each (see `FORCE_FIELDS.md`) |
| `FORCE_FIELDS.md` | which protein force field and water each set belongs to, and known issues |
| `BIBLIOGRAPHY.md`, `FeS_parameters.bib` | literature on Fe-S parameters |
| `example/` | seven test structures and the expected results |

## Installation

You need two conda environments, one for preparation and one for GROMACS. On Windows, use WSL (Ubuntu).

```bash
conda create -n ambertools -c conda-forge ambertools      # tleap, ParmEd, numpy
conda create -n gromacs    -c conda-forge "gromacs=2026.3=nompi*"
```

Tested with AmberTools 25 (Python 3.10, ParmEd 4.3) and GROMACS 2026.3 (and 2024.3).

## Quick start

```bash
conda activate ambertools
python fes_prep.py prep my_pdbs/*.pdb --out systems     # one folder per structure
conda activate gromacs
bash run_em.sh systems                                  # minimise everything
conda activate ambertools
python fes_prep.py report systems                       # writes systems/summary.csv
```

Each `systems/<name>/` folder contains:

| File | Content |
|---|---|
| `prep_report.json` | clusters found, parameter set used, atom mapping, warnings |
| `clean.pdb` | the structure as given to tleap |
| `leap.in`, `leap.log` | tleap input and log |
| `system.top`, `system.gro` | the solvated GROMACS system |
| `em.gro`, `em.log` | minimised structure and log |

Try it first on the example: `python fes_prep.py prep example/inputs/*.pdb --out test`.

## What the preparation does

1. **Cleaning.**
    - Keeps the first model and alternate location A.
    - Removes hydrogens (tleap adds them back) and waters.
    - Converts selenomethionine (MSE) to Met.
    - Keeps the ions listed in `config.ini`; everything else that is not protein or an Fe-S cluster is removed, with a warning.
    - Chain gaps (C–N distance > 2 Å) become chain ends.
    - Disulfides (SG–SG < 2.5 Å) are bonded.
    - Histidines get the name in `config.ini` (default HIE).
2. **Cluster detection.**
    - A cluster is a residue made of only Fe and S atoms (SF4, FES, FE, ...).
    - Ligands are Cys SG within 2.9 Å of an Fe and His N within 2.6 Å.
    - The environment is written as `<cluster>_<ligands>`, e.g. `4Fe4S_4Cys`, `4Fe4S_3Cys`, `2Fe2S_2Cys2His`, `1Fe_4Cys`.
3. **Parameter set.**
    - The default for each environment comes from `[parameter_sets]` in `config.ini`.
    - A cluster listed in `overrides.csv` gets the set named there.
4. **Mapping.**
    - The Fe and S atoms are matched to the set by topology: which Fe binds which S, and which Fe binds Cys or His.
    - The ligating residues are renamed to the set's residue names (for example C1J or H1E).
    - tleap gets explicit Fe–SG / Fe–ND1 bonds.
5. **tleap and ParmEd.**
    - The system is solvated in a rectangular box with a 10 Å margin and neutralised with Na+/Cl−.
    - tleap writes Amber files, which ParmEd converts to `system.top` / `system.gro`.

Every parameter set has its own atom types, so a protein with several cluster types (e.g. [2Fe-2S] and [4Fe-4S]) works.

## Choosing parameter sets and oxidation states

The oxidation state of a cluster cannot be read from a structure, so the choice is yours.

| Environment | Default set (state) | Alternatives |
|---|---|---|
| `4Fe4S_4Cys` | `4Fe4S_twoplus_4Cys_Teo2020` (2+) | `4Fe4S_twoplus_4Cys_CarvalhoSwart2014` (2+), `4Fe4S_oneplus_4Cys_CarvalhoSwart2014` (1+), `4Fe4S_threeplus_4Cys_Teo2020` (3+, HiPIP) |
| `4Fe4S_3Cys` (one free Fe) | `4Fe4S_twoplus_3Cys` (2+) | `4Fe4S_oneplus_3Cys` (1+) |
| `2Fe2S_4Cys` | `2Fe2S_twoplus_4Cys_Nakayoshi2022` (2+) | `2Fe2S_oneplus_4Cys_Nakayoshi2022` (1+), `2Fe2S_twoplus_4Cys_Carvalho2013` (2+) |
| `2Fe2S_2Cys2His` (Rieske) | `2Fe2S_twoplus_2Cys2His_Dawson2016` (2+) | `2Fe2S_oneplus_2Cys2His_Dawson2016` (1+) |
| `2Fe2S_3Cys1His` (mitoNEET) | `2Fe2S_oneplus_3Cys1His_Zuo2021` (1+) | none |
| `1Fe_4Cys` (rubredoxin) | `1Fe_threeplus_4Cys_Carvalho2013` (Fe3+) | none |

- **Change a default for all structures:** edit `[parameter_sets]` in `config.ini`.
- **Change one cluster:** add a line to `overrides.csv`, e.g. `4g9p,A,501,4Fe4S_oneplus_3Cys,`. Chain and residue number are those of the cluster in the original PDB; `*` means all clusters of that structure. The last column, `extra_ligands`, can be `unbound` (see "When a structure fails").
- **List all sets:** `python fes_prep.py sets`.

## Checking the results

Look at `systems/summary.csv` and, for every structure with warnings, at its `prep_report.json`.

| Column | What to look for |
|---|---|
| `prep` | `prepared`; `failed` comes with a reason in `error` |
| `warnings` | number of warnings in `prep_report.json`: removed ligands, chain breaks, bound ligands that were removed |
| `em_status` | `converged to Fmax < 1000 in N steps`; anything else needs a look |
| `fe_s_min_after_em`, `fe_s_max_after_em` | Fe–S(sulfide) bond range after minimisation; about 2.1–2.4 Å is normal (the `4Fe4S_*_3Cys` sets have a few bonds up to about 2.6 Å by construction) |
| `net_charge` | 0 |

Always look at a few minimised structures in VMD/PyMOL/ChimeraX, especially around the clusters.

## When a structure fails

| Message in `error` / warning | Meaning | What to do |
|---|---|---|
| `no parameter set for 4Fe4S_3Cys+other`, `... an Fe also binds another protein atom` | an Fe binds e.g. Glu/Asp/Ser O in addition to Cys | if that bond may be left out, add `<structure>,*,*,,unbound` to `overrides.csv` (the ligand stays, non-bonded; example: 2Y0F, where Glu350 releases the Fe during catalysis). `treat_extra_ligands_as_unbound = yes` in `config.ini` does this for all structures. Otherwise skip the structure. |
| `no parameter set for 3Fe4S_3Cys` or another environment | no parameters exist for this cluster type ([3Fe-4S], Cys3His [4Fe-4S], ...) | skip, or ask your supervisor |
| `removed ... was bound to FE..` | a water/substrate/cofactor on a cluster Fe was removed | fine for minimisation of the apo site; note it |
| `removed non-standard residue/ligand ...` | a ligand/cofactor (FAD, heme, substrate) was removed | if it matters, it needs its own parameters (GAFF2/MCPB.py) |
| `cluster ligand ... is a chain terminus` | a ligating Cys is the first/last residue of a chain or sits next to a gap | not supported by the residue libraries |
| `His ... binds through NE2` | a His ligand uses NE2, the set expects ND1 | not supported |
| `no consistent mapping ...` | the Fe–S distances in the PDB do not form a complete cluster (low resolution, missing atoms) | check the cluster in the PDB |
| `tleap failed: ...` | usually a non-standard residue in the protein chain | see `leap.log`; rename or remove that residue |
| `chain break between ...` (warning) | missing residues; both ends become charged termini | fine for minimisation |

Input must be PDB format. Convert mmCIF first, e.g. `gemmi convert x.cif x.pdb`.

## Running many structures on a cluster

Preparation takes a few seconds per structure on one core, so run it on a login or interactive node. For
the minimisations, list the systems and submit an array job (adapt modules and partition):

```bash
ls -d systems/*/ > systems.txt
mkdir -p logs
sbatch --array=1-$(wc -l < systems.txt)%50 em_array.slurm
```

A rerun of `run_em.sh` skips systems that already have `em.gro` (`FORCE=1` redoes them).

## Points to keep in mind

- **Force fields:** the package uses ff19SB + OPC for everything. Most literature sets were made with ff14SB/ff99SB and
  TIP3P. For minimisation this does not matter; see `FORCE_FIELDS.md` before doing MD.
- **Protonation:** histidines are HIE unless they bind a cluster, and other residues get their standard states.
  If protonation matters, assign it beforehand (e.g. with PROPKA/pdb2pqr) and use HID/HIE/HIP/ASH/GLH names in the PDB.
- **Several Cys on one Fe:** when a set has residue-specific charges (e.g. Nakayoshi, the IspG sets), which Cys of
  a new protein gets which template is arbitrary. This is the usual compromise when reusing such parameters.
- **Removed cofactors:** removing cofactors and substrates changes the structure around the cluster. Report what was removed.
- **AlphaFold models:** fes_prep.py only finds clusters that are present in the file. AlphaFold 2 models have no
  cofactors. Either predict with AlphaFold 3 including the clusters as ligands (CCD codes SF4, FES, F3S, FE; convert
  the mmCIF to PDB), or copy the cluster from a homologous crystal structure after superposing the coordinating Cys.
  Check afterwards that every Cys SG is within 2.9 Å of its Fe. Otherwise the cluster is detected with the wrong ligands,
  or the script stops with "no consistent mapping". Low-pLDDT tails are kept as they are; for minimisation that is fine.
- **Unmodelled residues:** crystal structures often lack tails and flexible loops. fes_prep.py does **not** rebuild them:
    - missing N-/C-terminal tails are left out, and the first/last modelled residue becomes a charged terminus;
    - internal gaps become chain breaks with charged ends;
    - missing side-chain atoms are rebuilt by tleap from its templates.

  `prep_report.json` lists the unmodelled residues per chain (from REMARK 465/470 of the PDB). `summary.csv` has the
  counts (`missing_internal`, `missing_termini`, `residues_missing_atoms`). For minimisation of the cluster region
  this is usually acceptable. If gaps sit near a cluster, or you need complete chains, model them first and give
  the completed PDB to fes_prep.py, e.g. with MODELLER, PDBFixer (OpenMM) or an AlphaFold model aligned onto the
  crystal structure. Keep the Fe-S cluster records from the crystal structure.
