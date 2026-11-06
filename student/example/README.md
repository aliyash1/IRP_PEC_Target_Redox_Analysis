# Example: seven structures, one of each cluster environment

```bash
conda activate ambertools
python ../fes_prep.py prep inputs/*.pdb --out systems
conda activate gromacs
bash ../run_em.sh systems
conda activate ambertools
python ../fes_prep.py report systems
```

| PDB | Protein | Cluster environment | Parameter set (default config) |
|---|---|---|---|
| [1FDN](https://www.rcsb.org/structure/1FDN) | 2[4Fe-4S] ferredoxin, *C. acidurici* | 2 × `4Fe4S_4Cys` | `4Fe4S_twoplus_4Cys_Teo2020` |
| [4G9P](https://www.rcsb.org/structure/4G9P) | IspG with MEcPP, *T. thermophilus* | `4Fe4S_3Cys` (free Fe) | `4Fe4S_twoplus_3Cys` |
| [2Y0F](https://www.rcsb.org/structure/2Y0F) | IspG, *T. thermophilus*, 4 chains | `4Fe4S_3Cys` + Glu350 on the fourth Fe | `4Fe4S_twoplus_3Cys`, Glu350 non-bonded (`overrides.csv`) |
| [1RFK](https://www.rcsb.org/structure/1RFK) | plant-type ferredoxin, *M. laminosus*, 2 chains | 2 × `2Fe2S_4Cys` | `2Fe2S_twoplus_4Cys_Nakayoshi2022` |
| [1RIE](https://www.rcsb.org/structure/1RIE) | Rieske protein, bovine | `2Fe2S_2Cys2His` | `2Fe2S_twoplus_2Cys2His_Dawson2016` |
| [2QH7](https://www.rcsb.org/structure/2QH7) | mitoNEET, human, 2 chains | 2 × `2Fe2S_3Cys1His` | `2Fe2S_oneplus_3Cys1His_Zuo2021` |
| [1IRO](https://www.rcsb.org/structure/1IRO) | rubredoxin, *C. pasteurianum* | `1Fe_4Cys` | `1Fe_threeplus_4Cys_Carvalho2013` |

## Expected results (AmberTools 25, GROMACS 2026.3, ff19SB + OPC)

| PDB | Atoms | EM | Fe–S(sulfide) after EM (Å) | Warnings |
|---|---|---|---|---|
| 1FDN | 13 155 | converged, 374 steps | 2.21–2.33 | none |
| 4G9P | 78 175 | converged, 1249 steps | 2.25–2.51 | MEcPP (CDI) and MES removed |
| 2Y0F | 249 300 | converged, 2843 steps | 2.26–2.50 | Glu350 kept non-bonded (4 chains), chain breaks |
| 1RFK | 39 810 | converged, 1075 steps | 2.12–2.22 | chain break A:9–A:13 |
| 1RIE | 27 516 | converged, 331 steps | 2.15–2.19 | disulfide A:144–A:160 |
| 2QH7 | 29 983 | converged, 432 steps | 2.14–2.28 | none |
| 1IRO | 13 167 | converged, 389 steps | (no sulfides) | none |

Unmodelled residues in the PDB files (not rebuilt; listed in `prep_report.json` and counted in `summary.csv`):
1RFK 3 internal (A:10–12) and 1 terminal; 2QH7 22 terminal and 7 residues with missing atoms; 2Y0F 30 internal
(loops of 9–11 residues in chains B–D), 12 terminal and 10 residues with missing atoms; 1IRO, 1RIE and 4G9P 1–2 terminal.

Step counts and energies differ slightly between runs and machines. The logs of this run (`prep_report.json`,
`leap.in`, `leap.log`, `grompp.log`, `em.log` per structure, and `summary.csv`) are in `expected_output/`.

- **2Y0F:** Glu350 binds the fourth Fe of the substrate-free enzyme. That bond is released on substrate binding and
  during the catalytic cycle, so `overrides.csv` contains `2y0f,*,*,,unbound`. The 3-Cys set is used and Glu350
  stays non-bonded, interacting with the cluster only through electrostatics and van der Waals. Without that line the
  script stops for 2Y0F with an explanation.
- **ff14SB + TIP3P:** the same run with `protein = ff14SB` and `water = tip3p` was also tested (1FDN, 1RIE, 4G9P;
  all converged).
