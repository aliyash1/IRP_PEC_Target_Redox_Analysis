> **Student package copy.** Atom types were renamed so that every set has unique types: M1 -> E1, M2 -> E2, Y1 -> E3, Y2 -> E4, Y3 -> E5, Y5 -> E6, Y6 -> E7. The published charges summed to -2.000200; +0.000200 e was spread over the cluster core atoms to make the total an exact integer. Parameters are otherwise identical to the original folder. The tleap snippet below uses the original type names; fes_prep.py takes care of all of this automatically.

# 2Fe2S_oneplus_3Cys1His_Zuo2021

mitoNEET / NAF-1-type [2Fe-2S]+ cluster with 3 Cys + 1 His (Ndelta1), **reduced** (Fe(III)/Fe(II)), His neutral (Nepsilon-H). Total charge of cluster + 3 Cys + His = -2.

Source: K. Zuo, H.-B. Marjault, K. L. Bren, G. Rossetti, R. Nechushtai, P. Carloni, *J. Biol. Inorg. Chem.* **2021**, 26, 763, doi:10.1007/s00775-021-01890-8 (MCPB, Seminario, RESP at uB3LYP/6-311++G(2d,2p); AMBER99SB-ILDN in GROMACS units). Bonds: main-text Table 1 (B3LYP, His protonated, reduced); angles, charges, Fe LJ: SI Table S4.

## Files

| File | Content |
|---|---|
| `FSG_core.mol2` | cluster core: Fe and bridging S atoms (types, charges) |
| `C1G.lib` | Cys bound to FE1 via SG |
| `C2G.lib` | Cys bound to FE1 via SG |
| `C3G.lib` | Cys bound to FE2 via SG |
| `H1G.lib` | His bound to FE2 via ND1 |
| `FeS_oneplus.frcmod` | bonded + nonbonded parameters |

## Charge check

| Unit | Charge |
|---|---|
| FSG (cluster core) | -0.1960 |
| C1G | -0.6070 |
| C2G | -0.4665 |
| C3G | -0.6621 |
| H1G | -0.0686 |
| **total** | **-2.0002** |

## Notes / assumptions (read before use)

- Reduced state only. Zuo et al. used the **reduced** bonded set for MD; their oxidised RESP charges (SI Table S4) sum to about -0.69 instead of -1 and no oxidised angles are given, so no oxidised folder was built (see top-level README).
- GROMACS -> AMBER conversion: K_bond(AMBER, kcal/mol/A2) = kb / 4.184 / 100 / 2; K_angle = ktheta / 4.184 / 2; R*/2 = sigma * 2^(1/6) / 2; eps / 4.184. Original values are in the frcmod comments.
- Water model: the converted Fe LJ values equal the Li/Merz IOD sets for **TIP3P** (FeX = Fe3+ 1.386/0.01357, FeY = Fe2+ 1.409/0.01721; AmberTools frcmod.ions234lm_iod_tip3p), consistent with the TIP3P water used in the paper. For OPC, replace them with the OPC IOD values (Fe3+ 1.400/0.01570749, Fe2+ 1.391/0.01430674; frcmod.ionslm_iod_opc).
- FeX (FE1, type M1) carries Cys1 and Cys2; FeY (FE2, M2) carries Cys3 and His (Ndelta1). In mitoNEET: Cys72, Cys74 on FeX; Cys83, His87 on FeY.
- His is Nepsilon-protonated (neutral, HIE-like) and coordinated through ND1. ND1 keeps the standard type NB (as in the Dawson Rieske set), so all His ring terms stay standard and the ring torsions X-CC-NB-X / X-CR-NB-X also act on Fe (keeps Fe in the ring plane). MCPB.py would instead create a new ND1 type with copied NB terms.
- The Cbeta-S bond and the HB/CA-CB-SG angles are not given; AMBER defaults used (as MCPB.py).
- S vdW not given in the paper; standard values used. Torsions involving Fe set to zero (as in the paper).
- Template coordinates are from a Rieske site (PDB 2DE5, Dawson 2016) and only serve as templates.

## tleap usage (sketch)

```
source leaprc.protein.ff19SB   # CA type XC as in your MCPB folders; for ff14SB replace XC by CX in libs+frcmod
addAtomTypes { { "M1" "Fe" "sp3" } { "M2" "Fe" "sp3" } { "Y1" "S" "sp3" } { "Y2" "S" "sp3" } { "Y3" "S" "sp3" } { "Y5" "S" "sp3" } { "Y6" "S" "sp3" } }
loadamberparams FeS_oneplus.frcmod
FSG = loadmol2 FSG_core.mol2
loadoff C1G.lib
loadoff C2G.lib
loadoff C3G.lib
loadoff H1G.lib
mol = loadpdb protein_renamed.pdb   # rename the coordinating residues/cluster in the PDB first
# <X> = residue number of X in your PDB
bond mol.<C1G>.SG  mol.<FSG>.FE1
bond mol.<C2G>.SG  mol.<FSG>.FE1
bond mol.<C3G>.SG  mol.<FSG>.FE2
bond mol.<H1G>.ND1  mol.<FSG>.FE2
```
