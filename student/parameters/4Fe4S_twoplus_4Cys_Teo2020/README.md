> **Student package copy.** Atom types were renamed so that every set has unique types: M1 -> 41, M2 -> 42, M3 -> 43, M4 -> 44, Y1 -> 45, Y2 -> 46, Y3 -> 47, Y4 -> 48, Y5 -> 49, Y6 -> 4a, Y7 -> 4b, Y8 -> 4c. The published charges summed to -1.999990; -0.000010 e was spread over the cluster core atoms to make the total an exact integer. Parameters are otherwise identical to the original folder. The tleap snippet below uses the original type names; fes_prep.py takes care of all of this automatically.

# 4Fe4S_twoplus_4Cys_Teo2020

[Fe4S4]2+ cluster with 4 Cys (reduced partner of the HiPIP-like 3+ set; 2Fe(III)/2Fe(II)). Total charge of cluster + 4 Cys = -2.

Source: R. D. Teo, A. Migliore, D. N. Beratan et al., *Chem. Sci.* **2020**, 11, 7076, doi:10.1039/d0sc02245d (B3LYP/6-31G*, Seminario-type, RESP; ff14SB).

## Files

| File | Content |
|---|---|
| `FSJ_core.mol2` | cluster core: Fe and bridging S atoms (types, charges) |
| `C1J.lib` | Cys bound to FE1 via SG |
| `C2J.lib` | Cys bound to FE2 via SG |
| `C3J.lib` | Cys bound to FE3 via SG |
| `C4J.lib` | Cys bound to FE4 via SG |
| `FeS_twoplus.frcmod` | bonded + nonbonded parameters |

## Charge check

| Unit | Charge |
|---|---|
| FSJ (cluster core) | +0.2365 |
| C1J | -0.5546 |
| C2J | -0.5574 |
| C3J | -0.5897 |
| C4J | -0.5348 |
| **total** | **-2.0000** |

## Notes / assumptions (read before use)

- All bonded terms are the "average" columns of the SI tables (S6/S7 bonds, S8/S9 angles) over the optimised geometries; charges are the "average" column of S10. Atom-specific: Teo Fe_i -> FE_i; sulfides Teo S1/S3/S5/S7 -> core S4/S3/S1/S2; thiolates S2/S4/S6/S8 -> SG of C1J..C4J (Cys i on FE_i).
- Charges were fitted on [Fe4S4(SCH3)4] with one methyl H fixed at -0.0987 (= sum of the ff14SB CYM backbone N,H,CA,HA,C,O). So: CB = methyl C, HB2/HB3 = the two free methyl H, SG = thiolate S, and the backbone uses CYM charges (CA -0.0351, HA 0.0508).
- S-CB bond from the SI (S-C ~1.85 A). HB-CB-SG and CA-CB-SG angles: AMBER defaults (Teo's S-C-H methyl angles not used).
- Fe vdW from the main text: interpolated between the Li/Merz IOD Fe2+ and Fe3+ sets for **TIP3P** water, which the paper also used. For OPC, interpolate the OPC IOD values instead (Fe3+ 1.400/0.01570749, Fe2+ 1.391/0.01430674; frcmod.ionslm_iod_opc).
- S vdW not stated; standard values used. No torsions derived (zero).

## tleap usage (sketch)

```
source leaprc.protein.ff19SB   # CA type XC as in your MCPB folders; for ff14SB replace XC by CX in libs+frcmod
addAtomTypes { { "M1" "Fe" "sp3" } { "M2" "Fe" "sp3" } { "M3" "Fe" "sp3" } { "M4" "Fe" "sp3" } { "Y1" "S" "sp3" } { "Y2" "S" "sp3" } { "Y3" "S" "sp3" } { "Y4" "S" "sp3" } { "Y5" "S" "sp3" } { "Y6" "S" "sp3" } { "Y7" "S" "sp3" } { "Y8" "S" "sp3" } }
loadamberparams FeS_twoplus.frcmod
FSJ = loadmol2 FSJ_core.mol2
loadoff C1J.lib
loadoff C2J.lib
loadoff C3J.lib
loadoff C4J.lib
mol = loadpdb protein_renamed.pdb   # rename the coordinating residues/cluster in the PDB first
# <X> = residue number of X in your PDB
bond mol.<C1J>.SG  mol.<FSJ>.FE1
bond mol.<C2J>.SG  mol.<FSJ>.FE2
bond mol.<C3J>.SG  mol.<FSJ>.FE3
bond mol.<C4J>.SG  mol.<FSJ>.FE4
```
