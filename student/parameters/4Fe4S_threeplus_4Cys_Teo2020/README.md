> **Student package copy.** Atom types were renamed so that every set has unique types: M1 -> 51, M2 -> 52, M3 -> 53, M4 -> 54, Y1 -> 55, Y2 -> 56, Y3 -> 57, Y4 -> 58, Y5 -> 59, Y6 -> 5a, Y7 -> 5b, Y8 -> 5c. The published charges summed to -1.000010; +0.000010 e was spread over the cluster core atoms to make the total an exact integer. Parameters are otherwise identical to the original folder. The tleap snippet below uses the original type names; fes_prep.py takes care of all of this automatically.

# 4Fe4S_threeplus_4Cys_Teo2020

HiPIP-like [Fe4S4]3+ cluster with 4 Cys (oxidised, 3Fe(III)/1Fe(II)); from the DNA primase p58c. Total charge of cluster + 4 Cys = -1.

Source: R. D. Teo, A. Migliore, D. N. Beratan et al., *Chem. Sci.* **2020**, 11, 7076, doi:10.1039/d0sc02245d (B3LYP/6-31G*, Seminario-type, RESP; ff14SB).

## Files

| File | Content |
|---|---|
| `FSI_core.mol2` | cluster core: Fe and bridging S atoms (types, charges) |
| `C1I.lib` | Cys bound to FE1 via SG |
| `C2I.lib` | Cys bound to FE2 via SG |
| `C3I.lib` | Cys bound to FE3 via SG |
| `C4I.lib` | Cys bound to FE4 via SG |
| `FeS_threeplus.frcmod` | bonded + nonbonded parameters |

## Charge check

| Unit | Charge |
|---|---|
| FSI (cluster core) | +0.5660 |
| C1I | -0.3899 |
| C2I | -0.3890 |
| C3I | -0.4170 |
| C4I | -0.3701 |
| **total** | **-1.0000** |

## Notes / assumptions (read before use)

- All bonded terms are the "average" columns of the SI tables (S1/S2 bonds, S3/S4 angles) over the optimised geometries; charges are the "average" column of S5. Atom-specific: Teo Fe_i -> FE_i; sulfides Teo S1/S3/S5/S7 -> core S4/S3/S1/S2; thiolates S2/S4/S6/S8 -> SG of C1I..C4I (Cys i on FE_i).
- Charges were fitted on [Fe4S4(SCH3)4] with one methyl H fixed at -0.0987 (= sum of the ff14SB CYM backbone N,H,CA,HA,C,O). So: CB = methyl C, HB2/HB3 = the two free methyl H, SG = thiolate S, and the backbone uses CYM charges (CA -0.0351, HA 0.0508).
- S-CB bond from the SI (S-C ~1.85 A). HB-CB-SG and CA-CB-SG angles: AMBER defaults (Teo's S-C-H methyl angles not used).
- Fe vdW from the main text: interpolated between the Li/Merz IOD Fe2+ and Fe3+ sets for **TIP3P** water, which the paper also used. For OPC, interpolate the OPC IOD values instead (Fe3+ 1.400/0.01570749, Fe2+ 1.391/0.01430674; frcmod.ionslm_iod_opc).
- S vdW not stated; standard values used. No torsions derived (zero).
- The 3+ tables (S1-S5) have some shifted cells in the text layer of the PDF; the "average" column (last number of each row) was used. Spot-check against the PDF if you rely on individual angles.

## tleap usage (sketch)

```
source leaprc.protein.ff19SB   # CA type XC as in your MCPB folders; for ff14SB replace XC by CX in libs+frcmod
addAtomTypes { { "M1" "Fe" "sp3" } { "M2" "Fe" "sp3" } { "M3" "Fe" "sp3" } { "M4" "Fe" "sp3" } { "Y1" "S" "sp3" } { "Y2" "S" "sp3" } { "Y3" "S" "sp3" } { "Y4" "S" "sp3" } { "Y5" "S" "sp3" } { "Y6" "S" "sp3" } { "Y7" "S" "sp3" } { "Y8" "S" "sp3" } }
loadamberparams FeS_threeplus.frcmod
FSI = loadmol2 FSI_core.mol2
loadoff C1I.lib
loadoff C2I.lib
loadoff C3I.lib
loadoff C4I.lib
mol = loadpdb protein_renamed.pdb   # rename the coordinating residues/cluster in the PDB first
# <X> = residue number of X in your PDB
bond mol.<C1I>.SG  mol.<FSI>.FE1
bond mol.<C2I>.SG  mol.<FSI>.FE2
bond mol.<C3I>.SG  mol.<FSI>.FE3
bond mol.<C4I>.SG  mol.<FSI>.FE4
```
