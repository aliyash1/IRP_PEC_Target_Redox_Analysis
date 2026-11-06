> **Student package copy.** Atom types were renamed so that every set has unique types: M1 -> 71, M2 -> 72, M3 -> 73, M4 -> 74, Y1 -> 75, Y2 -> 76, Y3 -> 77, Y4 -> 78, Y5 -> 79, Y6 -> 7a, Y7 -> 7b, Y8 -> 7c. Parameters are otherwise identical to the original folder. The tleap snippet below uses the original type names; fes_prep.py takes care of all of this automatically.

# 4Fe4S_twoplus_4Cys_CarvalhoSwart2014

[Fe4S4]2+ cluster with 4 Cys (2Fe(III)/2Fe(II), S=0). Total charge of cluster + 4 Cys = -2.

Source: A. T. P. Carvalho, M. Swart, *J. Chem. Inf. Model.* **2014**, 54, 613, doi:10.1021/ci400718m; Correction *J. Chem. Inf. Model.* **2015**, doi:10.1021/acs.jcim.5b00390 (OPBE/TZ2P DFT, Seminario, MDC-d charges). Also distributed in the Bryce AMBER parameter database ("Fe4S4", NapA).

## Files

| File | Content |
|---|---|
| `FSH_core.mol2` | cluster core: Fe and bridging S atoms (types, charges) |
| `C1H.lib` | Cys bound to FE1 via SG |
| `C2H.lib` | Cys bound to FE2 via SG |
| `C3H.lib` | Cys bound to FE3 via SG |
| `C4H.lib` | Cys bound to FE4 via SG |
| `FeS_twoplus.frcmod` | bonded + nonbonded parameters |

## Charge check

| Unit | Charge |
|---|---|
| FSH (cluster core) | -0.5284 |
| C1H | -0.3679 |
| C2H | -0.3679 |
| C3H | -0.3679 |
| C4H | -0.3679 |
| **total** | **-2.0000** |

## Notes / assumptions (read before use)

- Bonded parameters: SI Table S1 (bonds) and Table S3 (angles) of the 2015 Correction (doi:10.1021/acs.jcim.5b00390), [Fe4S4]2+ block. The cluster is treated as D2d: one Fe pair (FU) with its two shared sulfides (SU) and the other pair (FD/SD). Here FU = FE1/FE2, SU = S3/S4, FD = FE3/FE4, SD = S1/S2 (the FU-SD and FD-SU bonds are the short 2.22 A ones). Which pair is "up" is arbitrary - in a real protein assign it from the crystal structure if the split matters to you.
- Charges: MDC-d charges of Fig. S2 ([Fe4S4(SCH3)4]2-, S=0) averaged over the 4 Fe, 4 sulfides and 4 thiolates (their spread is <0.01 e and the SI does not map the numbering onto the D2d labels). The SCH3 methyl group charge (0.7642 in total) is only given as a sum, so I put 1/4 of it on each CB(+2HB) and kept the backbone at the CYX values your MCPB libs use (HB = 0.0910, CB fitted). Each Cys = -0.3679, which reproduces the Bryce-database implementation of this set.
- Fe vdW: not given in the SI nor in the main text of the paper (checked). The Bryce AMBER parameter database version of this set (contributed by Szaleniec/Swart) uses UFF Fe (R*=1.456, eps=0.013) and S 2.0/0.25; those are used here.
- No torsions derived (zero).

## tleap usage (sketch)

```
source leaprc.protein.ff19SB   # CA type XC as in your MCPB folders; for ff14SB replace XC by CX in libs+frcmod
addAtomTypes { { "M1" "Fe" "sp3" } { "M2" "Fe" "sp3" } { "M3" "Fe" "sp3" } { "M4" "Fe" "sp3" } { "Y1" "S" "sp3" } { "Y2" "S" "sp3" } { "Y3" "S" "sp3" } { "Y4" "S" "sp3" } { "Y5" "S" "sp3" } { "Y6" "S" "sp3" } { "Y7" "S" "sp3" } { "Y8" "S" "sp3" } }
loadamberparams FeS_twoplus.frcmod
FSH = loadmol2 FSH_core.mol2
loadoff C1H.lib
loadoff C2H.lib
loadoff C3H.lib
loadoff C4H.lib
mol = loadpdb protein_renamed.pdb   # rename the coordinating residues/cluster in the PDB first
# <X> = residue number of X in your PDB
bond mol.<C1H>.SG  mol.<FSH>.FE1
bond mol.<C2H>.SG  mol.<FSH>.FE2
bond mol.<C3H>.SG  mol.<FSH>.FE3
bond mol.<C4H>.SG  mol.<FSH>.FE4
```
