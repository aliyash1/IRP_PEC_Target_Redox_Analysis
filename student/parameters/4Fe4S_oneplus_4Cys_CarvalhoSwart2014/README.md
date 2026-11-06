> **Student package copy.** Atom types were renamed so that every set has unique types: M1 -> 91, M2 -> 92, M3 -> 93, M4 -> 94, Y1 -> 95, Y2 -> 96, Y3 -> 97, Y4 -> 98, Y5 -> 99, Y6 -> 9a, Y7 -> 9b, Y8 -> 9c. Parameters are otherwise identical to the original folder. The tleap snippet below uses the original type names; fes_prep.py takes care of all of this automatically.

# 4Fe4S_oneplus_4Cys_CarvalhoSwart2014

[Fe4S4]+ cluster with 4 Cys (reduced, 1Fe(III)/3Fe(II), S=1/2). Total charge of cluster + 4 Cys = -3.

Source: A. T. P. Carvalho, M. Swart, *J. Chem. Inf. Model.* **2014**, 54, 613, doi:10.1021/ci400718m; Correction *J. Chem. Inf. Model.* **2015**, doi:10.1021/acs.jcim.5b00390 (OPBE/TZP, COSMO eps=4, Seminario, MDC-d charges).

## Files

| File | Content |
|---|---|
| `FSK_core.mol2` | cluster core: Fe and bridging S atoms (types, charges) |
| `C1K.lib` | Cys bound to FE1 via SG |
| `C2K.lib` | Cys bound to FE2 via SG |
| `C3K.lib` | Cys bound to FE3 via SG |
| `C4K.lib` | Cys bound to FE4 via SG |
| `FeS_oneplus.frcmod` | bonded + nonbonded parameters |

## Charge check

| Unit | Charge |
|---|---|
| FSK (cluster core) | -0.9539 |
| C1K | -0.5115 |
| C2K | -0.5115 |
| C3K | -0.5115 |
| C4K | -0.5115 |
| **total** | **-3.0000** |

## Notes / assumptions (read before use)

- Bonded parameters: [Fe4S4]+ block of SI Tables S1/S3 of the 2015 Correction (doi:10.1021/acs.jcim.5b00390). Types F4, FY (Fe) and S3, SY (sulfides) are only defined as "Fe in [Fe4S4]+" / "S connected to three irons" in Table 3 of the paper; the atom mapping is not given anywhere.
- Mapping used (the only one for which every angle that occurs in the cube has a parameter in Table S3): the F4 pair (FE1, FE2) shares the two SY sulfides (S3, S4), the FY pair (FE3, FE4) shares the two S3-type sulfides (S1, S2), and the four FY-SY cross bonds are the short ones (76.4 kcal/mol/A2, 2.23 A). This mirrors the D2d split of the 2+ set (FU/SU, FD/SD). All Fe-S-Fe and S-Fe-S angle parameters are identical across types, so only the two short bonds depend on this choice. Which Fe pair is "FY" in a real protein is arbitrary - assign it from the crystal structure if it matters.
- Table S3 also lists S3-F4-S3 and F4-S3-F4, which do not occur with this mapping; they have the same values as the angles that do occur.
- Charges: MDC-d charges of Fig. S2 ([Fe4S4(SCH3)4]3-, S=1/2) averaged over the 4 Fe (0.234-0.282), 4 sulfides and 4 thiolates, because the numbered atoms cannot be mapped onto the types. In the paper (Table 2) Fe1/Fe2 and Fe3/Fe4 form the two antiparallel spin pairs. The methyl charge (0.6650 in total) was split equally over CB(+2HB) with the CYX backbone of your MCPB libs.
- Fe vdW: not given in the paper or SI (checked both). UFF values as in the Bryce implementation of the 2+ set are used; S 2.0/0.25.
- No torsions derived (zero).

## tleap usage (sketch)

```
source leaprc.protein.ff19SB   # CA type XC as in your MCPB folders; for ff14SB replace XC by CX in libs+frcmod
addAtomTypes { { "M1" "Fe" "sp3" } { "M2" "Fe" "sp3" } { "M3" "Fe" "sp3" } { "M4" "Fe" "sp3" } { "Y1" "S" "sp3" } { "Y2" "S" "sp3" } { "Y3" "S" "sp3" } { "Y4" "S" "sp3" } { "Y5" "S" "sp3" } { "Y6" "S" "sp3" } { "Y7" "S" "sp3" } { "Y8" "S" "sp3" } }
loadamberparams FeS_oneplus.frcmod
FSK = loadmol2 FSK_core.mol2
loadoff C1K.lib
loadoff C2K.lib
loadoff C3K.lib
loadoff C4K.lib
mol = loadpdb protein_renamed.pdb   # rename the coordinating residues/cluster in the PDB first
# <X> = residue number of X in your PDB
bond mol.<C1K>.SG  mol.<FSK>.FE1
bond mol.<C2K>.SG  mol.<FSK>.FE2
bond mol.<C3K>.SG  mol.<FSK>.FE3
bond mol.<C4K>.SG  mol.<FSK>.FE4
```
