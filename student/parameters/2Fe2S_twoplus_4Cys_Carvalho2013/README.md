> **Student package copy.** Atom types were renamed so that every set has unique types: M1 -> K1, M2 -> K2, Y1 -> K3, Y2 -> K4, Y3 -> K5, Y4 -> K6, Y5 -> K7, Y6 -> K8. Parameters are otherwise identical to the original folder. The tleap snippet below uses the original type names; fes_prep.py takes care of all of this automatically.

# 2Fe2S_twoplus_4Cys_Carvalho2013

Plant-type ferredoxin [2Fe-2S]2+ cluster with 4 Cys (both Fe(III), antiferromagnetically coupled, S=0). Total charge of cluster + 4 Cys = -2.

Source: A. T. P. Carvalho, A. F. S. Teixeira, M. J. Ramos, *J. Comput. Chem.* **2013**, 34, 1540, doi:10.1002/jcc.23287 (B3LYP/6-311++G(2d,2p), Seminario, RESP; ff99SB). Model: *M. laminosus* Fd (PDB 1RFK).

## Files

| File | Content |
|---|---|
| `FSB_core.mol2` | cluster core: Fe and bridging S atoms (types, charges) |
| `C1B.lib` | Cys bound to FE1 via SG |
| `C2B.lib` | Cys bound to FE1 via SG |
| `C3B.lib` | Cys bound to FE2 via SG |
| `C4B.lib` | Cys bound to FE2 via SG |
| `FeS_twoplus.frcmod` | bonded + nonbonded parameters |

## Charge check

| Unit | Charge |
|---|---|
| FSB (cluster core) | +0.1680 |
| C1B | -0.5420 |
| C2B | -0.5420 |
| C3B | -0.5420 |
| C4B | -0.5420 |
| **total** | **-2.0000** |

## Notes / assumptions (read before use)

- Charges from Fig. 3 (model system 3) of the paper; the paper gives one set for all four Cys and one for both Fe / both S. Site total = -2.000 exactly.
- Fig. 3 prints the amide H of the Cys as **+0.0443** and HA as -0.0443 (instead of the usual +0.2719). Kept as printed because the site then sums to exactly -2.000.
- **The Fe-S(Cys)-CB angle is NOT given in Table 1.** A placeholder (27.6 kcal/mol/rad2, 101 deg, the rubredoxin CT-SF-FE value of the same paper) is in the frcmod and marked PLACEHOLDER. Replace it, e.g. with your own Seminario value or ask the authors.
- Bond lengths are only given to one decimal (2.2 A) in the paper.
- The CT-SG bond is not given in the paper; AMBER default (237.0, 1.810) used.
- Coordinates in the mol2/libs are taken from the CrFd1 crystal geometry (Nakayoshi 2022 SI) - only used as templates by tleap.
- Fe vdW: Giammona (ref. 30 in the paper). S vdW not given; ff99SB value used. No torsions (zero).

## tleap usage (sketch)

```
source leaprc.protein.ff19SB   # CA type XC as in your MCPB folders; for ff14SB replace XC by CX in libs+frcmod
addAtomTypes { { "M1" "Fe" "sp3" } { "M2" "Fe" "sp3" } { "Y1" "S" "sp3" } { "Y2" "S" "sp3" } { "Y3" "S" "sp3" } { "Y4" "S" "sp3" } { "Y5" "S" "sp3" } { "Y6" "S" "sp3" } }
loadamberparams FeS_twoplus.frcmod
FSB = loadmol2 FSB_core.mol2
loadoff C1B.lib
loadoff C2B.lib
loadoff C3B.lib
loadoff C4B.lib
mol = loadpdb protein_renamed.pdb   # rename the coordinating residues/cluster in the PDB first
# <X> = residue number of X in your PDB
bond mol.<C1B>.SG  mol.<FSB>.FE1
bond mol.<C2B>.SG  mol.<FSB>.FE1
bond mol.<C3B>.SG  mol.<FSB>.FE2
bond mol.<C4B>.SG  mol.<FSB>.FE2
```
