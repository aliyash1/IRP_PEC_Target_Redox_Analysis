> **Student package copy.** Atom types were renamed so that every set has unique types: M1 -> L1, Y1 -> L2, Y2 -> L3, Y3 -> L4, Y4 -> L5. Parameters are otherwise identical to the original folder. The tleap snippet below uses the original type names; fes_prep.py takes care of all of this automatically.

# 1Fe_threeplus_4Cys_Carvalho2013

Mononuclear Fe(III) site of rubredoxin (Fe(SCys)4, high-spin d5, S=5/2), oxidised. Total charge of Fe + 4 Cys = -1.

Source: A. T. P. Carvalho, A. F. S. Teixeira, M. J. Ramos, *J. Comput. Chem.* **2013**, 34, 1540, doi:10.1002/jcc.23287 (B3LYP/6-311++G(2d,2p), Seminario, RESP; ff99SB).

## Files

| File | Content |
|---|---|
| `FRA_core.mol2` | cluster core: Fe and bridging S atoms (types, charges) |
| `C1A.lib` | Cys bound to FE1 via SG |
| `C2A.lib` | Cys bound to FE1 via SG |
| `C3A.lib` | Cys bound to FE1 via SG |
| `C4A.lib` | Cys bound to FE1 via SG |
| `Fe_rubredoxin.frcmod` | bonded + nonbonded parameters |

## Charge check

| Unit | Charge |
|---|---|
| FRA (cluster core) | +1.0500 |
| C1A | -0.5125 |
| C2A | -0.5125 |
| C3A | -0.5125 |
| C4A | -0.5125 |
| **total** | **-1.0000** |

## Notes / assumptions (read before use)

- Bonded terms (Table 1) and RESP charges (Fig. 3, model system 1) from the paper; all four Cys share identical charges and parameters.
- Table 1 lists the angle as "C-FE-SF" (27.6 kcal/mol/rad2, 101 deg). Physically this can only be the Fe-S-Cbeta angle (the text quotes Fe-S-CT = 101 deg in QM), so it is used as CT-SG-Fe.
- Fig. 3 prints the backbone amide H as **-0.2649** (negative). That is chemically odd, but with it the site sums to exactly -1.000; with +0.2649 it would not. Kept as printed - check the figure if in doubt.
- The CT-SG bond is not given in the paper; the AMBER default used by MCPB.py (237.0, 1.810) is used.
- Fe vdW: Giammona (ref. 30 in the paper, R*=1.20, eps=0.05 - the classic AMBER heme Fe). S vdW not stated; ff99SB value used.
- No torsions were derived (paper sets them to zero).

## tleap usage (sketch)

```
source leaprc.protein.ff19SB   # CA type XC as in your MCPB folders; for ff14SB replace XC by CX in libs+frcmod
addAtomTypes { { "M1" "Fe" "sp3" } { "Y1" "S" "sp3" } { "Y2" "S" "sp3" } { "Y3" "S" "sp3" } { "Y4" "S" "sp3" } }
loadamberparams Fe_rubredoxin.frcmod
FRA = loadmol2 FRA_core.mol2
loadoff C1A.lib
loadoff C2A.lib
loadoff C3A.lib
loadoff C4A.lib
mol = loadpdb protein_renamed.pdb   # rename the coordinating residues/cluster in the PDB first
# <X> = residue number of X in your PDB
bond mol.<C1A>.SG  mol.<FRA>.FE1
bond mol.<C2A>.SG  mol.<FRA>.FE1
bond mol.<C3A>.SG  mol.<FRA>.FE1
bond mol.<C4A>.SG  mol.<FRA>.FE1
```
