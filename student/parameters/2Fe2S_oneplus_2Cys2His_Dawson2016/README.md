> **Student package copy.** Atom types were renamed so that every set has unique types: FE -> D1, Sb -> D2. The published charges summed to -0.999990; -0.000010 e was spread over the cluster core atoms to make the total an exact integer. Parameters are otherwise identical to the original folder. The tleap snippet below uses the original type names; fes_prep.py takes care of all of this automatically.

# 2Fe2S_oneplus_2Cys2His_Dawson2016

Rieske [2Fe-2S] cluster with 2 Cys (on FE1) + 2 His (Ndelta1, on FE2), His neutral (Nepsilon-H), reduced ([2Fe-2S]+, the His-side Fe2 is Fe(II); total charge -1).

Source: W. K. Dawson, R. Jono, T. Terada, K. Shimizu, *PLoS One* **2016**, 11, e0162031, doi:10.1371/journal.pone.0162031 - carbazole 1,9a-dioxygenase / ferredoxin (PDB 2DE5). Files taken from their S4 File (`amber_qm1S4.prepin`, `oxygenase.frcmod`).

**This set is not MCPB-style.** It was hand-built ("estimated from QC calc and experiment") with generic atom types: FE (both Fe), Sb (bridging S), SH (Cys SG), His ND1 keeps the standard type NB. It includes real (non-zero) torsions and a ring improper. Their residues were renamed for uniqueness: FES -> FSF, CYR -> C1F/C2F (identical charges), HR1 -> H1F (His71), HR2 -> H2F (His93, different charges).

## Files

| File | Content |
|---|---|
| `FSF_core.mol2` | FE1, S1, S2, FE2 (types FE/Sb) |
| `C1F.lib`, `C2F.lib` | Cys bound to FE1 (SG type SH) |
| `H1F.lib`, `H2F.lib` | His bound to FE2 via ND1 (type NB) |
| `FeS_oneplus.frcmod` | their frcmod, Rieske lines only |

## Charge check

| Unit | Charge |
|---|---|
| FSF (FES: FE1 = Cys side, FE2 = His side) | -0.0767 |
| C1F (CYR, 2DE5 Cys69) | -0.5469 |
| C2F (CYR, 2DE5 Cys90) | -0.5469 |
| H1F (HR1, 2DE5 His71) | +0.1905 |
| H2F (HR2, 2DE5 His93) | -0.0201 |
| **total** | **-1.0000** |

## Notes

- Atom types are ff99SB-style, as in the original: CA = CT, Cys/His HB = HC. That works with ff14SB/ff19SB as well, but these residues get no CMAP/backbone-specific terms.
- The type **FE** also exists in the standard parm10 (heme iron, mass 55.00); loading this frcmod changes its mass/vdW for the whole system. If you also have heme, rename FE (e.g. to M1/M2) in mol2 and frcmod. **SH** is the ff Cys-thiol type; vdW identical, so harmless.
- The frcmod has a FE-FE bond term, but their tleap script does not create an Fe-Fe bond, so it is unused (same for the dihedral lines with the undefined type "SF", a typo in the original).
- Fe vdW: 1.20 / 0.05 ("taken from heme"); S 2.0 / 0.25.
- The S4 File contains a second reduced charge set, `amber_qm1S6.prepin` (another spin coupling; FE1 0.79968, S -0.77958, FE2 0.65086). The paper does not say which was used for the production MD; qm1S4 is used here. The bonded terms are the same for both redox states.
- For the separate MCPB alternative of an oxidised Rieske site, see Pabis et al. JCTC 2014 (top-level README).

## tleap usage (sketch)

```
source leaprc.protein.ff14SB   # original work: ff99SB
addAtomTypes { { "FE" "Fe" "sp3" } { "Sb" "S" "sp3" } }
loadamberparams FeS_oneplus.frcmod
FSF = loadmol2 FSF_core.mol2
loadoff C1F.lib
loadoff C2F.lib
loadoff H1F.lib
loadoff H2F.lib
mol = loadpdb protein_renamed.pdb
bond mol.<C1F>.SG  mol.<FSF>.FE1
bond mol.<C2F>.SG  mol.<FSF>.FE1
bond mol.<H1F>.ND1 mol.<FSF>.FE2
bond mol.<H2F>.ND1 mol.<FSF>.FE2
```
