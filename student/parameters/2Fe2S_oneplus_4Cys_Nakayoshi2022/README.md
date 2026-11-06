> **Student package copy.** Atom types were renamed so that every set has unique types: M1 -> W1, M2 -> W2, Y1 -> W3, Y2 -> W4, Y5 -> W5, Y6 -> W6, Y3 -> W7, Y4 -> W8. Parameters are otherwise identical to the original folder. The tleap snippet below uses the original type names; fes_prep.py takes care of all of this automatically.

# 2Fe2S_oneplus_4Cys_Nakayoshi2022

Plant-type ferredoxin [2Fe-2S] cluster with 4 Cys, reduced ([2Fe-2S]+, Fe(III)/Fe(II)); total charge -3. In the paper FE1 (Cys37/Cys42 side) is the reduced Fe(II).

Source: T. Nakayoshi, Y. Ohnishi, H. Tanaka, G. Kurisu et al., *Int. J. Mol. Sci.* **2022**, 23, 15913, doi:10.3390/ijms232415913 (MCPB.py; Seminario at LC-wPBE; RESP at HF/6-31G(d); ff14SB). Model: *C. reinhardtii* Fd1 (PDB 6LK1). Files converted **verbatim** from SI Table S6 (their mol2 + frcmod), which is literally MCPB.py output - the closest match to your own workflow.

## Files

| File | Content |
|---|---|
| `FSD_core.mol2` | FE1 (M1), FE2 (M2), bridging S1 (Y3), S2 (Y4) |
| `C1D.lib` | Cys37, SG type Y1, bound to FE1 |
| `C2D.lib` | Cys42, SG type Y2, bound to FE1 |
| `C3D.lib` | Cys45, SG type Y5, bound to FE2 |
| `C4D.lib` | Cys75, SG type Y6, bound to FE2 |
| `FeS_oneplus.frcmod` | MCPB.py frcmod from the SI |

## Charge check

| Unit | Charge |
|---|---|
| FSD (cluster core) | +0.1272 |
| C1D (Cys37, SG type Y1) | -0.9441 |
| C2D (Cys42, SG type Y2) | -0.9212 |
| C3D (Cys45, SG type Y5) | -0.5948 |
| C4D (Cys75, SG type Y6) | -0.6670 |
| **total** | **-3.0000** |

## Notes

- CA type CX (ff14SB) was renamed XC (ff19SB) to match your folders; rename back for ff14SB.
- Charges were fitted on the CrFd1 site; residue-specific (Cys37/42/45/75). For another ferredoxin map the Cys by topology: two Cys on Fe1 (Y1, Y2) and two on Fe2 (Y5, Y6).
- The two bridging sulfides carry identical charges and the backbone N/H/C/O keep the ff values (MCPB restraint).
- Fe vdW: IOD Fe3+ set for OPC water (Li et al. JCTC 2021) for both Fe, also in the reduced state (as in the SI).
- **Water-model mismatch in the original work:** these Fe LJ values are fitted for OPC water, but the paper solvated the protein in TIP3P. For TIP3P the matching IOD values are Fe3+ 1.386 / 0.01357097 and Fe2+ 1.409 / 0.01721 (AmberTools `frcmod.ions234lm_iod_tip3p`); for OPC keep the values as they are.

## tleap usage (sketch)

```
source leaprc.protein.ff19SB
addAtomTypes { { "M1" "Fe" "sp3" } { "M2" "Fe" "sp3" } { "Y1" "S" "sp3" } { "Y2" "S" "sp3" } { "Y3" "S" "sp3" } { "Y4" "S" "sp3" } { "Y5" "S" "sp3" } { "Y6" "S" "sp3" } }
loadamberparams FeS_oneplus.frcmod
FSD = loadmol2 FSD_core.mol2
loadoff C1D.lib
loadoff C2D.lib
loadoff C3D.lib
loadoff C4D.lib
mol = loadpdb protein_renamed.pdb
bond mol.<C1D>.SG mol.<FSD>.FE1
bond mol.<C2D>.SG mol.<FSD>.FE1
bond mol.<C3D>.SG mol.<FSD>.FE2
bond mol.<C4D>.SG mol.<FSD>.FE2
```
