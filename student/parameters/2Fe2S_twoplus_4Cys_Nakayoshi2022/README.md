> **Student package copy.** Atom types were renamed so that every set has unique types: M1 -> J1, M2 -> J2, Y1 -> J3, Y2 -> J4, Y5 -> J5, Y6 -> J6, Y3 -> J7, Y4 -> J8. Parameters are otherwise identical to the original folder. The tleap snippet below uses the original type names; fes_prep.py takes care of all of this automatically.

# 2Fe2S_twoplus_4Cys_Nakayoshi2022

Plant-type ferredoxin [2Fe-2S] cluster with 4 Cys, oxidised ([2Fe-2S]2+, Fe(III)/Fe(III)); total charge -2.

Source: T. Nakayoshi, Y. Ohnishi, H. Tanaka, G. Kurisu et al., *Int. J. Mol. Sci.* **2022**, 23, 15913, doi:10.3390/ijms232415913 (MCPB.py; Seminario at LC-wPBE; RESP at HF/6-31G(d); ff14SB). Model: *C. reinhardtii* Fd1 (PDB 6LK1). Files converted **verbatim** from SI Table S5 (their mol2 + frcmod), which is literally MCPB.py output - the closest match to your own workflow.

## Files

| File | Content |
|---|---|
| `FSC_core.mol2` | FE1 (M1), FE2 (M2), bridging S1 (Y3), S2 (Y4) |
| `C1C.lib` | Cys37, SG type Y1, bound to FE1 |
| `C2C.lib` | Cys42, SG type Y2, bound to FE1 |
| `C3C.lib` | Cys45, SG type Y5, bound to FE2 |
| `C4C.lib` | Cys75, SG type Y6, bound to FE2 |
| `FeS_twoplus.frcmod` | MCPB.py frcmod from the SI |

## Charge check

| Unit | Charge |
|---|---|
| FSC (cluster core) | +0.3615 |
| C1C (Cys37, SG type Y1) | -0.6136 |
| C2C (Cys42, SG type Y2) | -0.6106 |
| C3C (Cys45, SG type Y5) | -0.5496 |
| C4C (Cys75, SG type Y6) | -0.5877 |
| **total** | **-2.0000** |

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
loadamberparams FeS_twoplus.frcmod
FSC = loadmol2 FSC_core.mol2
loadoff C1C.lib
loadoff C2C.lib
loadoff C3C.lib
loadoff C4C.lib
mol = loadpdb protein_renamed.pdb
bond mol.<C1C>.SG mol.<FSC>.FE1
bond mol.<C2C>.SG mol.<FSC>.FE1
bond mol.<C3C>.SG mol.<FSC>.FE2
bond mol.<C4C>.SG mol.<FSC>.FE2
```
