# Parameter sets and the force fields they belong to

Every Fe-S parameter set was derived or tested together with a particular protein force field and water
model. The table lists them. The package builds every system with one protein force field and one water
model, set in `config.ini`:

- Default: ff19SB with OPC.
- Alternative: ff14SB with TIP3P.

## Overview

| Parameter set (folder in `parameters/`) | Cluster, ligands, state | Total charge (cluster + ligands) | Published with: protein FF | Water | Fe van der Waals | Source |
|---|---|---|---|---|---|---|
| `4Fe4S_twoplus_3Cys` | [4Fe-4S]2+, 3 Cys, one free Fe | -1 | ff19SB | OPC | Li/Merz IOD Fe3+ (OPC) | Loonen et al., bioRxiv 2026 (MCPB.py) |
| `4Fe4S_oneplus_3Cys` | [4Fe-4S]+, 3 Cys, one free Fe | -2 | ff19SB | OPC | Li/Merz IOD Fe3+ (OPC) | Loonen et al., bioRxiv 2026 (MCPB.py) |
| `4Fe4S_twoplus_4Cys_Teo2020` | [4Fe-4S]2+, 4 Cys | -2 | ff14SB | TIP3P | interpolated Li/Merz IOD (TIP3P) | Teo et al., Chem. Sci. 2020 |
| `4Fe4S_threeplus_4Cys_Teo2020` | [4Fe-4S]3+ (HiPIP), 4 Cys | -1 | ff14SB | TIP3P | interpolated Li/Merz IOD (TIP3P) | Teo et al., Chem. Sci. 2020 |
| `4Fe4S_twoplus_4Cys_CarvalhoSwart2014` | [4Fe-4S]2+, 4 Cys | -2 | ff99SB | not stated | UFF | Carvalho & Swart, JCIM 2014 |
| `4Fe4S_oneplus_4Cys_CarvalhoSwart2014` | [4Fe-4S]+, 4 Cys | -3 | ff99SB | not stated | UFF | Carvalho & Swart, JCIM 2014 |
| `2Fe2S_twoplus_4Cys_Nakayoshi2022` | [2Fe-2S]2+, 4 Cys | -2 | ff14SB | TIP3P | Li/Merz IOD Fe3+ (**OPC**) | Nakayoshi et al., IJMS 2022 (MCPB.py) |
| `2Fe2S_oneplus_4Cys_Nakayoshi2022` | [2Fe-2S]+, 4 Cys | -3 | ff14SB | TIP3P | Li/Merz IOD Fe3+ (**OPC**) | Nakayoshi et al., IJMS 2022 (MCPB.py) |
| `2Fe2S_twoplus_4Cys_Carvalho2013` | [2Fe-2S]2+, 4 Cys | -2 | ff99SB | TIP3P | Giammona heme Fe | Carvalho, Teixeira, Ramos, JCC 2013 |
| `2Fe2S_twoplus_2Cys2His_Dawson2016` | Rieske [2Fe-2S]2+, 2 Cys + 2 His | 0 | ff99SB | TIP3P | Giammona heme Fe | Dawson et al., PLoS One 2016 |
| `2Fe2S_oneplus_2Cys2His_Dawson2016` | Rieske [2Fe-2S]+, 2 Cys + 2 His | -1 | ff99SB | TIP3P | Giammona heme Fe | Dawson et al., PLoS One 2016 |
| `2Fe2S_oneplus_3Cys1His_Zuo2021` | [2Fe-2S]+, 3 Cys + 1 His | -2 | AMBER99SB-ILDN | TIP3P | Li/Merz IOD Fe3+/Fe2+ (TIP3P) | Zuo et al., JBIC 2021 |
| `1Fe_threeplus_4Cys_Carvalho2013` | Fe(III), 4 Cys (rubredoxin) | -1 | ff99SB | TIP3P | Giammona heme Fe | Carvalho, Teixeira, Ramos, JCC 2013 |

`python fes_prep.py sets` prints the same list from the `set.json` files.

## What this means in practice

- **Mixing is unavoidable.** No protein force field or water model matches every set, so the package uses one
  combination for everything. The cluster terms (bonds, angles, charges) do not depend on the water model; only
  the Fe van der Waals terms do, and Fe is buried in the cluster. For energy minimisation the effect is negligible.
- **ff19SB + OPC (default):** matches the `4Fe4S_*_3Cys` sets exactly and the Nakayoshi Fe terms (which are OPC
  values). The ff14SB/ff99SB/TIP3P sets work without changes.
- **ff14SB + TIP3P:** matches Teo, Nakayoshi (protein side), Zuo, Carvalho and Dawson more closely. The
  `4Fe4S_*_3Cys` sets were made for OPC; their Fe terms are then slightly off. fes_prep.py renames the backbone
  CA type of the cluster ligands (XC in ff19SB, CX in ff14SB) automatically.
- **Production MD:** if you go beyond minimisation, consider replacing the Fe NONB lines of a set with the IOD values
  for your water model (AmberTools `dat/leap/parm/frcmod.ionslm_iod_opc`, `frcmod.ions234lm_iod_tip3p`):
  Fe3+ 1.400 / 0.01570749 (OPC), 1.386 / 0.01357097 (TIP3P); Fe2+ 1.391 / 0.01430674 (OPC), 1.409 / 0.01721 (TIP3P).

## Known issues per set

- **`4Fe4S_oneplus_3Cys`:** the RESP charges were fitted to the electrostatic potential of a different broken-symmetry
  state (15 kcal/mol above the optimised one). Per-atom differences are small (at most about 0.08 e), the total charge is
  correct. Fine for minimisation; a corrected charge set may follow.
- **`4Fe4S_*_3Cys`:** the cluster core is somewhat expanded compared to crystal structures (longest Fe-S bonds about
  2.55 A instead of at most 2.3 A). The free Fe must really be unligated: if a structure has a ligand on that Fe
  (water, substrate, Glu/Asp), the ligand is removed or kept non-bonded (see README).
- **`4Fe4S_oneplus_4Cys_CarvalhoSwart2014`:** the paper does not say which atoms carry which type; the mapping was inferred.
- **`2Fe2S_twoplus_4Cys_Carvalho2013`:** the Fe-S(Cys)-CB angle is not in the paper; a placeholder from the same paper is used.
- **`2Fe2S_oneplus_3Cys1His_Zuo2021`:** reduced state only; parameters were converted from GROMACS units.
- **Sets derived on one protein:** charges and bonded terms are atom-specific (e.g. Nakayoshi Cys37/42/45/75, the
  own IspG sets). In another protein the Cys are assigned to these templates by topology; which Cys gets which
  template is arbitrary. This is the standard compromise when reusing MCPB.py parameters.

Each set folder has a README with the full source, charge check and assumptions. `BIBLIOGRAPHY.md` and
`FeS_parameters.bib` list all literature.
