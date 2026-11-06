> **Student package copy.** Atom types were renamed so that every set has unique types: M1 -> 11, M2 -> 12, M3 -> 13, M4 -> 14, Y5 -> 15, Y1 -> 16, Y7 -> 17, Y8 -> 18, Y6 -> 19, Y2 -> 1a, Y3 -> 1b, Y4 -> 1c, Y9 -> 1d. The published charges summed to -1.999839; -0.000161 e was spread over the cluster core atoms to make the total an exact integer. Parameters are otherwise identical to the original folder. The tleap snippet below uses the original type names; fes_prep.py takes care of all of this automatically.

# 4Fe4S_oneplus_3Cys

Own MCPB.py parameters for a [4Fe-4S]+ cluster with 3 Cys; FE4 (type M4) has no protein ligand.
Total charge of cluster + 3 Cys = -2.

**Water model: OPC.** The Fe van der Waals terms (M1-M4: Rmin/2 = 1.400 A, eps = 0.01570749 kcal/mol) are the
Li/Merz IOD set for Fe3+ in OPC water (Li et al., *JCTC* 2021, 17, 2342; AmberTools `frcmod.ionslm_iod_opc`).
They are applied to all four Fe, including the reduced ones. They are only balanced for OPC; do not combine
these files with TIP3P or other water models without replacing the Fe NONB lines with the IOD set for that
water model.

Other nonbonded terms:

- Cys SG and sulfides: gaff2 S (1.9825, 0.2824).
- Y8/Y9 (O types) are MCPB.py leftovers that no atom in these files uses.
