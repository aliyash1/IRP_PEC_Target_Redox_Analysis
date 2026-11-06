> **Student package copy.** Atom types were renamed so that every set has unique types: M1 -> 01, M2 -> 02, M3 -> 03, M4 -> 04, Y5 -> 05, Y1 -> 06, Y7 -> 07, Y8 -> 08, Y6 -> 09, Y2 -> 0a, Y3 -> 0b, Y4 -> 0c, Y9 -> 0d. The published charges summed to -0.999854; -0.000146 e was spread over the cluster core atoms to make the total an exact integer. Parameters are otherwise identical to the original folder. The tleap snippet below uses the original type names; fes_prep.py takes care of all of this automatically.

# 4Fe4S_twoplus_3Cys

Own MCPB.py parameters for a [4Fe-4S]2+ cluster with 3 Cys; FE4 (type M4) has no protein ligand.
Total charge of cluster + 3 Cys = -1.

**Water model: OPC.** The Fe van der Waals terms (M1-M4: Rmin/2 = 1.400 A, eps = 0.01570749 kcal/mol) are the
Li/Merz IOD set for Fe3+ in OPC water (Li et al., *JCTC* 2021, 17, 2342; AmberTools `frcmod.ionslm_iod_opc`).
They are only balanced for OPC; do not combine these files with TIP3P or other water models without
replacing the Fe NONB lines with the IOD set for that water model.

Other nonbonded terms:

- Cys SG and sulfides: gaff2 S (1.9825, 0.2824).
- Y8/Y9 (O types) are MCPB.py leftovers that no atom in these files uses.
