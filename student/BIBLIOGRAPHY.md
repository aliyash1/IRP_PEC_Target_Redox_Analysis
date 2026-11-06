# Fe-S cluster force-field parameters: literature found

Every source with Fe-S cluster parameters that turned up in the search (October 2026), grouped by what is
usable. All DOIs and bibliographic data were checked against Crossref. BibTeX entries are in
`FeS_parameters.bib` (keys in square brackets).

Abbreviations: ox/red = oxidised/reduced; Sem. = Seminario method; MCPB = MCPB.py.

## 1. Used to build the folders in this directory

| Ref | Cluster / ligands / states | Method | Folder(s) |
|---|---|---|---|
| [Carvalho2013] Carvalho, Teixeira, Ramos, *J. Comput. Chem.* **2013**, 34, 1540–1548. doi:10.1002/jcc.23287 | Rubredoxin Fe(III)(Cys)4; desulforedoxin Fe(Cys)4; plant Fd [2Fe-2S]2+(Cys)4 | B3LYP, Sem., RESP; ff99SB | `1Fe_threeplus_4Cys_Carvalho2013`, `2Fe2S_twoplus_4Cys_Carvalho2013` (desulforedoxin charges not integer, not built) |
| [CarvalhoSwart2014] Carvalho, Swart, *J. Chem. Inf. Model.* **2014**, 54, 613–620. doi:10.1021/ci400718m | [2Fe-2S]2+/+, [3Fe-4S]+/0, [4Fe-4S]2+/+, all (Cys) | OPBE/TZP, Sem., MDC-d | `4Fe4S_twoplus_4Cys_CarvalhoSwart2014`, `4Fe4S_oneplus_4Cys_CarvalhoSwart2014` |
| [CarvalhoSwart2015] Carvalho, Swart, *J. Chem. Inf. Model.* **2015**, 55, 1508 (Correction). doi:10.1021/acs.jcim.5b00390 | Corrected SI with all angle parameters and charges | – | (same) |
| [Nakayoshi2022] Nakayoshi, Ohnishi, Tanaka, Kurisu, Kondo, Takano, *Int. J. Mol. Sci.* **2022**, 23, 15913. doi:10.3390/ijms232415913 | Plant Fd [2Fe-2S]2+/+ (Cys)4 | MCPB, Sem. (LC-wPBE), RESP; ff14SB. **Note:** Fe LJ is the IOD Fe3+ set for OPC, but the MD used TIP3P water | `2Fe2S_twoplus_4Cys_Nakayoshi2022`, `2Fe2S_oneplus_4Cys_Nakayoshi2022` |
| [Teo2020] Teo, Migliore, Beratan, *Chem. Sci.* **2020**, 11, 7076–7085. doi:10.1039/d0sc02245d | p58c (HiPIP-like) [4Fe-4S]3+/2+ (Cys)4 | B3LYP, Sem., RESP; ff14SB | `4Fe4S_threeplus_4Cys_Teo2020`, `4Fe4S_twoplus_4Cys_Teo2020` |
| [Dawson2016] Dawson, Jono, Terada, Shimizu, *PLoS One* **2016**, 11, e0162031. doi:10.1371/journal.pone.0162031 | Rieske [2Fe-2S]2+/+ (Cys)2(His)2 (+ non-heme Fe site) | Hand-fitted to QM/PDB, incl. torsions; files in S4 File | `2Fe2S_twoplus_2Cys2His_Dawson2016`, `2Fe2S_oneplus_2Cys2His_Dawson2016` |
| [Zuo2021] Zuo, Marjault, Bren, Rossetti, Nechushtai, Carloni, *J. Biol. Inorg. Chem.* **2021**, 26, 763–774. doi:10.1007/s00775-021-01890-8 | mitoNEET/NAF-1 [2Fe-2S]+ (Cys)3(His), ox charges also given | MCPB, Sem., RESP; AMBER99SB-ILDN (GROMACS) | `2Fe2S_oneplus_3Cys1His_Zuo2021` (ox not built) |

## 2. Other AMBER-compatible sets with open numbers (not built)

| Ref | Cluster / ligands / states | Content / problem |
|---|---|---|
| [Pabis2014] Pabis, Geronimo, York, Paneth, *J. Chem. Theory Comput.* **2014**, 10, 2246–2254. doi:10.1021/ct500205z | Rieske [2Fe-2S]2+ (Cys)2(His)2, ox (nitrobenzene dioxygenase) | MCPB; full charges, bonds, angles in SI (incl. Fe-Fe bond). His ND1 types N0/N1 used inconsistently. `_sources/amber_nbdo` |
| [Pesce2017] Pesce, Calandrini, Marjault, Lipper, Rossetti, Mittler, Jennings, Bauer et al., *J. Phys. Chem. B* **2017**, 121, 10648–10656. doi:10.1021/acs.jpcb.7b10584 | mitoNEET/NAF-1 [2Fe-2S] (Cys)3(His), His protonated/deprotonated | MCPB; bonds/angles (GROMACS units) + charges in SI Tables S1–S3; charges do not sum to integers. Predecessor of [Zuo2021]. `_sources/jpcb2017_NEET` |
| [Blachly2013] Blachly, de Oliveira, Williams, McCammon, *PLoS Comput. Biol.* **2013**, 9, e1003395. doi:10.1371/journal.pcbi.1003395 | IspH [4Fe-4S]2+ (Cys)3 + substrate/water on unique Fe | RESP charges for cluster + 3 Cys, Fe/S LJ; **no bonded terms** (geometry restrained). `_sources/amber_isph` |
| [Tendwa2021] Tendwa, Chebon-Bore, Lobb, Musyoka, Tastan Bishop, *Molecules* **2021**, 26, 2929. doi:10.3390/molecules26102929 | DPD [4Fe-4S]2+ (Cys)4 and (Cys)3(Gln) | MCPB bonds/angles + VFFDT set; charges only for Fe/S/ligating atoms. `_sources/molecules2021` |
| [ChebonBore2021] Chebon-Bore, Sanyanga, Manyumwa, Khairallah, Tastan Bishop, *Int. J. Mol. Sci.* **2021**, 22, 2138. doi:10.3390/ijms22042138 | Rieske [2Fe-2S] (Cys)2(His)2 in P. falciparum cytb-ISP | MCPB frcmod/mol2 printed in SI, but garbled (atom names mixed with types, Fe labels swapped between tables). |
| [Kaukonen2008] Kaukonen, Söderhjelm, Heimdal, Ryde, *J. Chem. Theory Comput.* **2008**, 4, 985–1001. doi:10.1021/ct700347h | [NiFe]-hydrogenase [4Fe-4S]2+ (Cys)4 | QM bonded terms + ESP charges; files posted on the AMBER mailing list (archive.ambermd.org/200906/att-0647/). `_sources/amberlist/kaukonen` |
| [CaseFe4S4] D. A. Case lab, `frcmod.fe4s4`, `FS4.off`, `CYF.off`, AMBER mailing list, June 2009 (archive.ambermd.org/200906/att-0507/fe4s4.tar.gz); charges from [Mouesca1994] | [4Fe-4S]2+ (Cys)4 (HiPIP 1IUA) | Simple complete set; README says force constants "not fully optimized". `_sources/amberlist/case_fe4s4` |
| [BryceDB] R. Bryce, AMBER Parameter Database, "Fe4S4" entry (Szaleniec/Swart), University of Manchester (amber.manchester.ac.uk, via web.archive.org) | [4Fe-4S]2+ (Cys)4 (NapA) | Ready-made implementation of [CarvalhoSwart2014]; Cys C/O charges look swapped. `_sources/bryce` |
| [Yadav2026] Yadav, Pantazis, *Chem. Sci.* **2026**, doi:10.1039/d6sc04231g; data: Edmond doi:10.17617/3.AUOKHW | Ni-CODH B- and D-cluster [4Fe-4S] (Cys)4; C-cluster [NiFe4S4] | MCPB frcmod/mol2 deposited; hand-edited, cluster charges non-integer, redox states not stated. |
| [Tang2026] Tang et al., *Nat. Commun.* **2026**, 17, 8422. doi:10.1038/s41467-026-73842-2; data: Zenodo doi:10.5281/zenodo.18410016 | NreB [2Fe-2S] (Cys)3 + terminal SH | GROMACS itp from Sobtop; zero charges, bonded terms "guessed". Not recommended. |

## 3. AMBER-compatible sets that exist but could not be retrieved (paywalled SI)

| Ref | Cluster / ligands | Where the numbers are |
|---|---|---|
| [Smith2012] Smith, Xiong, Straatsma, Rosso, Squier, *J. Chem. Theory Comput.* **2012**, 8, 2103–2114. doi:10.1021/ct300185u | [NiFe]-hydrogenase proximal [4Fe-4S] (Cys)4, medial [3Fe-4S], distal [4Fe-4S] (Cys)3(His), Ni-Fe site | SI ct300185u_si_001. Main AMBER source for [3Fe-4S] and Cys3His [4Fe-4S]; reused for Photosystem I and MUTYH. |
| [Smith2014] Smith, Danyal, Raugei, Seefeldt, *Biochemistry* **2014**, 53, 2278–2285. doi:10.1021/bi401313j | Nitrogenase FeMo-co, P-cluster | SI. |

## 4. CHARMM / GROMOS sets (not directly usable with AMBER)

| Ref | Cluster / ligands / states | Availability |
|---|---|---|
| [ChangKim2009] Chang, Kim, *J. Chem. Theory Comput.* **2009**, 5, 1137–1145. doi:10.1021/ct800342w | [4Fe-4S] (Cys)4 and (Cys)3(His), [2Fe-2S] (Cys)4, H-cluster; +/2+ | Paywalled. Most-used CHARMM source. |
| [KailaLab] Kaila lab, `ff_parameters`, github.com/KailaLab/ff_parameters (MIT); refs [Muehlbauer2021], [Sirohiwal2024] | Rieske [2Fe-2S] (Cys)2(His)2 in 8 redox/protonation states; [4Fe-4S] (Cys)4 ox/red | Open, complete (CHARMM36). `_sources/charmm_kailalab_github` |
| [Muehlbauer2021] Mühlbauer, Gamiz-Hernandez, Kaila, *J. Am. Chem. Soc.* **2021**, 143, 20873–20883. doi:10.1021/jacs.1c09356 | [4Fe-4S] (Cys)4 | Reference for the KailaLab [4Fe-4S] file. |
| [Sirohiwal2024] Sirohiwal, Gamiz-Hernandez, Kaila, *J. Am. Chem. Soc.* **2024**, 146, 18019–18031. doi:10.1021/jacs.4c04476 | [4Fe-4S] (Cys)4 | Reference for the KailaLab [4Fe-4S] file. |
| [Sharma2015] Sharma et al., *Proc. Natl. Acad. Sci. USA* **2015**, 112, 11571–11576. doi:10.1073/pnas.1503761112 | Complex I Fe-S clusters, ox (own charges) / red (from [ChangKim2009]) | No numbers in SI. |
| [Kaszuba2013] Kaszuba, Postila, Cramariuc, Sarewicz, Osyczka, Vattulainen, Róg, *Theor. Chem. Acc.* **2013**, 132, 1370. doi:10.1007/s00214-013-1370-8 | Rieske [2Fe-2S] (Cys)2(His)2 ox/red (bc1 complex) | Open (green OA); charges in Fig. 4, bonded terms in ESM; no Fe/S LJ. |
| [Postila2013] Postila, Kaszuba, Sarewicz, Osyczka, Vattulainen, Róg, *Biochim. Biophys. Acta Bioenerg.* **2013**, 1827, 761–768. doi:10.1016/j.bbabio.2013.02.005 | Rieske (uses [Kaszuba2013]) | – |
| [Barragan2015] Barragan, Crofts, Schulten, Solov'yov, *J. Phys. Chem. B* **2015**, 119, 433–447. doi:10.1021/jp510022w | Rieske with deprotonated His, ox | Charges/topology in SI (Kaszuba types). |
| [Camilo2025] Camilo, Arantes, *J. Chem. Inf. Model.* **2025**, 65, 6184–6197. doi:10.1021/acs.jcim.5c00655 | Rieske, ox, His deprotonated | Charges in SI Table S2 and Zenodo (doi:10.5281/zenodo.15084037); bonded terms from [ChangKim2009]. |
| [Botticelli2024] Botticelli, La Penna, Minicozzi, Stellato, Morante, Rossi, Faraloni, *Int. J. Mol. Sci.* **2024**, 25, 3663. doi:10.3390/ijms25073663 | [4Fe-4S]2+ (Cys)4 + [2Fe]H | Open CHARMM36 files in SI (adapted [ChangKim2009], averaged charges, modified LJ). |
| [Seki2026] Seki, Ohnuki, Okazaki, *J. Chem. Inf. Model.* **2026**, 66, 4878–4887. doi:10.1021/acs.jcim.6c00347 | Na+-NQR [2Fe-2S] (Cys)4 ox/red | Charges in SI Fig. S4; bonded "adapted from Kaila group". |
| [Teixeira2006] Teixeira, Baptista, Soares, *Biophys. J.* **2006**, 91, 2035–2045. doi:10.1529/biophysj.106.084376 | [NiFe]-hydrogenase [4Fe-4S]2+ ×2, [3Fe-4S]+ | GROMOS 43A1; charges in online supplement (not retrievable). |
| [Yelle1995] Yelle, Park, Ichiye, *Proteins* **1995**, 22, 154–167. doi:10.1002/prot.340220208 | Rubredoxin Fe(Cys)4 | CHARMM; cited in [Carvalho2013]. |

## 5. Charge sets for electrostatics, and nonbonded (LJ) sources used above

| Ref | Content |
|---|---|
| [Mouesca1994] Mouesca, Chen, Noodleman, Bashford, Case, *J. Am. Chem. Soc.* **1994**, 116, 11898–11914. doi:10.1021/ja00105a033 | DFT charges of Fe-S clusters (origin of [CaseFe4S4] charges). |
| [Ullmann2002] Ullmann, Noodleman, Case, *J. Biol. Inorg. Chem.* **2002**, 7, 632–639. doi:10.1007/s00775-002-0342-6 | Rieske charges for continuum electrostatics (supplement). |
| [Torres2003] Torres, Lovell, Noodleman, Case, *J. Am. Chem. Soc.* **2003**, 125, 1923–1936. doi:10.1021/ja0211104 | [4Fe-4S] charges/redox calculations. |
| [Giammona1984] D. A. Giammona, PhD thesis, University of California, Davis, **1984** | Fe LJ (R* 1.20, ε 0.05) used in [Carvalho2013], [Dawson2016], [Blachly2013], [Kaukonen2008]. |
| [Rappe1992] Rappé, Casewit, Colwell, Goddard, Skiff, *J. Am. Chem. Soc.* **1992**, 114, 10024–10035. doi:10.1021/ja00051a040 | UFF Fe LJ, used for the Carvalho & Swart folders (via [BryceDB]). |
| [Li2013] Li, Roberts, Chakravorty, Merz, *J. Chem. Theory Comput.* **2013**, 9, 2733–2748. doi:10.1021/ct400146w | IOD Fe2+ LJ (used in [Tendwa2021]). |
| [Li2021] Li, Song, Li, Merz, *J. Chem. Theory Comput.* **2021**, 17, 2342–2354. doi:10.1021/acs.jctc.0c01320 | IOD Fe3+ LJ for OPC water (used by MCPB.py in your own folders and [Nakayoshi2022]). |

No open parameters were found for: radical-SAM clusters with SAM on the unique Fe, aconitase,
Asp/Ser-ligated [4Fe-4S], the [4Fe-3S] proximal cluster, IscU/glutaredoxin [2Fe-2S].
