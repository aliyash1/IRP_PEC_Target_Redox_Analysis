#!/usr/bin/env python3
"""
AlphaFill screen of cofactor binding sites, input = cofactor_residues_FeS_only.json

Input JSON (one entry per protein; one site per cofactor label):
    [{"protein": "FER_ECOLI",
      "cofactors": [{"name": "[2Fe-2S] cluster", "label": 1, "residues": [42, 48, 51, 87]}]}, ...]

Web services used
  UniProt   GET /uniprotkb/{entry_name}.json?fields=accession
            Only to translate "FER_ECOLI" -> "P0A9R4" (AlphaFill is keyed by accession).
            Cached in accession_cache.json; skipped for names in --accessions FILE
            ({"FER_ECOLI": "P0A9R4", ...}) or an "accession" key in the JSON.
  AlphaFill GET /v1/aff/{accession}/json        -> transplants + identity / local RMSD / TCS / residues
  AlphaFill GET /v1/aff/{accession}             -> full AlphaFill structure (mmCIF), only used to compute
                                                   the residues coordinating each transplanted cluster
  AlphaFill GET /v1/aff/{accession}/stripped/{asymlist}
            -> AlphaFill structure with only the listed transplanted compounds (accepted cofactors)

Per site (ACCEPT): correct cofactor + identity + local RMSD + TCS (clash.score) within thresholds
                   AND the AlphaFill residue set covers the table residues.

Why sites were rejected before: "clash.distances" in the JSON only lists atoms near a clash, so
coordinating residues were often "missing_from_alphafill". Two remedies are built in:
  --residue-source union   (default) residues = JSON clash residues  UNION  residues with any atom
                           within 3.0 A of an Fe atom of that transplanted cluster, computed from the
                           AlphaFill structure (matched through the transplant's asym_id).
                           'metadata' = JSON only, 'structure' = structure only.
  --match-mode partial     (default) at least --min-residue-fraction (0.75) of the table residues must be
                           present, i.e. one missing residue is tolerated for a 4-residue site, none for
                           3 residues. Also: 'superset' (all present) and 'exact' (identical sets).

Outputs (in <project>/Results/1_alphafill/):
  all_sites_screening.csv   every site, status, reason, missing/extra residues, diagnostics
  accepted_sites.csv/.xlsx  every ACCEPT site with its structure (protein_all_sites_accepted says whether
                            all screened sites of that protein were accepted)
  accepted_proteins.csv/.xlsx  rows of proteins where every screened site is ACCEPT
  relaxation_sweep.csv      accepted counts under looser thresholds (no new queries)
  structures/               stripped protein + accepted-cofactor structures

Usage
    python alphafill_screen_json.py
    python alphafill_screen_json.py --dry-run
    python alphafill_screen_json.py --residue-source metadata --match-mode superset   # old behaviour
    python alphafill_screen_json.py --min-identity 0.5 --structures fully_accepted

Install: pip install pandas openpyxl requests biopython numpy

Not verified against the live service: JSON field names, that the transplant 'asym_id' equals
label_asym_id in the full mmCIF, the comma-separated asym list for /stripped/, and the compound
codes used with --cofactors all. If a structure cannot be matched the script falls back to the JSON
residues and says so in the 'residue_source_used' column.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests

# --- project layout: <project>/<a>/<b>/this_script.py ---------------------------------
base_dir = Path(__file__).resolve().parent
project_dir = base_dir.parent.parent
INPUT_FILE = "cofactor_residues_FeS_only.json"
input_dir = project_dir / "Inputs" / "json_files"
output_dir = project_dir / "Results" / "1_alphafill"
ALPHAFILL = "https://alphafill.eu/v1/aff"
UNIPROT = "https://rest.uniprot.org/uniprotkb"

MIN_IDENTITY, MAX_LOCAL_RMSD, MAX_TCS = 0.70, 0.92, 0.64
MIN_RESIDUE_FRACTION = 0.75
COORD_CUTOFF = 3.0      # A: protein atom to a cluster Fe atom

# Fe-S clusters (default screen): table name -> AlphaFill compound codes
FES_CODES = {"[4Fe-4S] cluster": {"SF4"}, "[2Fe-2S] cluster": {"FES"}, "[3Fe-4S] cluster": {"F3S"},
             "iron-sulfur cluster": {"SF4", "FES", "F3S"}}
# Other cofactors (--cofactors all), best effort; the hybrid [4Fe-2O-2S] cluster has no code here
OTHER_CODES = {"Fe cation": {"FE", "FE2"}, "Fe(2+)": {"FE2", "FE"}, "Mg(2+)": {"MG"}, "Zn(2+)": {"ZN"},
               "FMN": {"FMN"}, "FAD": {"FAD"}, "NAD(+)": {"NAD"}, "NADP(+)": {"NAP"}, "ATP": {"ATP"},
               "GTP": {"GTP"}, "AMP": {"AMP"}, "S-adenosyl-L-methionine": {"SAM"},
               "pyridoxal 5'-phosphate": {"PLP"}, "thiamine diphosphate": {"TPP"}, "heme b": {"HEM"}}


# ---------------------------------------------------------------------------
# 1. table: JSON -> one row per site
# ---------------------------------------------------------------------------
def read_table(path: Path) -> pd.DataFrame:
    rows = []
    for p in json.load(path.open()):
        for c in p["cofactors"]:
            rows.append({"uniProtkbId": p["protein"], "primaryAccession": p.get("accession"),
                         "cofactor": c["name"], "ligand_label": c["label"],
                         "previous_residues": {int(r) for r in c["residues"]}})
    return pd.DataFrame(rows)


def select_sites(df: pd.DataFrame, which: str) -> pd.DataFrame:
    codes = dict(FES_CODES, **(OTHER_CODES if which == "all" else {}))
    df = df.assign(codes=df["cofactor"].map(lambda n: codes.get(n)))
    return df[df["codes"].notna()].reset_index(drop=True)


def fill_accessions(df: pd.DataFrame, session, out: Path, mapping_file: Path | None) -> pd.DataFrame:
    cache_path = out / "accession_cache.json"
    cache = json.load(cache_path.open()) if cache_path.exists() else {}
    if mapping_file:
        cache.update(json.load(mapping_file.open()))
    todo = [n for n in df.loc[df["primaryAccession"].isna(), "uniProtkbId"].unique() if not cache.get(n)]
    for name in todo:                                   # GET UniProt entry -> accession
        try:
            r = session.get(f"{UNIPROT}/{name}.json", params={"fields": "accession"}, timeout=60)
            cache[name] = r.json().get("primaryAccession") if r.ok else None
        except (requests.RequestException, ValueError):
            cache[name] = None
        time.sleep(0.2)
    cache_path.write_text(json.dumps(cache, indent=1))
    df["primaryAccession"] = df["primaryAccession"].fillna(df["uniProtkbId"].map(cache))
    return df


# ---------------------------------------------------------------------------
# 2. AlphaFill
# ---------------------------------------------------------------------------
def get(session, url, *, as_json, retries=4):
    """None on 404; RuntimeError if the service keeps failing."""
    err = None
    for i in range(retries):
        try:
            r = session.get(url, timeout=180)
            if r.status_code == 404:
                return None
            if r.status_code == 429 or r.status_code >= 500:
                time.sleep(min(float(r.headers.get("Retry-After", 2 ** i)), 30))
                continue
            r.raise_for_status()
            return r.json() if as_json else r.content
        except (requests.RequestException, ValueError) as e:
            err = e
            time.sleep(min(2 ** i, 30))
    raise RuntimeError(f"{url}: {err}")


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def transplants(data: dict) -> list[dict]:
    out = []
    for hit in data.get("hits") or []:
        ident = num((hit.get("alignment") or {}).get("identity"))
        if ident is not None and ident > 1:                 # percent -> fraction
            ident /= 100
        for t in hit.get("transplants") or []:
            clash = t.get("clash") or {}
            res = {int(d["poly_atom"]["seq_id"]) for d in clash.get("distances") or []
                   if (d.get("poly_atom") or {}).get("seq_id") is not None}
            out.append({
                "pdb_id": hit.get("pdb_id"), "asym_id": t.get("asym_id"),
                "codes": {str(t.get("compound_id") or "").upper(), str(t.get("analogue_id") or "").upper()},
                "identity": ident, "local_rmsd": num(t.get("local_rmsd")), "tcs": num(clash.get("score")),
                "metadata_residues": res, "alphafill_residues": set(res), "source_used": "metadata"})
    return out


# ---------------------------------------------------------------------------
# 3. residues coordinating a transplanted cluster, from the AlphaFill mmCIF
# ---------------------------------------------------------------------------
def load_atoms(path: Path) -> dict:
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    d = MMCIF2Dict(str(path))

    def col(k):
        v = d["_atom_site." + k]
        return [v] if isinstance(v, str) else v

    return {"group": np.array(col("group_PDB")), "asym": np.array(col("label_asym_id")),
            "comp": np.char.upper(np.array(col("label_comp_id"))),
            "elem": np.char.upper(np.array(col("type_symbol"))),
            "seq": np.array(col("label_seq_id")),
            "xyz": np.array([col("Cartn_x"), col("Cartn_y"), col("Cartn_z")], dtype=float).T}


def coordination_residues(atoms: dict, asym_id, codes: set) -> set[int] | None:
    """Residue numbers with any atom within COORD_CUTOFF of an Fe atom of ligand `asym_id`."""
    lig = (atoms["asym"] == asym_id) & (atoms["group"] == "HETATM") & np.isin(atoms["comp"], list(codes))
    if not lig.any():
        return None
    fe = lig & (atoms["elem"] == "FE")
    probe = atoms["xyz"][fe if fe.any() else lig]
    prot = (atoms["group"] == "ATOM") & (atoms["seq"] != ".")
    dist = np.linalg.norm(atoms["xyz"][prot][:, None, :] - probe[None, :, :], axis=2).min(axis=1)
    return {int(s) for s in atoms["seq"][prot][dist < COORD_CUTOFF]}


def add_structure_residues(ts, acc, wanted_codes, source, session, out, delay):
    """Replace/extend each relevant transplant's residue set using the full AlphaFill structure."""
    todo = [t for t in ts if t["codes"] & wanted_codes]   # all right-cofactor transplants (keeps the sweep consistent)
    if not todo:
        return
    path = out / "alphafill_full" / f"{acc}.cif"
    try:
        if not path.exists():
            time.sleep(delay)
            content = get(session, f"{ALPHAFILL}/{acc}", as_json=False)           # GET /v1/aff/{id}
            if content is None:
                raise RuntimeError("no AlphaFill structure (404)")
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(content)
        atoms = load_atoms(path)
    except Exception as e:                                # network / parse problems -> keep JSON residues
        for t in todo:
            t["source_used"] = f"metadata (structure unavailable: {str(e)[:80]})"
        return
    for t in todo:
        coord = coordination_residues(atoms, t["asym_id"], t["codes"])
        if coord is None:
            t["source_used"] = "metadata (ligand asym_id not found in structure)"
            continue
        t["structure_residues"] = coord
        t["alphafill_residues"] = coord if source == "structure" else (t["metadata_residues"] | coord)
        t["source_used"] = source


# ---------------------------------------------------------------------------
# 4. site screening
# ---------------------------------------------------------------------------
def residues_match(prev: set, a: set, mode: str, frac: float) -> bool:
    if mode == "exact":
        return a == prev
    if mode == "superset":
        return prev <= a
    return len(prev & a) / len(prev) >= frac - 1e-9


def failures(t) -> list[str]:
    f = []
    if t["identity"] is None or t["identity"] < MIN_IDENTITY:
        f.append(f"identity<{MIN_IDENTITY:.0%}")
    if t["local_rmsd"] is None or t["local_rmsd"] >= MAX_LOCAL_RMSD:
        f.append(f"local_RMSD>={MAX_LOCAL_RMSD}")
    if t["tcs"] is None or t["tcs"] >= MAX_TCS:
        f.append(f"TCS>={MAX_TCS}")
    return f


def screen_site(site, ts, mode, frac, used) -> dict:
    prev = site["previous_residues"]
    right = [t for t in ts if t["codes"] & site["codes"]]
    if not right:
        return {"status": "REJECT", "reason": "no transplant with the expected cofactor"}
    ok = lambda a: residues_match(prev, a, mode, frac)
    hits = sorted((t for t in right if not failures(t) and ok(t["alphafill_residues"])),
                  key=lambda t: ((t["pdb_id"], t["asym_id"]) in used, t["local_rmsd"], t["tcs"]))
    if hits:
        t, status, reason = hits[0], "ACCEPT", "all criteria pass"
        used.add((t["pdb_id"], t["asym_id"]))
    else:
        t = sorted(right, key=lambda t: (len(failures(t)), len(prev - t["alphafill_residues"])))[0]
        why = failures(t)
        if not ok(t["alphafill_residues"]):
            why.append(f"residues {len(prev & t['alphafill_residues'])}/{len(prev)} present")
        status, reason = "REVIEW", "; ".join(why)
    a = t["alphafill_residues"]
    single = {"identity": [], "rmsd": [], "tcs": [], "residues": []}
    for c in right:                                     # transplants failing exactly ONE criterion
        bad = [k for k, v in (("identity", c["identity"] is not None and c["identity"] >= MIN_IDENTITY),
                              ("rmsd", c["local_rmsd"] is not None and c["local_rmsd"] < MAX_LOCAL_RMSD),
                              ("tcs", c["tcs"] is not None and c["tcs"] < MAX_TCS),
                              ("residues", ok(c["alphafill_residues"]))) if not v]
        if len(bad) == 1:
            single[bad[0]].append(c)
    return {"status": status, "reason": reason, "pdb_id": t["pdb_id"], "asym_id": t["asym_id"],
            "identity": t["identity"], "local_rmsd": t["local_rmsd"], "tcs": t["tcs"],
            "alphafill_residues": sorted(a), "residue_source_used": t["source_used"],
            "missing_from_alphafill": sorted(prev - a),     # previous_residues - alphafill_residues
            "extra_in_alphafill": sorted(a - prev),         # alphafill_residues - previous_residues
            "n_transplants_correct_cofactor": len(right),
            "best_identity_any": max((c["identity"] for c in right if c["identity"] is not None), default=None),
            "accepted_if_only_relaxed": ",".join(k for k, v in single.items() if v) if status != "ACCEPT" else "",
            "identity_needed": max((c["identity"] for c in single["identity"] if c["identity"] is not None),
                                   default=None),
            "local_rmsd_needed": min((c["local_rmsd"] for c in single["rmsd"] if c["local_rmsd"] is not None),
                                     default=None),
            "tcs_needed": min((c["tcs"] for c in single["tcs"] if c["tcs"] is not None), default=None)}


def relaxation_sweep(site_cands, frac) -> pd.DataFrame:
    """Re-apply the filters with looser settings (uses the residue sets already computed)."""
    rows = []
    for mode in ("exact", "superset", "partial"):
        for ident in (0.70, 0.60, 0.50, 0.40, 0.30, 0.0):
            for rmsd, tcs in ((0.92, 0.64), (1.5, 1.0), (2.0, 1.5)):
                n_sites, per_prot = 0, {}
                for name, prev, right in site_cands:
                    hit = any(c["identity"] is not None and c["local_rmsd"] is not None and c["tcs"] is not None
                              and c["identity"] >= ident and c["local_rmsd"] < rmsd and c["tcs"] < tcs
                              and residues_match(prev, c["alphafill_residues"], mode, frac) for c in right)
                    n_sites += hit
                    per_prot[name] = per_prot.get(name, True) and hit
                rows.append({"residue_match": mode, "min_identity": ident, "max_local_rmsd": rmsd,
                             "max_tcs": tcs, "sites_accepted": n_sites, "of_sites": len(site_cands),
                             "proteins_accepted": sum(per_prot.values()), "of_proteins": len(per_prot)})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 5. main
# ---------------------------------------------------------------------------
def main():
    global ALPHAFILL, MIN_IDENTITY, MAX_LOCAL_RMSD, MAX_TCS
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("table", nargs="?", type=Path, default=input_dir / INPUT_FILE,
                    help="JSON file, or a folder containing %s (default: %%(default)s)" % INPUT_FILE)
    ap.add_argument("--out", type=Path, default=output_dir, help="output folder (default: %(default)s)")
    ap.add_argument("--match-mode", choices=["exact", "superset", "partial"], default="partial")
    ap.add_argument("--min-residue-fraction", type=float, default=MIN_RESIDUE_FRACTION,
                    help="partial mode: fraction of table residues that must be present (default %(default)s)")
    ap.add_argument("--residue-source", choices=["metadata", "structure", "union"], default="union")
    ap.add_argument("--structures", choices=["any_accepted", "fully_accepted"], default="any_accepted",
                    help="download stripped structures for proteins with >=1 accepted site (default) "
                         "or only where all screened sites are accepted")
    ap.add_argument("--cofactors", choices=["fes", "all"], default="fes",
                    help="fes = Fe-S clusters only (default); all = also metals/FMN/FAD/NAD/SAM... (best effort)")
    ap.add_argument("--accessions", type=Path, help='JSON {"FER_ECOLI": "P0A9R4", ...} to skip UniProt lookups')
    ap.add_argument("--min-identity", type=float, default=MIN_IDENTITY, help="fraction, default %(default)s")
    ap.add_argument("--max-local-rmsd", type=float, default=MAX_LOCAL_RMSD, help="A, default %(default)s")
    ap.add_argument("--max-tcs", type=float, default=MAX_TCS, help="default %(default)s")
    ap.add_argument("--base-url", default=ALPHAFILL)
    ap.add_argument("--delay", type=float, default=0.5)
    ap.add_argument("--dry-run", action="store_true", help="only read/group the table and exit")
    args = ap.parse_args()
    ALPHAFILL = args.base_url
    MIN_IDENTITY = args.min_identity
    MAX_LOCAL_RMSD, MAX_TCS = args.max_local_rmsd, args.max_tcs
    print(f"Thresholds: identity >= {MIN_IDENTITY:.0%}, local RMSD < {MAX_LOCAL_RMSD} A, TCS < {MAX_TCS}; "
          f"residues: {args.match_mode}"
          f"{f' (>= {args.min_residue_fraction:.0%} of table residues)' if args.match_mode == 'partial' else ''}"
          f", source = {args.residue_source}")
    if args.table.is_dir():
        args.table = args.table / INPUT_FILE
    if not args.table.exists():
        raise SystemExit(f"Input file not found: {args.table}")
    for d in ("json", "structures"):
        (args.out / d).mkdir(parents=True, exist_ok=True)
    session = requests.Session()
    session.headers["User-Agent"] = "alphafill-screen-json/2.0"

    df = read_table(args.table)
    sites = select_sites(df, args.cofactors)
    print(f"{len(df)} cofactor sites in file; {len(sites)} to screen "
          f"({sites['uniProtkbId'].nunique()} proteins; --cofactors {args.cofactors})")
    if args.dry_run:
        print(sites["cofactor"].value_counts().to_string())
        return

    sites = fill_accessions(sites, session, args.out, args.accessions)
    results, site_cands, missing = [], [], {"identity": 0, "local_rmsd": 0, "tcs": 0, "total": 0}
    for name, group in sites.groupby("uniProtkbId", sort=False):
        acc = group["primaryAccession"].iloc[0]
        base = [s.drop("codes").to_dict() for _, s in group.iterrows()]
        if not isinstance(acc, str):
            results += [{**b, "status": "REVIEW", "reason": "UniProt accession not found"} for b in base]
            print(f"{name}: no accession")
            continue
        jpath = args.out / "json" / f"{acc}.json"
        try:
            if jpath.exists():
                data = json.load(jpath.open())
            else:
                time.sleep(args.delay)
                data = get(session, f"{ALPHAFILL}/{acc}/json", as_json=True)       # GET /v1/aff/{id}/json
                if data is not None:
                    jpath.write_text(json.dumps(data))
        except RuntimeError as e:
            results += [{**b, "status": "REVIEW", "reason": f"AlphaFill query failed (not a reject): {e}"}
                        for b in base]
            print(f"{name} {acc}: query failed")
            continue
        if data is None:
            results += [{**b, "status": "REJECT", "reason": "no AlphaFill entry (404)"} for b in base]
            print(f"{name} {acc}: no AlphaFill entry")
            continue
        ts, used, n0 = transplants(data), set(), len(results)
        missing["total"] += len(ts)
        for t in ts:
            for k in ("identity", "local_rmsd", "tcs"):
                missing[k] += t[k] is None
        if args.residue_source != "metadata":
            wanted = set().union(*group["codes"])
            add_structure_residues(ts, acc, wanted, args.residue_source, session, args.out, args.delay)
        for (_, s), b in zip(group.iterrows(), base):
            results.append({**b, **screen_site(s, ts, args.match_mode, args.min_residue_fraction, used)})
            site_cands.append((name, s["previous_residues"], [t for t in ts if t["codes"] & s["codes"]]))
        print(f"{name} {acc}: " + ", ".join(r["status"] for r in results[n0:]))

    if missing["total"] and any(missing[k] for k in ("identity", "local_rmsd", "tcs")):
        print(f"\nWARNING: of {missing['total']} transplants, values missing for identity={missing['identity']}, "
              f"local_rmsd={missing['local_rmsd']}, tcs={missing['tcs']}. Missing values always fail the filters; "
              "if these are large, the JSON field names differ from what the script expects - inspect json/.")
    sweep = relaxation_sweep(site_cands, args.min_residue_fraction)
    sweep.to_csv(args.out / "relaxation_sweep.csv", index=False)
    print("\nProteins accepted if identity is relaxed (RMSD/TCS as set, residue_match as column):")
    sub = sweep[(sweep["max_local_rmsd"] == 0.92) & (sweep["max_tcs"] == 0.64)]
    print(sub.pivot(index="min_identity", columns="residue_match", values="proteins_accepted")
             .sort_index(ascending=False).to_string())

    res = pd.DataFrame(results)
    res["previous_residues"] = res["previous_residues"].map(sorted)
    res.to_csv(args.out / "all_sites_screening.csv", index=False)
    acc_sites = res[res["status"] == "ACCEPT"]
    per_prot = res.groupby("uniProtkbId")["status"].agg(
        n_sites="size", n_accepted=lambda s: (s == "ACCEPT").sum())
    per_prot["all"] = per_prot["n_accepted"] == per_prot["n_sites"]
    keep = per_prot[per_prot["all"] if args.structures == "fully_accepted" else per_prot["n_accepted"] > 0].index

    rows = []
    for name in keep:
        sub = acc_sites[acc_sites["uniProtkbId"] == name]
        acc, asyms = sub["primaryAccession"].iloc[0], sorted(set(sub["asym_id"]))
        status, path = "no_structure", ""
        try:                                              # GET /v1/aff/{id}/stripped/{asymlist}
            time.sleep(args.delay)
            content = get(session, f"{ALPHAFILL}/{acc}/stripped/{','.join(asyms)}", as_json=False)
            if content:
                p = args.out / "structures" / f"{acc}_{'_'.join(asyms)}_stripped.cif"
                p.write_bytes(content)
                path = str(p)
                status = "downloaded" if content[:100].lstrip().startswith(b"data_") else "check_format"
        except RuntimeError as e:
            status = f"failed: {e}"
        for _, r in sub.iterrows():
            rows.append({"uniProtkbId": name, "primaryAccession": acc, "cofactor": r["cofactor"],
                         "ligand_label": r["ligand_label"], "table_residues": r["previous_residues"],
                         "alphafill_residues": r["alphafill_residues"],
                         "missing_from_alphafill": r["missing_from_alphafill"],
                         "residue_source_used": r["residue_source_used"], "pdb_id": r["pdb_id"],
                         "identity": r["identity"], "local_rmsd": r["local_rmsd"], "tcs": r["tcs"],
                         "transplant_asym_id": r["asym_id"],
                         "sites_accepted_in_protein": f"{per_prot.loc[name, 'n_accepted']}/{per_prot.loc[name, 'n_sites']}",
                         "protein_all_sites_accepted": bool(per_prot.loc[name, "all"]),
                         "alphafill_structure": path, "download_status": status})
    accepted = pd.DataFrame(rows)
    accepted.to_csv(args.out / "accepted_sites.csv", index=False)
    accepted.to_excel(args.out / "accepted_sites.xlsx", index=False)
    full = accepted[accepted["protein_all_sites_accepted"]] if len(accepted) else accepted
    full.to_csv(args.out / "accepted_proteins.csv", index=False)
    full.to_excel(args.out / "accepted_proteins.xlsx", index=False)

    print("\nSite status counts:\n" + res["status"].value_counts().to_string())
    print(f"\nProteins with >=1 accepted site: {int((per_prot['n_accepted'] > 0).sum())}; "
          f"with ALL screened sites accepted: {int(per_prot['all'].sum())} of {len(per_prot)}")
    if len(accepted):
        with pd.option_context("display.width", 220, "display.max_columns", None):
            print(accepted[["uniProtkbId", "primaryAccession", "cofactor", "ligand_label", "pdb_id",
                            "identity", "sites_accepted_in_protein", "residue_source_used"]].to_string(index=False))
    print(f"\nWrote tables and structures to {args.out}/")


if __name__ == "__main__":
    main()