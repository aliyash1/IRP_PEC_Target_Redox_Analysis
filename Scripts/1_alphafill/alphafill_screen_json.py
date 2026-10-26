#!/usr/bin/env python3
"""
AlphaFill screen of cofactor binding sites, input = cofactor_residues_FeS_only.json

Input JSON (one entry per protein; one site per cofactor label):
    [{"protein": "FER_ECOLI",
      "cofactors": [{"name": "[2Fe-2S] cluster", "label": 1, "residues": [42, 48, 51, 87]}]}, ...]

Which web services are used
  UniProt   GET /uniprotkb/{entry_name}.json?fields=accession
            Only to translate "FER_ECOLI" -> "P0A9R4" (AlphaFill is keyed by accession).
            Cached in accession_cache.json; skipped for names listed in --accessions FILE
            ({"FER_ECOLI": "P0A9R4", ...}) or given as an "accession" key in the JSON.
  AlphaFill GET /v1/aff/{accession}/json       -> all transplants + identity / local RMSD / TCS / residues
  AlphaFill GET /v1/aff/{accession}/stripped/{asymlist}
            -> AlphaFill structure with only the listed transplanted compounds (accepted cofactors)

Per site: correct cofactor + identity >= 70 % + local RMSD < 0.92 A + TCS (clash.score) < 0.64
          + residue set  {d["poly_atom"]["seq_id"] for d in transplant["clash"]["distances"]}
            equal to the table residues (--match-mode exact) or containing them (superset).
A protein is ACCEPTED when every screened site is ACCEPT.

Usage
    python alphafill_screen_json.py                      # reads <project>/Inputs/json_files/cofactor_residues_FeS_only.json
    python alphafill_screen_json.py --dry-run
    python alphafill_screen_json.py --match-mode superset
Outputs go to <project>/Results/1_alphafill/ (override with --out; give another file as the first argument).

Install: pip install pandas openpyxl requests

Not verified against the live service: JSON field names, the comma-separated asym list for
/stripped/, and the compound codes used with --cofactors all.
"""
from __future__ import annotations

import argparse
import json
import re
import time
from pathlib import Path

import pandas as pd
import requests

# --- settings---------------------------------
base_dir = Path(__file__).resolve().parent
project_dir = base_dir.parent.parent
INPUT_FILE = "cofactor_residues_FeS_only.json"
input_dir = project_dir / "Inputs" / "json_files"
output_dir = project_dir / "Results" / "1_alphafill"
ALPHAFILL = "https://alphafill.eu/v1/aff"
UNIPROT = "https://rest.uniprot.org/uniprotkb"

MIN_IDENTITY, MAX_LOCAL_RMSD, MAX_TCS = 0.50, 0.92, 0.64

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
            r = session.get(url, timeout=120)
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
            out.append({
                "pdb_id": hit.get("pdb_id"), "asym_id": t.get("asym_id"),
                "codes": {str(t.get("compound_id") or "").upper(), str(t.get("analogue_id") or "").upper()},
                "identity": ident, "local_rmsd": num(t.get("local_rmsd")), "tcs": num(clash.get("score")),
                "alphafill_residues": {int(d["poly_atom"]["seq_id"]) for d in clash.get("distances") or []
                                       if (d.get("poly_atom") or {}).get("seq_id") is not None}})
    return out


def failures(t) -> list[str]:
    f = []
    if t["identity"] is None or t["identity"] < MIN_IDENTITY:
        f.append("identity<70%")
    if t["local_rmsd"] is None or t["local_rmsd"] >= MAX_LOCAL_RMSD:
        f.append("local_RMSD>=0.92")
    if t["tcs"] is None or t["tcs"] >= MAX_TCS:
        f.append("TCS>=0.64")
    return f


def screen_site(site, ts, mode, used) -> dict:
    prev = site["previous_residues"]
    right = [t for t in ts if t["codes"] & site["codes"]]
    if not right:
        return {"status": "REJECT", "reason": "no transplant with the expected cofactor"}
    ok = (lambda a: a == prev) if mode == "exact" else (lambda a: prev <= a)
    hits = sorted((t for t in right if not failures(t) and ok(t["alphafill_residues"])),
                  key=lambda t: ((t["pdb_id"], t["asym_id"]) in used, t["local_rmsd"], t["tcs"]))
    if hits:
        t, status, reason = hits[0], "ACCEPT", "all criteria pass"
        used.add((t["pdb_id"], t["asym_id"]))
    else:
        t = sorted(right, key=lambda t: (len(failures(t)), len(prev ^ t["alphafill_residues"])))[0]
        why = failures(t) + ([] if ok(t["alphafill_residues"]) else ["residue_set_mismatch"])
        status, reason = "REVIEW", "; ".join(why)
    a = t["alphafill_residues"]
    return {"status": status, "reason": reason, "pdb_id": t["pdb_id"], "asym_id": t["asym_id"],
            "identity": t["identity"], "local_rmsd": t["local_rmsd"], "tcs": t["tcs"],
            "alphafill_residues": sorted(a),
            "missing_from_alphafill": sorted(prev - a),     # previous_residues - alphafill_residues
            "extra_in_alphafill": sorted(a - prev),         # alphafill_residues - previous_residues
            "n_transplants_correct_cofactor": len(right)}


# ---------------------------------------------------------------------------
# 3. main
# ---------------------------------------------------------------------------
def main():
    global ALPHAFILL
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("table", nargs="?", type=Path, default=input_dir / INPUT_FILE,
                    help="JSON file, or a folder containing %s (default: %%(default)s)" % INPUT_FILE)
    ap.add_argument("--out", type=Path, default=output_dir, help="output folder (default: %(default)s)")
    ap.add_argument("--match-mode", choices=["exact", "superset"], default="exact")
    ap.add_argument("--cofactors", choices=["fes", "all"], default="fes",
                    help="fes = Fe-S clusters only (default); all = also metals/FMN/FAD/NAD/SAM... (best effort)")
    ap.add_argument("--accessions", type=Path, help='JSON {"FER_ECOLI": "P0A9R4", ...} to skip UniProt lookups')
    ap.add_argument("--base-url", default=ALPHAFILL)
    ap.add_argument("--delay", type=float, default=0.5)
    ap.add_argument("--dry-run", action="store_true", help="only read/group the table and exit")
    args = ap.parse_args()
    ALPHAFILL = args.base_url
    if args.table.is_dir():
        args.table = args.table / INPUT_FILE
    if not args.table.exists():
        raise SystemExit(f"Input file not found: {args.table}")
    for d in ("json", "structures"):
        (args.out / d).mkdir(parents=True, exist_ok=True)
    session = requests.Session()
    session.headers["User-Agent"] = "alphafill-screen-json/1.0"

    df = read_table(args.table)
    print(df.columns.tolist())
    print(df.head())
    sites = select_sites(df, args.cofactors)
    print(f"{len(df)} cofactor sites in file; {len(sites)} to screen "
          f"({sites['uniProtkbId'].nunique()} proteins; --cofactors {args.cofactors})")
    if args.dry_run:
        print(sites["cofactor"].value_counts().to_string())
        return

    sites = fill_accessions(sites, session, args.out, args.accessions)
    results = []
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
        for (_, s), b in zip(group.iterrows(), base):
            results.append({**b, **screen_site(s, ts, args.match_mode, used)})
        print(f"{name} {acc}: " + ", ".join(r["status"] for r in results[n0:]))

    res = pd.DataFrame(results)
    res["previous_residues"] = res["previous_residues"].map(sorted)
    all_ok = res.groupby("uniProtkbId")["status"].agg(lambda s: set(s) == {"ACCEPT"})

    rows = []
    for name in all_ok[all_ok].index:
        sub = res[res["uniProtkbId"] == name]
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
                         "alphafill_residues": r["alphafill_residues"], "pdb_id": r["pdb_id"],
                         "identity": r["identity"], "local_rmsd": r["local_rmsd"], "tcs": r["tcs"],
                         "transplant_asym_id": r["asym_id"], "alphafill_structure": path,
                         "download_status": status})
    accepted = pd.DataFrame(rows)

    res.to_csv(args.out / "all_sites_screening.csv", index=False)
    accepted.to_csv(args.out / "accepted_proteins.csv", index=False)
    accepted.to_excel(args.out / "accepted_proteins.xlsx", index=False)
    print("\nSite status counts:\n" + res["status"].value_counts().to_string())
    print(f"\nAccepted proteins: {all_ok.sum()} of {len(all_ok)}")
    if len(accepted):
        with pd.option_context("display.width", 220, "display.max_columns", None):
            print(accepted[["uniProtkbId", "primaryAccession", "cofactor", "ligand_label", "pdb_id",
                            "identity", "local_rmsd", "tcs", "alphafill_structure"]].to_string(index=False))
    print(f"\nWrote tables and structures to {args.out}/")


if __name__ == "__main__":
    main()