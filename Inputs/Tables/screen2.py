from __future__ import annotations
#!/usr/bin/env python3
"""
Permissive AlphaFill screen for protein/cofactor binding sites.

Inputs:
  Original long Excel table: columns such as primaryAccession, uniProtkbId,
    description, position, and optionally ligand_label.
   Project Excel table: Protein, Cofactor (label) / Cofactor(s), Binding-site residues.
   JSON list: records with protein/accession and cofactors [{name, label, residues}, ...].

APIs:
  AlphaFill metadata: GET https://alphafill.eu/v1/aff/{accession}/json
  AlphaFill stripped structure: GET https://alphafill.eu/v1/aff/{accession}/stripped/{asymlist}
  UniProt entry lookup: GET https://rest.uniprot.org/uniprotkb/{entry_name}.json?fields=accession
  UniProt search fallback: GET https://rest.uniprot.org/uniprotkb/search

"""


'''Input workbook: Inputs/Tables/Binding_sites_with_FeS_cofactor.xlsx

Script location: Scripts/1_alphafill/

Results output: Results/1_alphafill/

Cached AlphaFill JSON: Results/1_alphafill/json/

Downloaded structures: Results/1_alphafill/structures/
'''

import argparse
import json
import re
import time
from pathlib import Path
from typing import Any

import pandas as pd
import requests

ALPHAFILL_DEFAULT = "https://alphafill.eu/v1/aff"
UNIPROT = "https://rest.uniprot.org/uniprotkb"

# thresholds
MIN_IDENTITY = 0.50
MAX_LOCAL_RMSD = 3.10 
MAX_TCS = 2.00         

# AlphaFill analogue_id/compound_id codes for common Fe-S cluster descriptions.
ACCEPTED_CODES = {
    "SF4": {"SF4"},
    "FES": {"FES"},
    "F3S": {"F3S"},
    "FES_ANY": {"SF4", "FES", "F3S"},
}


def canonical_cofactor(name: Any) -> str:
    t = re.sub(r"[\s_\[\]\(\)]", "", str(name).upper())
    if "4FE-2O-2S" in t or "4FE2O2S" in t:
        return "UNSUPPORTED"
    if "4FE-4S" in t or "4FE4S" in t or t == "SF4":
        return "SF4"
    if "2FE-2S" in t or "2FE2S" in t or t == "FES":
        return "FES"
    if "3FE-4S" in t or "3FE4S" in t or t == "F3S":
        return "F3S"
    if "IRON-SULFUR" in t or "IRONSULFUR" in t or "FE-S" in t:
        return "FES_ANY"
    return "OTHER"


def safe_name(value: Any) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_") or "unknown"


def accession_lookup(session: requests.Session, entry_name: str) -> str | None:
    """Resolve a UniProt entry name to its primary accession."""
    name = str(entry_name).strip()
    if not name or name.lower() == "nan":
        return None

    try:
        r = session.get(f"{UNIPROT}/{name}.json", params={"fields": "accession"}, timeout=30)
        if r.ok:
            result = r.json().get("primaryAccession")
            if result:
                return result
    except (requests.RequestException, ValueError):
        pass

    try:
        r = session.get(
            f"{UNIPROT}/search",
            params={"query": f"id:{name}", "fields": "accession", "size": 1, "format": "json"},
            timeout=30,
        )
        if r.ok:
            results = r.json().get("results", [])
            if results:
                return results[0].get("primaryAccession")
    except (requests.RequestException, ValueError):
        pass
    return None


def parse_residues(value: Any) -> list[int]:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return []
    if isinstance(value, (list, tuple, set)):
        vals = []
        for item in value:
            vals.extend(parse_residues(item))
        return sorted(set(vals))
    # Residue cells may be e.g. "[12, 15, 18]" or "12, 15, 18".
    return sorted({int(x) for x in re.findall(r"\d+", str(value))})


def split_cofactor_label(value: Any) -> tuple[str, str | None]:
    """Split text like '[4Fe-4S] cluster (label 2)' into name and label."""
    text = str(value).strip()
    m = re.search(r"\s*\(\s*label\s+([^)]+)\s*\)\s*$", text, flags=re.I)
    if m:
        return text[:m.start()].strip(), m.group(1).strip()
    return text, None


def read_input(path: Path, session: requests.Session, delay: float) -> pd.DataFrame:
    """Normalise supported input files into one row per site residue."""
    if not path.exists():
        raise FileNotFoundError(
            f"Input file not found: {path}\n"
            f"Current working directory: {Path.cwd()}\n"
            "Pass the full path or run the script from the directory containing the file."
        )

    rows: list[dict[str, Any]] = []
    if path.suffix.lower() == ".json":
        with path.open(encoding="utf-8") as handle:
            obj = json.load(handle)
        if isinstance(obj, dict):
            # Accept either a top-level list under a common key or a mapping of proteins.
            for key in ("proteins", "data", "records"):
                if isinstance(obj.get(key), list):
                    obj = obj[key]
                    break
            else:
                obj = [{"protein": k, **(v if isinstance(v, dict) else {"cofactors": v})}
                       for k, v in obj.items()]
        if not isinstance(obj, list):
            raise ValueError("JSON input must be a list of protein records or a mapping of proteins.")
        for p in obj:
            if not isinstance(p, dict):
                continue
            protein = p.get("protein") or p.get("Protein") or p.get("uniProtkbId") or p.get("id") or p.get("name")
            accession = p.get("accession") or p.get("primaryAccession")
            cofactors = p.get("cofactors", [])
            # Also tolerate a single cofactor entry stored directly on a record.
            if not cofactors and any(k in p for k in ("residues", "binding_site_residues")):
                cofactors = [p]
            for c in cofactors if isinstance(cofactors, list) else []:
                if not isinstance(c, dict):
                    continue
                cofactor = c.get("name") or c.get("cofactor") or c.get("description") or "Unknown cofactor"
                label = c.get("label", c.get("ligand_label"))
                residues = parse_residues(c.get("residues", c.get("binding_site_residues", c.get("position"))))
                for residue in residues:
                    rows.append({
                        "protein_name": protein, "accession": accession, "cofactor": cofactor,
                        "label": str(label) if label is not None else None, "residue": residue,
                    })
    else:
        df = pd.read_excel(path, engine="openpyxl")
        # Case-insensitive alias mapping.
        colmap = {str(c).strip().lower(): c for c in df.columns}

        def find_col(*aliases):
            for alias in aliases:
                if alias.lower() in colmap:
                    return colmap[alias.lower()]
            return None

        accession_col = find_col("primaryAccession", "accession", "uniprot accession")
        name_col = find_col("uniProtkbId", "Protein", "protein", "protein name", "entry name")
        cofactor_col = find_col("description", "Cofactor", "Cofactor(s)", "cofactor (label)", "cofactor name")
        label_col = find_col("ligand_label", "label", "cofactor label")
        residue_col = find_col("position", "residue", "residues", "binding-site residues", "binding site residues")

        if not name_col and not accession_col:
            raise ValueError(f"Cannot identify a protein/accession column in {list(df.columns)}")
        if not cofactor_col or not residue_col:
            raise ValueError(
                f"Cannot identify cofactor/residue columns in {list(df.columns)}. "
                "Expected description/Cofactor and position/Binding-site residues."
            )

        for _, row in df.iterrows():
            protein = row[name_col] if name_col else None
            accession = row[accession_col] if accession_col else None
            if pd.isna(protein):
                protein = accession
            if pd.isna(accession):
                accession = None
            raw_cofactor = row[cofactor_col]
            if pd.isna(raw_cofactor):
                continue
            cofactor, embedded_label = split_cofactor_label(raw_cofactor)
            label = row[label_col] if label_col else embedded_label
            if label is not None and pd.isna(label):
                label = embedded_label
            residues = parse_residues(row[residue_col])
            for residue in residues:
                rows.append({
                    "protein_name": protein, "accession": accession, "cofactor": cofactor,
                    "label": str(label) if label is not None else None, "residue": residue,
                })

    out = pd.DataFrame(rows, columns=["protein_name", "accession", "cofactor", "label", "residue"])
    if out.empty:
        raise ValueError("No binding-site residues were read. Check the input columns and file contents.")

    # Resolve missing accessions from UniProt entry names, cache lookups in memory.
    missing = out["accession"].isna() | out["accession"].astype(str).str.strip().eq("")
    names = out.loc[missing, "protein_name"].dropna().astype(str).unique()
    lookup = {}
    for name in names:
        lookup[name] = accession_lookup(session, name)
        time.sleep(delay)
    out.loc[missing, "accession"] = out.loc[missing, "protein_name"].astype(str).map(lookup)
    unresolved = out["accession"].isna() | out["accession"].astype(str).str.strip().eq("")
    if unresolved.any():
        examples = out.loc[unresolved, "protein_name"].drop_duplicates().head(10).tolist()
        raise ValueError(
            f"Could not resolve UniProt accession for {unresolved.sum()} residue rows. "
            f"Examples: {examples}. Add primaryAccession to the input workbook."
        )
    out["accession"] = out["accession"].astype(str).str.strip()
    out["residue"] = out["residue"].astype(int)
    return out


def build_sites(df: pd.DataFrame, include_all: bool) -> pd.DataFrame:
    df = df.copy()
    df["canonical"] = df["cofactor"].map(canonical_cofactor)
    keys = ["accession", "protein_name", "cofactor", "canonical", "label"]
    sites = (
        df.groupby(keys, dropna=False)["residue"]
        .apply(lambda s: set(int(x) for x in s))
        .rename("expected_residues")
        .reset_index()
    )
    if not include_all:
        sites = sites[sites["canonical"].isin(ACCEPTED_CODES)]
    if sites.empty:
        raise ValueError("No Fe-S sites remain after cofactor filtering.")
    return sites.reset_index(drop=True)


def request(session: requests.Session, url: str, *, want_json: bool, retries: int = 4):
    last_error = None
    for attempt in range(retries):
        try:
            response = session.get(url, timeout=120)
            if response.status_code == 404:
                return None
            if response.status_code == 429 or response.status_code >= 500:
                wait = response.headers.get("Retry-After")
                try:
                    wait = float(wait) if wait is not None else 2 ** attempt
                except ValueError:
                    wait = 2 ** attempt
                time.sleep(min(wait, 30))
                continue
            response.raise_for_status()
            return response.json() if want_json else response.content
        except (requests.RequestException, ValueError) as exc:
            last_error = exc
            time.sleep(min(2 ** attempt, 30))
    raise RuntimeError(f"Request failed for {url}: {last_error}")


def as_float(value: Any) -> float | None:
    try:
        return float(value) if value is not None else None
    except (TypeError, ValueError):
        return None


def flatten_transplants(data: dict) -> list[dict]:
    output = []
    for hit in data.get("hits") or []:
        alignment = hit.get("alignment") or {}
        identity = as_float(alignment.get("identity"))
        if identity is not None and identity > 1:
            identity /= 100.0
        for transplant in hit.get("transplants") or []:
            clash = transplant.get("clash") or {}
            distances = clash.get("distances") or []
            residues = set()
            for distance in distances:
                atom = distance.get("poly_atom") or {}
                seq_id = atom.get("seq_id")
                if seq_id is not None:
                    try:
                        residues.add(int(float(seq_id)))
                    except (TypeError, ValueError):
                        pass
            output.append({
                "pdb_id": hit.get("pdb_id"),
                "pdb_asym_id": hit.get("pdb_asym_id"),
                "asym_id": transplant.get("asym_id"),
                "compound_id": str(transplant.get("compound_id") or "").upper(),
                "analogue_id": str(transplant.get("analogue_id") or "").upper(),
                "identity": identity,
                "local_rmsd": as_float(transplant.get("local_rmsd")),
                "tcs": as_float(clash.get("score")),
                "alphafill_residues": residues,
            })
    return output


def quality_failures(candidate: dict, min_identity: float, max_rmsd: float, max_tcs: float) -> list[str]:
    failures = []
    if candidate["identity"] is None or candidate["identity"] < min_identity:
        failures.append(f"identity<{min_identity:.2f}")
    if candidate["local_rmsd"] is None or candidate["local_rmsd"] > max_rmsd:
        failures.append(f"local_RMSD>{max_rmsd:.2f}")
    if candidate["tcs"] is None or candidate["tcs"] > max_tcs:
        failures.append(f"TCS>{max_tcs:.2f}")
    return failures


def residue_comparison(expected: set[int], observed: set[int], mode: str, min_overlap: float):
    intersection = expected & observed
    union = expected | observed
    recall = len(intersection) / len(expected) if expected else 0.0
    precision = len(intersection) / len(observed) if observed else 0.0
    jaccard = len(intersection) / len(union) if union else 0.0
    if mode == "none":
        okay = True
    elif mode == "exact":
        okay = expected == observed
    elif mode == "superset":
        okay = expected <= observed
    else:  # overlap
        # Require at least one shared residue and sufficient coverage of the expected set.
        okay = bool(intersection) and recall >= min_overlap
    return okay, recall, precision, jaccard, sorted(expected - observed), sorted(observed - expected)


def screen_site(site, candidates, args, used):
    expected = site["expected_residues"]
    codes = ACCEPTED_CODES.get(site["canonical"], {site["canonical"]})
    right = [c for c in candidates if c["analogue_id"] in codes or c["compound_id"] in codes]
    base = {"n_transplants_correct_cofactor": len(right)}
    if not right:
        return {**base, "status": "REJECT", "reason": "no AlphaFill transplant with matching cofactor code"}

    evaluated = []
    for candidate in right:
        match = residue_comparison(expected, candidate["alphafill_residues"], args.match_mode, args.min_overlap)
        failures = quality_failures(candidate, args.min_identity, args.max_local_rmsd, args.max_tcs)
        evaluated.append((candidate, match, failures))

    passing = [
        item for item in evaluated
        if not item[2] and item[1][0] and (item[0]["pdb_id"], item[0]["asym_id"]) not in used
    ]
    # Prefer better residue overlap, then lower local RMSD/TCS, then higher identity.
    passing.sort(key=lambda x: (-x[1][3], x[0]["local_rmsd"], x[0]["tcs"], -(x[0]["identity"] or 0)))
    if passing:
        candidate, match, failures = passing[0]
        status = "ACCEPT"
        reason = "cofactor, quality thresholds and residue matching pass"
        used.add((candidate["pdb_id"], candidate["asym_id"]))
    else:
        # Keep the best candidate for review; don't call a poor-metric candidate a biological reject.
        evaluated.sort(
            key=lambda x: (
                len(x[2]) + (0 if x[1][0] else 1),
                -x[1][3],
                x[0]["local_rmsd"] if x[0]["local_rmsd"] is not None else float("inf"),
                x[0]["tcs"] if x[0]["tcs"] is not None else float("inf"),
            )
        )
        candidate, match, failures = evaluated[0]
        if not match[0]:
            failures.append("residue_match_failed")
        status = "REVIEW"
        reason = "; ".join(failures) if failures else "candidate needs manual review"

    _, recall, precision, jaccard, missing, extra = match
    return {
        **base,
        "status": status,
        "reason": reason,
        "pdb_id": candidate["pdb_id"],
        "pdb_asym_id": candidate["pdb_asym_id"],
        "asym_id": candidate["asym_id"],
        "compound_id": candidate["compound_id"],
        "analogue_id": candidate["analogue_id"],
        "identity": candidate["identity"],
        "local_rmsd": candidate["local_rmsd"],
        "tcs": candidate["tcs"],
        "alphafill_residues": sorted(candidate["alphafill_residues"]),
        "residue_recall": recall,
        "residue_precision": precision,
        "residue_jaccard": jaccard,
        "missing_from_alphafill": missing,
        "extra_in_alphafill": extra,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("table", nargs="?", type=Path, default=Path("Binding_sites_with_FeS_cofactor.xlsx"))
    parser.add_argument("--out", type=Path, default=Path("alphafill_results"))
    parser.add_argument("--match-mode", choices=["overlap", "exact", "superset", "none"], default="overlap")
    parser.add_argument("--min-overlap", type=float, default=0.50,
                        help="Minimum fraction of expected residues found in AlphaFill contact residues (overlap mode).")
    parser.add_argument("--min-identity", type=float, default=MIN_IDENTITY,
                        help="Minimum sequence identity as a fraction (default 0.50 = 50%%).")
    parser.add_argument("--max-local-rmsd", type=float, default=MAX_LOCAL_RMSD,
                        help="Maximum local RMSD in Angstrom (default 3.10).")
    parser.add_argument("--max-tcs", type=float, default=MAX_TCS,
                        help="Maximum TCS score (default 2.00; provisional/permissive).")
    parser.add_argument("--include-all-cofactors", action="store_true")
    parser.add_argument("--base-url", default=ALPHAFILL_DEFAULT)
    parser.add_argument("--delay", type=float, default=0.5)
    args = parser.parse_args()

    if not 0 <= args.min_identity <= 1:
        parser.error("--min-identity must be between 0 and 1, e.g. 0.50 for 50%.")
    if not 0 <= args.min_overlap <= 1:
        parser.error("--min-overlap must be between 0 and 1.")
    args.out.mkdir(parents=True, exist_ok=True)
    json_dir = args.out / "json"
    structure_dir = args.out / "structures"
    json_dir.mkdir(parents=True, exist_ok=True)
    structure_dir.mkdir(parents=True, exist_ok=True)

    session = requests.Session()
    session.headers.update({"User-Agent": "alphafill-screen/2.0 (research use)"})
    df = read_input(args.table, session, args.delay)
    sites = build_sites(df, args.include_all_cofactors)
    print(f"Input: {args.table.resolve()}")
    print(f"Sites to screen: {len(sites)} across {sites['accession'].nunique()} proteins")
    print(f"Criteria: identity >= {args.min_identity:.2f}, local RMSD <= {args.max_local_rmsd:.2f} A, "
          f"TCS <= {args.max_tcs:.2f}; residue mode={args.match_mode}")

    results = []
    for accession, group in sites.groupby("accession", sort=True):
        cache_path = json_dir / f"{safe_name(accession)}.json"
        try:
            if cache_path.exists():
                with cache_path.open(encoding="utf-8") as handle:
                    metadata = json.load(handle)
            else:
                time.sleep(args.delay)
                metadata = request(session, f"{args.base_url.rstrip('/')}/{accession}/json", want_json=True)
                if metadata is not None:
                    cache_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
        except (RuntimeError, OSError, ValueError) as exc:
            for _, site in group.iterrows():
                results.append({**site.to_dict(), "status": "REVIEW",
                                "reason": f"AlphaFill API/cache error: {exc}"})
            print(f"{accession}: API error (review)")
            continue

        if metadata is None:
            for _, site in group.iterrows():
                results.append({**site.to_dict(), "status": "REVIEW",
                                "reason": "AlphaFill returned 404/no metadata; not a biological rejection"})
            print(f"{accession}: no AlphaFill metadata (review)")
            continue

        candidates = flatten_transplants(metadata)
        used = set()
        for _, site in group.iterrows():
            screened = screen_site(site, candidates, args, used)
            results.append({**site.to_dict(), **screened})
        print(f"{accession}: " + ", ".join(r["status"] for r in results[-len(group):]))

    res = pd.DataFrame(results)
    if res.empty:
        raise RuntimeError("No results were produced.")
    res["expected_residues"] = res["expected_residues"].map(
        lambda x: ",".join(map(str, sorted(x))) if isinstance(x, set) else x
    )
    for col in ("alphafill_residues", "missing_from_alphafill", "extra_in_alphafill"):
        if col in res.columns:
            res[col] = res[col].map(
                lambda x: ",".join(map(str, x)) if isinstance(x, (list, set, tuple)) else x
            )
    res["label"] = res["label"].where(res["label"].notna(), "")

    # Download stripped structures only for proteins whose every site was accepted.
    accepted_proteins = []
    status_by_accession = res.groupby("accession")["status"].agg(lambda s: set(s) == {"ACCEPT"})
    for accession in status_by_accession[status_by_accession].index:
        sub = res[res["accession"] == accession]
        asyms = sorted(set(str(x) for x in sub["asym_id"].dropna() if str(x).strip()))
        if not asyms:
            continue
        path = structure_dir / f"{safe_name(accession)}_{'_'.join(map(safe_name, asyms))}_stripped.cif"
        try:
            time.sleep(args.delay)
            content = request(
                session,
                f"{args.base_url.rstrip('/')}/{accession}/stripped/{','.join(asyms)}",
                want_json=False,
            )
            if not content:
                download_status = "no structure returned"
                structure_path = ""
            else:
                # AlphaFill structure endpoints are documented as mmCIF. Save response but flag unexpected content.
                text = content[:1000].decode("utf-8", errors="ignore").lstrip()
                if not (text.startswith("data_") or "_atom_site." in text):
                    download_status = "downloaded; check response format"
                else:
                    download_status = "downloaded"
                path.write_bytes(content)
                structure_path = str(path.resolve())
        except (RuntimeError, OSError) as exc:
            download_status = f"download failed: {exc}"
            structure_path = ""
        for _, row in sub.iterrows():
            accepted_proteins.append({
                "protein_name": row["protein_name"], "accession": accession,
                "cofactor": row["cofactor"], "label": row["label"],
                "expected_residues": row["expected_residues"],
                "alphafill_residues": row.get("alphafill_residues", ""),
                "pdb_id": row.get("pdb_id", ""), "identity": row.get("identity", ""),
                "local_rmsd": row.get("local_rmsd", ""), "tcs": row.get("tcs", ""),
                "asym_id": row.get("asym_id", ""), "structure_file": structure_path,
                "download_status": download_status,
            })

    res.to_csv(args.out / "all_sites_screening.csv", index=False)
    res.to_excel(args.out / "all_sites_screening.xlsx", index=False)
    accepted = pd.DataFrame(accepted_proteins)
    accepted.to_csv(args.out / "accepted_proteins.csv", index=False)
    accepted.to_excel(args.out / "accepted_proteins.xlsx", index=False)

    print("\nSite status counts:")
    print(res["status"].value_counts(dropna=False).to_string())
    print(f"\nProteins with all sites ACCEPTed: {len(accepted_proteins and accepted['accession'].nunique() or [])}")
    print(f"Results directory: {args.out.resolve()}")
    print("Wrote all_sites_screening.csv/.xlsx, accepted_proteins.csv/.xlsx, json/, structures/")


if __name__ == "__main__":
    main()
