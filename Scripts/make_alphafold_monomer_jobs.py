import json
import re
from pathlib import Path

PEC_FILE = Path("final_pec_dataset.json")
TARGET_FILE = Path("final_target_proteins_dataset.json")
OUTPUT_FILE = Path("all_monomer_jobs.json")

def load_records(path):
    data = json.loads(path.read_text(encoding="utf-8"))

    if isinstance(data, list):
        return data

    if isinstance(data, dict):
        for key in ("data", "results", "proteins", "targets", "pecs", "entries"):
            if isinstance(data.get(key), list):
                return data[key]

        # Some datasets are dictionaries keyed by protein name.
        return [
            {"name": key, **value} if isinstance(value, dict)
            else {"name": key, "sequence": value}
            for key, value in data.items()
        ]

    raise ValueError(f"Unsupported JSON format: {path}")


def get_value(record, keys):
    for key in keys:
        value = record.get(key)
        if value is not None and str(value).strip():
            return value
    return None


jobs = []
seen_sequences = set()
skipped = 0

for category, filename in [
    ("PEC", PEC_FILE),
    ("target", TARGET_FILE),
]:
    for i, record in enumerate(load_records(filename), start=1):
        name = get_value(record, [
            "protein_name", "proteinName", "uniProtkbId",
            "entry_name", "protein", "name", "id",
            "accession", "primaryAccession", "ID"
        ])

        sequence = get_value(record, [
            "sequence", "Sequence", "seq",
            "protein_sequence", "amino_acid_sequence"
        ])

        if not sequence:
            print(f"Skipping {category} record {i}: no sequence")
            skipped += 1
            continue

        sequence = re.sub(r"\s+", "", str(sequence)).upper()

        if not re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWYBXZJUO]+", sequence):
            print(f"Skipping {name or category + str(i)}: invalid sequence")
            skipped += 1
            continue

        # Remove duplicate sequences
        if sequence in seen_sequences:
            continue
        seen_sequences.add(sequence)

        name = str(name or f"{category}_{i:04d}")
        name = re.sub(r"[^A-Za-z0-9_.-]+", "_", name)

        jobs.append({
            "name": name,
            "modelSeeds": [],
            "sequences": [
                {
                    "proteinChain": {
                        "sequence": sequence,
                        "count": 1
                    }
                }
            ],
            "dialect": "alphafoldserver",
            "version": 1
        })

OUTPUT_FILE.write_text(
    json.dumps(jobs, indent=2) + "\n",
    encoding="utf-8"
)

print(f"Created {OUTPUT_FILE}: {len(jobs)} unique monomer jobs")
print(f"Skipped records: {skipped}")