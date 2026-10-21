# This script generates jobs.json for the AlphaFold server.
# the jobs.json file will later be submitted to the AlphaFold server to run the structure prediction for protein monomers without cofactors.
import json
from pathlib import Path


BASE_DIR = Path(__file__).resolve().parent
PROJECT_DIR = BASE_DIR.parent.parent
INPUT_DIR = PROJECT_DIR / "Inputs" / "json_files" # "combined_filtered_proteins.json" is in INPUT_DIR
OUTPUT_DIR = PROJECT_DIR / "Inputs" / "json_files" / "AlphaFold_jobs" # the output file "jobs.json" will be written to OUTPUT_DIR
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

jobs = []


for entry in json.loads((INPUT_DIR / "combined_filtered_proteins.json").read_text()):
    jobs.append({
        "name": entry["uniProtkbId"],
        "modelSeeds": [],
        "sequences": [{"proteinChain": {"sequence": entry["sequence"]["value"], "count": 1}}],
        "dialect": "alphafoldserver",
        "version": 1,
    })

(OUTPUT_DIR / "jobs.json").write_text(json.dumps(jobs, indent=2))
print(len(jobs), "jobs written")


'''Path("jobs.json").write_text(json.dumps(jobs, indent=2))
print(len(jobs), "jobs written")'''