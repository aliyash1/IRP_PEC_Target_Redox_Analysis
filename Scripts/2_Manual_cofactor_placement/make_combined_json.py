import json
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
PROJECT_DIR = BASE_DIR.parent.parent
INPUT_DIR = PROJECT_DIR / "Inputs" / "json_files"
OUTPUT_FILE = INPUT_DIR / "combined_filtered_proteins.json"

# Target proteins, in the desired output order
NAMES = """
FER NRFC YCCM YSAA YKGF YGFS NORV HCR YFAE BFD YFHL HYDN ALKB
ILVD TSAD ISCS DEF FNR PFLA ADHE FUCO ISCA END3 BFR FUMA ERPA MIAB SODF ISCR
NADA BIOB FUMB DMSA NARY NARZ HYPD ACNA CODA DING ROXA THIH MOAA HEMN NAPA
ACNB TAUD FHUF FADH LIPA ISPG ISPH NFUA HCP YCBX TTCA YOAA GLAH D2HDH HYBA
NAPF EDD RIMO MUTY THIC EPMB QUEG LUXS HEMW QUEE YDHV EUTG YNFE YNFF YEIL
YHCC CPDA NSRR PFLC AEGA YJJW PFLE HGLS YDEP UACF ASLB YJIL YBHJ YBIX
YDEM YGIQ
""".split()


SOURCE_FILES = ["final_target_proteins_dataset.json", "final_pec_dataset.json"]

# Merge both files, de-duplicating by accession (keep the entry with more fields,
# so the target-proteins version with cofactors/ec_number wins over the pec one)
merged = {}
for fname in SOURCE_FILES:
    for entry in json.loads((INPUT_DIR / fname).read_text()):
        acc = entry["primaryAccession"]
        if acc not in merged or len(entry) > len(merged[acc]):
            merged[acc] = entry

# Index by ID prefix (the part before _ECOLI)
by_name = {}
for entry in merged.values():
    by_name.setdefault(entry["uniProtkbId"].rsplit("_", 1)[0], []).append(entry)

# Keep only the listed proteins, in list order
out, missing = [], []

for name in NAMES:
    if name in by_name:
        out.extend(by_name[name])
    else:
        missing.append(name)

OUTPUT_FILE.write_text(json.dumps(out, indent=2))
print(len(out), "proteins written to", OUTPUT_FILE)
if missing:
    print("Not found:", missing)