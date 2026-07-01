# orchestrates workflow
import json
import os
import re
import yaml
import csv
from src.search import run_search
from src.storage import load_seen_ids, save_seen_ids

ARCHIVE_PATH = "data/gse_ids.csv"
REPORT_PATH = "outputs/weekly_report.csv"

# -----------------------
# FILTER RESULTS
# -----------------------
def normalize(text):
    text = text.lower()
    text = text.replace("-", " ")
    text = re.sub(r"[^\w\s]", " ", text)
    text = re.sub(r"\s+", " ", text)
    return text.strip()
    
def keep_study(study, config):
    
    text = normalize(
        study.get("title", "") + " " + study.get("summary", "") + " " + study.get("type", "")
    )

    if "sepsis" not in title:
        return False
   
    # hard exclusion first
    if any(k in text for k in config["exclude_keys"]):
        return False

    # DNA microarray
    if any(k in text for k in config["dna_keys"]):
        return True

    # bulk RNA-seq
    if any(k in text for k in config["rna_keys"]):
        return True
    return False

# -----------------------
# SAVE REPORT (ALL DATA + STATUS)
# -----------------------
def save_report(rows):
    os.makedirs(os.path.dirname(REPORT_PATH), exist_ok=True)

    with open(REPORT_PATH, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["GSE_ID", "LINK", "STATUS"])
        writer.writerows(rows)


# -----------------------
# MAIN PIPELINE
# -----------------------
def main():
    with open("configs/config.yaml", "r") as f:
        config = yaml.safe_load(f)

    print("🔍 Running GEO search...")

    # =========================================================
    # DNA SEARCHES (UPDATED 6-25-2026)
    # =========================================================
    archive_dna = run_search(config["dna_archive_search"], config["email"])
    recent_dna  = run_search(config["dna_weekly_search"], config["email"])

    # OLD CODE (COMMENTED OUT - PREVIOUS PIPELINE)
    # archive_dna_ids = set(run_search(config["dna_archive_search"], config["email"]))
    # recent_dna_ids  = set(run_search(config["dna_weekly_search"], config["email"]))

    # =========================================================
    # RNA SEARCHES (UPDATED 6-25-2026)
    # =========================================================
    archive_rna = run_search(config["rna_archive_search"], config["email"])
    recent_rna  = run_search(config["rna_weekly_search"], config["email"])

    # OLD CODE (COMMENTED OUT - PREVIOUS PIPELINE)
    # archive_rna_ids = set(run_search(config["rna_archive_search"], config["email"]))
    # recent_rna_ids  = set(run_search(config["rna_weekly_search"], config["email"]))

    # =========================================================
    # APPLY FILTER (TITLE + SUMMARY ONLY) | added 7-1-2026
    # =========================================================
    archive_dna = [x for x in archive_dna if keep_study(x, config)]
    recent_dna  = [x for x in recent_dna if keep_study(x, config)]

    archive_rna = [x for x in archive_rna if keep_study(x, config)]
    recent_rna  = [x for x in recent_rna if keep_study(x, config)]
    

    # =========================================================
    # COMBINE ALL STUDIES (UPDATED 6-25-2026)
    # =========================================================
    archive_ids = {
        x["gse"] for x in archive_dna + archive_rna
    }

    recent_ids = {
        x["gse"] for x in recent_dna + recent_rna
    }

    # OLD CODE (COMMENTED OUT - PREVIOUS PIPELINE)
    # archive_ids = archive_dna_ids.union(archive_rna_ids)
    # recent_ids  = recent_dna_ids.union(recent_rna_ids)

    # -----------------------
    # LOAD EXISTING ARCHIVE
    # -----------------------
    existing_ids = load_seen_ids()

    # -----------------------
    # UPDATE FULL ARCHIVE
    # -----------------------
    full_archive = existing_ids.union(archive_ids).union(recent_ids)
    save_seen_ids(full_archive)

    # -----------------------
    # BUILD REPORT
    # -----------------------
    rows = []

    for gse in full_archive:
        status = "new" if gse in recent_ids else "old"

        rows.append((
            gse,
            f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}",
            status
        ))

    save_report(rows)

    # =========================================================
    # BUILD TRUE DNA / RNA SETS
    # =========================================================
    dna_all = {
        x["gse"] for x in archive_dna + recent_dna
        if x["gse"] in full_archive
    }

    rna_all = {
        x["gse"] for x in archive_rna + recent_rna
        if x["gse"] in full_archive
    }

    # OLD CODE (COMMENTED OUT - PREVIOUS PIPELINE)
    # dna_detected = archive_dna_ids.union(recent_dna_ids)
    # rna_detected = archive_rna_ids.union(recent_rna_ids)
    # dna_all = full_archive.intersection(dna_detected)
    # rna_all = full_archive.intersection(rna_detected)

    # -----------------------
    # EXPORT JSON FOR R
    # -----------------------
    os.makedirs("outputs", exist_ok=True)

    dna_named = {gse: gse for gse in dna_all}
    rna_named = {gse: gse for gse in rna_all}

    with open("outputs/dna_gse_ids.json", "w") as f:
        json.dump(dna_named, f, indent=2)

    with open("outputs/rna_gse_ids.json", "w") as f:
        json.dump(rna_named, f, indent=2)

    print(f"✔️ Exported {len(dna_named)} DNA studies")
    print(f"✔️ Exported {len(rna_named)} RNA studies")

    # -----------------------
    # SUMMARY
    # -----------------------
    new_ids = recent_ids - existing_ids

    print("\n✔ Pipeline complete")
    print(f"Total archived studies: {len(full_archive)}")
    print(f"New studies (730d): {len(recent_ids)}")


if __name__ == "__main__":
    main()
