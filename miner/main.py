# orchestrates workflow

import yaml
import csv
from src.search import run_search
from src.storage import load_seen_ids, save_seen_ids

ARCHIVE_PATH = "data/gse_ids.csv"
REPORT_PATH = "outputs/weekly_report.csv"


# -----------------------
# SAVE REPORT (ALL DATA + STATUS)
# -----------------------
def save_report(rows):
    import os
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

    # FULL ARCHIVE SEARCH (all time)
    archive_ids = set(run_search(config["archive_search"], config["email"]))

    # RECENT SEARCH (last 730 days)
    recent_ids = set(run_search(config["weekly_search"], config["email"]))

    # -----------------------
    # LOAD EXISTING ARCHIVE
    # -----------------------
    existing_ids = load_seen_ids()

    # -----------------------
    # UPDATE FULL ARCHIVE
    # -----------------------
    full_archive = existing_ids.union(archive_ids).union(recent_ids)
    save_seen_ids(full_archive)   # <-- THIS updates gse_ids.csv

    # -----------------------
    # BUILD REPORT (LABELLED)
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

    # -----------------------
    # SUMMARY
    # -----------------------
    new_ids = recent_ids - existing_ids

    print("\n✔ Pipeline complete")
    print(f"Total archived studies: {len(full_archive)}")
    print(f"New studies (730d): {len(recent_ids)}")


if __name__ == "__main__":
    main()
