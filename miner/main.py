# orchestrates workflow
import json
import os
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

    # =========================
    # DNA SEARCHES
    # =========================
    archive_dna_ids = set(run_search(config["archive_search_dna"], config["email"]))
    recent_dna_ids  = set(run_search(config["weekly_search_dna"], config["email"]))

    # =========================
    # RNA SEARCHES
    # =========================
    archive_rna_ids = set(run_search(config["archive_search_rna"], config["email"]))
    recent_rna_ids  = set(run_search(config["weekly_search_rna"], config["email"]))

    # =========================
    # COMBINE ALL STUDIES
    # =========================
    archive_ids = archive_dna_ids.union(archive_rna_ids)
    recent_ids  = recent_dna_ids.union(recent_rna_ids)

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
    # BUILD REPORT (UNCHANGED)
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

    # =========================
    # SPLIT DNA vs RNA (FINAL SETS)
    # =========================
    dna_all = archive_dna_ids.union(recent_dna_ids)
    rna_all = archive_rna_ids.union(recent_rna_ids)

    # -----------------------
    # EXPORT JSON FOR R (NAMED LISTS)
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
