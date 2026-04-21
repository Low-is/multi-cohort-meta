import os
import csv

# =========================
# FILE PATHS
# =========================
GSE_PATH = "data/gse_ids.csv"
REPORT_PATH = "outputs/weekly_reports.csv"


# =========================
# LOAD GSE IDS (ARCHIVE)
# =========================
def load_seen_ids():
    if not os.path.exists(GSE_PATH):
        return set()

    with open(GSE_PATH, "r") as f:
        return set(line.strip() for line in f if line.strip())


# =========================
# SAVE FULL ARCHIVE (IDs ONLY)
# =========================
def save_seen_ids(all_ids):
    """
    FULL HISTORY STORE:
    - 2014 → present
    - no duplicates
    - always overwritten snapshot
    """

    os.makedirs(os.path.dirname(GSE_PATH), exist_ok=True)

    all_ids = sorted(set(all_ids))

    with open(GSE_PATH, "w") as f:
        for gse in all_ids:
            f.write(f"{gse}\n")


# =========================
# LOAD REPORT IDS (IDs ONLY)
# =========================
def load_report_ids():
    if not os.path.exists(REPORT_PATH):
        return set()

    with open(REPORT_PATH, "r", newline="") as f:
        reader = csv.reader(f)
        next(reader, None)  # skip header
        return set(row[0] for row in reader if row)


# =========================
# SAVE FULL REPORT (IDS + LINKS)
# =========================
def save_weekly_report(all_ids):
    """
    FULL REPORT STORE:
    - IDs + GEO links
    - no duplicates
    - overwritten snapshot
    """

    os.makedirs(os.path.dirname(REPORT_PATH), exist_ok=True)

    all_ids = sorted(set(all_ids))

    with open(REPORT_PATH, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["GSE_ID", "LINK"])

        for gse in all_ids:
            writer.writerow([
                gse,
                f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}"
            ])
