import os
import csv

REPORT_PATH = "outputs/weekly_report.csv"

def append_to_weekly_report(gse_ids):
    os.makedirs(os.path.dirname(REPORT_PATH), exist_ok=True)

    # Load existing rows to prevent duplicates
    existing = set()

    if os.path.isfile(REPORT_PATH):
        with open(REPORT_PATH, "r", newline="") as f:
            reader = csv.reader(f)
            next(reader, None)  # skip header
            for row in reader:
                if row:
                    existing.add(row[0])

    write_header = not os.path.isfile(REPORT_PATH)

    with open(REPORT_PATH, "a", newline="") as f:
        writer = csv.writer(f)

        if write_header:
            writer.writerow(["GSE_ID", "Link"])

        for gse_id in gse_ids:
            if gse_id in existing:
                continue

            link = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse_id}"
            writer.writerow([gse_id, link])
            existing.add(gse_id)
