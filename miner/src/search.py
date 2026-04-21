from Bio import Entrez

def run_search(search_cfg, email):

    Entrez.email = email

    # -------------------------
    # BUILD SEARCH PARAMETERS
    # -------------------------
    esearch_params = {
        "db": search_cfg["database"],
        "term": search_cfg["query"],
        "retmax": search_cfg["retmax"],
    }

    # ONLY ADD reldate IF IT EXISTS
    if "reldate" in search_cfg:
        esearch_params["reldate"] = search_cfg["reldate"]
        esearch_params["datetype"] = "pdat"

    # -------------------------
    # SEARCH GEO DATASETS
    # -------------------------
    handle = Entrez.esearch(**esearch_params)
    record = Entrez.read(handle)
    handle.close()

    gse_ids = record.get("IdList", [])

    if not gse_ids:
        return []

    # -------------------------
    # FETCH SUMMARIES
    # -------------------------
    handle = Entrez.esummary(
        db=search_cfg["database"],
        id=",".join(gse_ids)
    )

    summaries = Entrez.read(handle)
    handle.close()

    # -------------------------
    # PARSE RESULTS
    # -------------------------
    gse_list = []

    if isinstance(summaries, dict) and "DocumentSummarySet" in summaries:
        docs = summaries["DocumentSummarySet"]["DocumentSummary"]
    else:
        docs = summaries

    for doc in docs:
        acc = doc.get("Accession", "")
        if acc.startswith("GSE"):
            gse_list.append(acc)

    return gse_list
