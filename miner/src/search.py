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
    
    # Made changes 6-25-2026
    if isinstance(summaries, list):
        docs = summaries
    else:
        docs = summaries.get("DocumentSummarySet", {}).get("DocumentSummary", [])

    gse_list = []
    for doc in docs:
        acc = doc.get("Accession", "") or doc.get("accession", "")
        if not str(acc).startswith("GSE"):
            continue

        gse_list.append({
            "gse": acc,
            "title": doc.get("title", "") or doc.get("Title", ""),
            "summary": doc.get("summary", "") or doc.get("Summary", ""),
            "type": doc.get("type", "") or doc.get("Type", ""),
            "overall_design": doc.get("overall_design", "") or doc.get("Overall_Design", "")
        })

    #if isinstance(summaries, dict) and "DocumentSummarySet" in summaries:
        #docs = summaries["DocumentSummarySet"]["DocumentSummary"]
    #else:
        #docs = summaries

    #for doc in docs:
        #acc = doc.get("Accession", "")
        #if acc.startswith("GSE"):
            #gse_list.append(acc)

    return gse_list
