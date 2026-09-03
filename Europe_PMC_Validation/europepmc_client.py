"""
europepmc_client.py
====================
Europe PMC REST API client — the drop-in replacement for the PubMed/NCBI
E-utilities calls (Literature_Search.pubmed_search / run_network.py's
_ncbi_get, search_pubmed, fetch_abstracts, _parse_xml_abstracts) used by the
Edge_cosine_met_reference_network / Edge_cosine_general_reference_network
pipelines. This module exists so the query-construction and search/fetch
layer can differ from those pipelines while every downstream step (bootstrap
resampling, same-sentence cross-mention detection, Otsuka-Ochiai scoring,
elbow threshold, key-player extraction, dashboard assembly, LLM
descriptions) stays byte-identical — they only ever consume
{pmid, title, abstract, year} paper dicts, the same shape the PubMed
pipeline produced.

What's different about Europe PMC, and why:
  - One paginated JSON call (resultType=core) returns title + abstract +
    year + id directly. There is no separate esearch -> efetch two-step and
    no XML parsing — search_europepmc() below returns full paper dicts in
    one pass, unlike search_pubmed()+fetch_abstracts().
  - Pagination uses an opaque cursorMark continuation token (starting at
    "*"), not an offset/retstart.
  - No API key, "tool", or "email" identification parameter is required.
  - Query field syntax differs: MESH:"term" instead of term[MeSH Terms],
    (TITLE:"term" OR ABSTRACT:"term") instead of term[Title/Abstract], and
    a FIRST_PDATE:[YYYY TO 3000] clause embedded in the query string instead
    of separate mindate/datetype request parameters.
  - Europe PMC does not publish a hard rate limit for this endpoint, but a
    conservative delay between page requests is kept anyway, consistent
    with the PubMed pipeline's approach.

Reference: https://europepmc.org/RestfulWebService
"""

import time
from datetime import datetime

import requests

EUROPEPMC_SEARCH_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
PAGE_SIZE = 1000  # Europe PMC's maximum page size per request


def _epmc_get(params: dict, retries: int = 3):
    for attempt in range(retries):
        try:
            resp = requests.get(EUROPEPMC_SEARCH_URL, params=params, timeout=30)
            resp.raise_for_status()
            return resp
        except requests.RequestException as e:
            wait = 2 ** attempt
            print(f"    [!] HTTP error (attempt {attempt+1}/{retries}): {e} — retry in {wait}s")
            time.sleep(wait)
    return None


def search_europepmc(query: str, max_results: int, years_back: int,
                     delay: float = 0.4) -> list[dict]:
    """
    Query Europe PMC and return up to max_results papers as
    {"pmid", "title", "abstract", "year"} dicts.

    `query` should already be a fully-formed Europe PMC query string (built
    by the alias/keyword clause helpers below) — this function only adds
    the publication-date filter and handles pagination.
    """
    min_year = datetime.now().year - years_back
    full_query = f"({query}) AND FIRST_PDATE:[{min_year} TO 3000]"

    papers: list[dict] = []
    cursor = "*"
    while len(papers) < max_results:
        page_size = min(PAGE_SIZE, max_results - len(papers))
        params = {
            "query":      full_query,
            "format":     "json",
            "resultType": "core",
            "pageSize":   page_size,
            "cursorMark": cursor,
        }
        resp = _epmc_get(params)
        if resp is None:
            break
        data    = resp.json()
        results = data.get("resultList", {}).get("result", [])
        if not results:
            break

        for r in results:
            pmid = r.get("pmid") or r.get("id", "")
            if not pmid:
                continue
            papers.append({
                "pmid":     str(pmid),
                "title":    r.get("title", "") or "",
                "abstract": r.get("abstractText", "") or "",
                "year":     str(r.get("pubYear", "") or ""),
            })

        next_cursor = data.get("nextCursorMark")
        if not next_cursor or next_cursor == cursor or len(results) < page_size:
            break  # no more pages
        cursor = next_cursor
        time.sleep(delay)

    return papers[:max_results]


# ── Query clause builders (Europe PMC field syntax) ────────────────────────

def alias_clause(organ: str, organ_aliases: dict, organ_mesh: dict) -> str:
    """
    Europe PMC equivalent of the PubMed pipeline's _alias_clause(): MeSH
    term OR each alias, searched across title+abstract.
    organ_mesh values look like 'Liver[MeSH Terms]' (PubMed field-tag
    syntax) — only the bare term before '[' is reused here, since Europe
    PMC's MESH: field takes a plain term instead.
    """
    mesh_term   = organ_mesh.get(organ, organ).split("[")[0].strip()
    aliases     = organ_aliases.get(organ, [organ])
    alias_terms = " OR ".join(f'TITLE:"{a}" OR ABSTRACT:"{a}"' for a in aliases)
    return f'(MESH:"{mesh_term}" OR {alias_terms})'


def keyword_clause(keywords: list[str]) -> str:
    """
    Europe PMC equivalent of the PubMed pipeline's _keyword_clause(): at
    least one keyword must appear in title or abstract.
    """
    parts = [f'TITLE:"{kw}" OR ABSTRACT:"{kw}"' for kw in keywords]
    return "(" + " OR ".join(parts) + ")"
