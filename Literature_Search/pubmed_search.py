"""
PubMed literature search for organ-organ metabolic connections.

Uses NCBI E-utilities (no API key required, but email recommended).
Results are saved as JSON per edge and can be resumed incrementally.

Search strategy (cascade — moves to next level if too few results):
  1. Both organ names in Title with metabolic context keywords
  2. Both names anywhere in Title/Abstract
  3. MeSH terms for both organs
  4. All organ aliases, no metabolic filter — just co-occurrence
  5. PMC full-text search (broader coverage)
"""

import json
import time
import re
from pathlib import Path
from datetime import datetime, timedelta

import requests

# ---------------------------------------------------------------------------
# Organ name aliases — maps canonical name → list of PubMed search terms
# ---------------------------------------------------------------------------

ORGAN_ALIASES: dict[str, list[str]] = {
    "Adrenal Glands":   ["adrenal gland", "adrenal cortex", "adrenal medulla",
                         "adrenocortical", "HPA axis", "hypothalamic-pituitary-adrenal"],
    "Bone Marrow":      ["bone marrow", "hematopoietic", "haematopoietic",
                         "hematopoiesis", "bone marrow niche"],
    "Brain":            ["brain", "cerebral", "hypothalamus", "hypothalamic",
                         "central nervous system", "CNS", "neuronal"],
    "Colon":            ["colon", "large intestine", "colorectal",
                         "colonic", "gut microbiota", "microbiome"],
    "Heart":            ["heart", "cardiac", "myocardial", "cardiomyocyte",
                         "left ventricle", "myocardium"],
    "Kidney":           ["kidney", "renal", "nephron", "glomerular",
                         "tubular", "nephrology"],
    "Liver":            ["liver", "hepatic", "hepatocyte", "hepatocellular",
                         "nonalcoholic fatty liver", "NAFLD", "NASH"],
    "Lung":             ["lung", "pulmonary", "alveolar", "bronchial",
                         "respiratory", "pneumocyte"],
    "Muscle":           ["skeletal muscle", "muscle", "myocyte", "myofiber",
                         "sarcopenia", "muscular"],
    "Pancreas":         ["pancreas", "pancreatic", "islet", "beta cell",
                         "beta-cell", "insulin secretion", "exocrine pancreas"],
    "Small Intestine":  ["small intestine", "duodenum", "jejunum", "ileum",
                         "intestinal", "gut", "enterocyte", "Peyer's patches",
                         "gut-brain axis", "enteroendocrine"],
    "Spleen":           ["spleen", "splenic", "splenomegaly",
                         "lymphoid organ", "immune organ"],
    "Thyroid":          ["thyroid", "thyroid gland", "thyroid hormone",
                         "T3", "T4", "TSH", "hypothyroidism", "hyperthyroidism"],
    "WAT":              ["adipose tissue", "white adipose tissue", "fat tissue",
                         "adipocyte", "adipogenesis", "visceral fat",
                         "subcutaneous fat", "WAT", "lipogenesis"],
}

# MeSH terms for organ MeSH-based queries
ORGAN_MESH: dict[str, str] = {
    "Adrenal Glands":  "Adrenal Glands[MeSH Terms]",
    "Bone Marrow":     "Bone Marrow[MeSH Terms]",
    "Brain":           "Brain[MeSH Terms]",
    "Colon":           "Colon[MeSH Terms]",
    "Heart":           "Heart[MeSH Terms]",
    "Kidney":          "Kidney[MeSH Terms]",
    "Liver":           "Liver[MeSH Terms]",
    "Lung":            "Lung[MeSH Terms]",
    "Muscle":          "Muscle, Skeletal[MeSH Terms]",
    "Pancreas":        "Pancreas[MeSH Terms]",
    "Small Intestine": "Intestine, Small[MeSH Terms]",
    "Spleen":          "Spleen[MeSH Terms]",
    "Thyroid":         "Thyroid Gland[MeSH Terms]",
    "WAT":             "Adipose Tissue[MeSH Terms]",
}

METABOLIC_CONTEXT = (
    "(metabolism[Title/Abstract] OR metabolic[Title/Abstract] "
    "OR crosstalk[Title/Abstract] OR signaling[Title/Abstract] "
    "OR communication[Title/Abstract] OR axis[Title/Abstract] "
    "OR interaction[Title/Abstract] OR regulation[Title/Abstract] "
    "OR hormone[Title/Abstract] OR glucose[Title/Abstract] "
    "OR insulin[Title/Abstract] OR lipid[Title/Abstract] "
    "OR energy[Title/Abstract])"
)

# ---------------------------------------------------------------------------
# Curated biomedical vocabulary for key-player extraction
# ---------------------------------------------------------------------------

# Small molecules, energy substrates, lipids, organic acids — NOT hormones/proteins
METABOLITES = {
    # Carbohydrates & glycolytic intermediates
    "glucose", "fructose", "galactose", "glycogen", "lactate", "pyruvate",
    "phosphoenolpyruvate", "glucose-6-phosphate", "fructose-1,6-bisphosphate",
    # TCA cycle intermediates
    "citrate", "isocitrate", "alpha-ketoglutarate", "succinate", "fumarate",
    "malate", "oxaloacetate", "acetyl-coa",
    # Fatty acids & lipids
    "fatty acid", "fatty acids", "free fatty acid", "free fatty acids",
    "triglyceride", "triglycerides", "diacylglycerol", "monoacylglycerol",
    "cholesterol", "phospholipid", "sphingolipid", "ceramide",
    "palmitate", "oleate", "linoleate", "stearate", "arachidonate",
    "ffa", "nefa",
    # Ketone bodies
    "ketone body", "ketone bodies", "beta-hydroxybutyrate", "acetoacetate",
    # Short-chain fatty acids
    "acetate", "propionate", "butyrate", "valerate",
    # Glycerol & alcohol sugars
    "glycerol", "glycerol-3-phosphate",
    # Energy currency
    "atp", "adp", "amp", "nadh", "nadph", "nad+", "fadh2",
    # Amino acids & nitrogen metabolites
    "amino acid", "amino acids", "glutamine", "glutamate", "alanine",
    "leucine", "isoleucine", "valine", "branched-chain amino acid",
    "arginine", "serine", "glycine", "cysteine", "methionine",
    "urea", "ammonia", "creatinine", "creatine", "uric acid",
    # Bile acids (small molecules, not proteins/hormones)
    "bile acid", "bile acids", "cholic acid", "deoxycholic acid",
    "chenodeoxycholic acid", "lithocholic acid", "ursodeoxycholic acid",
    # Other organic acids & signaling lipids
    "bilirubin", "biliverdin",
    "prostaglandin", "thromboxane", "leukotriene", "eicosanoid",
    "lysophosphatidylcholine", "sphingosine-1-phosphate",
    # Lipoproteins (lipid-carrying particles, not hormones)
    "hdl", "ldl", "vldl", "lipoprotein", "chylomicron",
}

# Secreted signalling molecules (peptide/steroid hormones, cytokines, growth factors)
HORMONES = {
    # Pancreatic hormones
    "insulin", "glucagon", "somatostatin", "amylin",
    # Adrenal hormones
    "cortisol", "corticosterone", "aldosterone", "dhea",
    "adrenaline", "epinephrine", "noradrenaline", "norepinephrine",
    "glucocorticoid", "mineralocorticoid",
    # Thyroid hormones
    "t3", "triiodothyronine", "t4", "thyroxine", "tsh", "thyrotropin",
    "thyroid hormone", "thyroid hormones",
    # Pituitary / hypothalamic hormones
    "growth hormone", "igf-1", "igf1", "igf-2",
    "acth", "crh", "gnrh", "lh", "fsh",
    "prolactin", "oxytocin", "vasopressin", "adh",
    # Adipokines
    "leptin", "adiponectin", "resistin", "visfatin", "chemerin",
    "omentin", "apelin",
    # Gut hormones
    "ghrelin", "glp-1", "glp1", "gip", "pyy", "cck",
    "secretin", "motilin", "gastrin", "neurotensin",
    # Hepatokines
    "fgf21", "fgf19", "fetuin-a", "hepassocin", "selenoprotein p",
    # Myokines
    "irisin", "meteorin-like", "il-6", "il6", "interleukin-6",
    "fndc5",
    # Sex steroids
    "testosterone", "estrogen", "estradiol", "estrone", "progesterone",
    "androstenedione", "dhea-s",
    # Cardiovascular / renal hormones
    "angiotensin", "angiotensin ii", "renin", "erythropoietin", "epo",
    "natriuretic peptide", "bnp", "anp", "cnp",
    # Bone-derived
    "osteocalcin", "fgf23", "klotho",
    # Cytokines with systemic metabolic roles
    "tnf-alpha", "tnf", "tumor necrosis factor",
    "il-1beta", "il1b", "interleukin-1",
    "il-10", "il-4", "tgf-beta", "tgf-b",
    # Vitamin D (acts as hormone)
    "calcitriol", "vitamin d",
}

# Enzymes, transporters, receptors, and intracellular signalling proteins
PROTEINS = {
    # Glucose transporters
    "glut1", "glut2", "glut3", "glut4", "glut5",
    "sglt1", "sglt2",
    # Insulin / nutrient signalling kinases
    "ampk", "mtor", "mtorc1", "mtorc2",
    "pi3k", "akt", "pdk1",
    "irs-1", "irs-2", "irs1", "irs2",
    "mapk", "erk1", "erk2", "jnk", "p38",
    # Transcription factors
    "ppar-alpha", "ppar-gamma", "ppar-delta", "ppara", "pparg",
    "pgc-1alpha", "pgc1a", "pgc-1a", "pgc-1beta",
    "srebp", "srebp-1c", "chrebp", "foxo1", "foxo",
    "hif-1alpha", "hif1a", "nf-kb", "nfkb",
    "lxr", "fxr", "rxr", "tgr5",
    # Lipid metabolism enzymes & transporters
    "cd36", "fatp", "fatp1", "fatp4", "fabp",
    "cpt1", "cpt2", "acc", "fas", "fatty acid synthase",
    "lpl", "lipoprotein lipase",
    "hsl", "hormone-sensitive lipase",
    "atgl", "adipose triglyceride lipase",
    "dgat1", "dgat2",
    # Gluconeogenesis / glycolysis enzymes
    "pepck", "g6pase", "glucokinase", "hexokinase",
    "hk1", "hk2", "pfk", "aldolase", "gapdh", "pkm2",
    "pdk", "pdk4", "pdhc",
    # Mitochondrial proteins
    "ucp1", "ucp2", "ucp3",
    "sirt1", "sirt3", "sirt4",
    "pgc-1", "tfam",
    # Receptor proteins
    "insulin receptor", "glucagon receptor",
    "leptin receptor", "adiponectin receptor",
    "glp-1 receptor", "glp1r",
    # Angiopoietin-like proteins
    "angptl4", "fiaf",
    # Bile acid synthesis
    "cyp7a1", "cyp27a1",
    # Inflammation-related enzymes
    "cox-2", "cox2", "nos2", "inos",
}

NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


# ---------------------------------------------------------------------------
# HTTP helpers
# ---------------------------------------------------------------------------

def _ncbi_get(endpoint: str, params: dict, retries: int = 3) -> "requests.Response | None":
    params.setdefault("tool", "MetabolicReferenceNetwork")
    params.setdefault("email", "research@metabolic-network.org")
    for attempt in range(retries):
        try:
            resp = requests.get(NCBI_BASE + endpoint, params=params, timeout=25)
            resp.raise_for_status()
            return resp
        except requests.RequestException as e:
            wait = 2 ** attempt
            print(f"    [!] HTTP error (attempt {attempt+1}/{retries}): {e} — retrying in {wait}s")
            time.sleep(wait)
    return None


def search_pubmed_raw(query: str, max_results: int, years_back: int,
                      db: str = "pubmed") -> list[str]:
    """Return PMIDs from a single PubMed/PMC query."""
    min_date = (datetime.now() - timedelta(days=365 * years_back)).strftime("%Y/%m/%d")
    params = {
        "db": db,
        "term": query,
        "retmax": max_results,
        "retmode": "json",
        "mindate": min_date,
        "datetype": "pdat",
    }
    resp = _ncbi_get("esearch.fcgi", params)
    if resp is None:
        return []
    data = resp.json()
    ids = data.get("esearchresult", {}).get("idlist", [])

    # PMC returns PMC IDs — convert to PubMed IDs via elink
    if db == "pmc" and ids:
        ids = _pmc_to_pubmed(ids)
    return ids


def _pmc_to_pubmed(pmc_ids: list[str]) -> list[str]:
    """Convert PMC IDs to PubMed IDs via elink."""
    params = {
        "dbfrom": "pmc",
        "db": "pubmed",
        "id": ",".join(pmc_ids),
        "retmode": "json",
        "cmd": "neighbor",
    }
    resp = _ncbi_get("elink.fcgi", params)
    if resp is None:
        return []
    try:
        data = resp.json()
        linksets = data.get("linksets", [])
        pmids = []
        for ls in linksets:
            for lsdb in ls.get("linksetdbs", []):
                if lsdb.get("dbto") == "pubmed":
                    pmids.extend(lsdb.get("links", []))
        return list(dict.fromkeys(pmids))  # deduplicate, preserve order
    except Exception:
        return []


def fetch_abstracts(pmids: list[str], batch_size: int = 200, delay: float = 0.4) -> list[dict]:
    """
    Fetch title + abstract + keywords for a list of PMIDs.

    PMIDs are fetched in batches of `batch_size` (PubMed efetch is unreliable
    above ~200 IDs per request). Results from all batches are combined.
    """
    if not pmids:
        return []

    def clean(text: str) -> str:
        text = re.sub(r"<[^>]+>", " ", text)
        return re.sub(r"\s+", " ", text).strip()

    def parse_xml_batch(xml: str) -> list[dict]:
        papers = []
        article_blocks = re.findall(r"<PubmedArticle>(.*?)</PubmedArticle>", xml, re.DOTALL)
        for block in article_blocks:
            pmid_m  = re.search(r"<PMID[^>]*>(\d+)</PMID>", block)
            title_m = re.search(r"<ArticleTitle>(.*?)</ArticleTitle>", block, re.DOTALL)
            year_m  = re.search(r"<PubDate>.*?<Year>(\d{4})</Year>", block, re.DOTALL)
            doi_m   = re.search(r'<ArticleId IdType="doi">(.*?)</ArticleId>', block)

            abstract_parts = re.findall(r"<AbstractText[^>]*>(.*?)</AbstractText>", block, re.DOTALL)
            abstract = " ".join(clean(p) for p in abstract_parts)

            kw_matches = re.findall(r"<DescriptorName[^>]*>(.*?)</DescriptorName>", block, re.DOTALL)
            keywords = [clean(k) for k in kw_matches]

            papers.append({
                "pmid":     pmid_m.group(1) if pmid_m else "",
                "title":    clean(title_m.group(1)) if title_m else "",
                "abstract": abstract,
                "keywords": keywords,
                "year":     year_m.group(1) if year_m else "",
                "doi":      doi_m.group(1).strip() if doi_m else "",
            })
        return papers

    all_papers: list[dict] = []
    batches = [pmids[i:i + batch_size] for i in range(0, len(pmids), batch_size)]

    for i, batch in enumerate(batches):
        if i > 0:
            time.sleep(delay)
        params = {
            "db":      "pubmed",
            "id":      ",".join(batch),
            "rettype": "abstract",
            "retmode": "xml",
        }
        resp = _ncbi_get("efetch.fcgi", params)
        if resp is not None:
            all_papers.extend(parse_xml_batch(resp.text))

    return all_papers


# ---------------------------------------------------------------------------
# Query builders — cascade from strict to broad
# ---------------------------------------------------------------------------

def _alias_clause(organ: str, fields: str = "Title/Abstract") -> str:
    """Build an OR clause for all aliases of an organ."""
    aliases = ORGAN_ALIASES.get(organ, [organ])
    terms = [f'"{a}"[{fields}]' for a in aliases]
    return "(" + " OR ".join(terms) + ")"


def build_queries(organ1: str, organ2: str) -> list[tuple[str, str]]:
    """
    Return a list of (label, query_string) in order of decreasing strictness.
    The search cascade tries them in order and stops when enough results are found.
    """
    a1 = _alias_clause(organ1, "Title")
    a2 = _alias_clause(organ2, "Title")
    ab1 = _alias_clause(organ1, "Title/Abstract")
    ab2 = _alias_clause(organ2, "Title/Abstract")
    mesh1 = ORGAN_MESH.get(organ1, f'"{organ1}"[MeSH Terms]')
    mesh2 = ORGAN_MESH.get(organ2, f'"{organ2}"[MeSH Terms]')

    return [
        # Strategy 1: both organs in Title + metabolic context
        ("title+context",
         f"{a1} AND {a2} AND {METABOLIC_CONTEXT}"),

        # Strategy 2: both organs in Title/Abstract + metabolic context
        ("abstract+context",
         f"{ab1} AND {ab2} AND {METABOLIC_CONTEXT}"),

        # Strategy 3: MeSH terms + metabolic context
        ("mesh+context",
         f"({mesh1} OR {ab1}) AND ({mesh2} OR {ab2}) AND {METABOLIC_CONTEXT}"),

        # Strategy 4: MeSH terms alone (broad)
        ("mesh+aliases",
         f"({mesh1} OR {ab1}) AND ({mesh2} OR {ab2})"),

        # Strategy 5: PMC full-text search (uses pmc db in caller)
        ("pmc-fulltext",
         f"{ab1} AND {ab2} AND {METABOLIC_CONTEXT}"),
    ]


# ---------------------------------------------------------------------------
# Main search function with cascade
# ---------------------------------------------------------------------------

def search_with_cascade(
    organ1: str,
    organ2: str,
    max_results: int = 25,
    years_back: int = 10,
    min_papers: int = 5,
    delay: float = 0.4,
) -> tuple[list[str], str, str]:
    """
    Try query strategies in order until `min_papers` PMIDs are found.

    Returns:
        (pmids, strategy_used, query_used)
    """
    queries = build_queries(organ1, organ2)
    best_pmids: list[str] = []
    best_strategy = ""
    best_query = ""

    for label, query in queries:
        db = "pmc" if label == "pmc-fulltext" else "pubmed"
        print(f"    [?] Strategy '{label}' (db={db})…", end=" ", flush=True)
        time.sleep(delay)

        pmids = search_pubmed_raw(query, max_results=max_results,
                                  years_back=years_back, db=db)
        print(f"{len(pmids)} results")

        # Keep the best result so far (most papers)
        if len(pmids) > len(best_pmids):
            best_pmids = pmids
            best_strategy = label
            best_query = query

        if len(pmids) >= min_papers:
            break  # good enough — stop cascade

    return best_pmids, best_strategy, best_query


# ---------------------------------------------------------------------------
# Key player extraction
# ---------------------------------------------------------------------------

def extract_key_players(papers: list[dict]) -> dict:
    """
    Scan titles, abstracts, and MeSH keywords for known biomedical terms.

    Returns a dict with two parallel structures per category:
      - "<cat>":        list of names sorted by descending mention count
      - "<cat>_counts": dict of {name: mention_count}
    """
    combined_text = " ".join(
        (
            p.get("title", "") + " " +
            p.get("abstract", "") + " " +
            " ".join(p.get("keywords", []))
        ).lower()
        for p in papers
    )

    def find_terms(vocab: set[str]) -> tuple[list[str], dict[str, int]]:
        counts: dict[str, int] = {}
        for term in vocab:
            pattern = r"\b" + re.escape(term) + r"\b"
            n = len(re.findall(pattern, combined_text))
            if n > 0:
                counts[term] = n
        ranked = [t for t, _ in sorted(counts.items(), key=lambda x: -x[1])]
        return ranked, counts

    metabolites, metabolites_counts = find_terms(METABOLITES)
    hormones,    hormones_counts    = find_terms(HORMONES)
    proteins,    proteins_counts    = find_terms(PROTEINS)

    return {
        "metabolites":        metabolites,
        "metabolites_counts": metabolites_counts,
        "hormones":           hormones,
        "hormones_counts":    hormones_counts,
        "proteins":           proteins,
        "proteins_counts":    proteins_counts,
    }


def _infer_connection_type(key_players: dict[str, list[str]]) -> str:
    hormones    = set(key_players.get("hormones", []))
    metabolites = set(key_players.get("metabolites", []))
    cytokines   = {"il-6", "il6", "tnf-alpha", "tnf", "il-1beta", "il1b",
                   "tgf-beta", "il-10"}

    if hormones & {"insulin", "glucagon", "cortisol", "corticosterone",
                   "leptin", "adiponectin", "tsh", "t3", "t4",
                   "thyroid hormone", "thyroid hormones",
                   "ghrelin", "glp-1", "glp1", "fgf21", "igf-1", "igf1",
                   "growth hormone"}:
        if metabolites & {"glucose", "fatty acid", "fatty acids", "lactate",
                          "triglyceride", "triglycerides", "ketone body",
                          "ketone bodies", "beta-hydroxybutyrate"}:
            return "Hormonal-Metabolic"
        return "Hormonal"
    if hormones & cytokines:
        return "Inflammatory/Immune"
    if metabolites & {"glucose", "fatty acid", "fatty acids", "lactate",
                      "ketone body", "ketone bodies", "triglyceride",
                      "bile acid", "bile acids"}:
        return "Metabolic Substrate Exchange"
    if hormones:
        return "Hormonal"
    if metabolites:
        return "Metabolic"
    return "Unknown"


# ---------------------------------------------------------------------------
# Public run function
# ---------------------------------------------------------------------------

def run_literature_search(
    organ_pairs: list[tuple[str, str]],
    output_path: "str | Path",
    max_results_per_pair: int = 25,
    years_back: int = 10,
    min_papers: int = 5,
    delay: float = 0.4,
    resume: bool = True,
    force_research_empty: bool = False,
) -> dict:
    """
    Run PubMed searches for all organ pairs with a multi-strategy cascade.

    Parameters
    ----------
    organ_pairs : list of (organ1, organ2)
    output_path : JSON file for storing results (auto-saved after each pair)
    max_results_per_pair : max PMIDs to retrieve per strategy per pair
    years_back : how many years back to search
    min_papers : cascade stops when this many papers are found
    delay : seconds between HTTP requests
    resume : skip already-searched pairs (cached in output_path)
    force_research_empty : re-search pairs that previously returned 0 papers
    """
    output_path = Path(output_path)
    results: dict = {}

    if resume and output_path.exists():
        with open(output_path, encoding="utf-8") as f:
            results = json.load(f)
        n_cached = len(results)
        n_empty  = sum(1 for v in results.values() if v.get("n_papers_found", 0) == 0)
        print(f"[i] Resuming: {n_cached} cached ({n_empty} empty).")
        if force_research_empty:
            # Delete empty entries so they get re-searched
            empty_keys = [k for k, v in results.items() if v.get("n_papers_found", 0) == 0]
            for k in empty_keys:
                del results[k]
            print(f"[i] Cleared {len(empty_keys)} empty entries for re-search.")

    total = len(organ_pairs)
    for idx, (o1, o2) in enumerate(organ_pairs):
        edge_key = f"{o1}|{o2}"
        sym_key  = f"{o2}|{o1}"

        if edge_key in results or sym_key in results:
            existing = results.get(edge_key) or results.get(sym_key, {})
            n = existing.get("n_papers_found", 0)
            print(f"  [{idx+1}/{total}] Cached ({n} papers): {o1} <-> {o2}")
            continue

        print(f"\n  [{idx+1}/{total}] Searching: {o1} <-> {o2}")
        pmids, strategy, query = search_with_cascade(
            o1, o2,
            max_results=max_results_per_pair,
            years_back=years_back,
            min_papers=min_papers,
            delay=delay,
        )
        time.sleep(delay)

        papers = fetch_abstracts(pmids, delay=delay)
        time.sleep(delay)

        key_players     = extract_key_players(papers)
        connection_type = _infer_connection_type(key_players)

        n = len(papers)
        summary = (
            f"    => {n} papers via strategy '{strategy}'"
            f" | hormones: {len(key_players['hormones'])}"
            f" | metabolites: {len(key_players['metabolites'])}"
            f" | proteins: {len(key_players['proteins'])}"
        )
        print(summary)

        results[edge_key] = {
            "organ1":          o1,
            "organ2":          o2,
            "pubmed_query":    query,
            "strategy_used":   strategy,
            "n_papers_found":  n,
            "papers":          papers,
            "key_players":     key_players,
            "connection_type": connection_type,
            "search_date":     datetime.now().isoformat(),
        }

        with open(output_path, "w", encoding="utf-8") as f:
            json.dump(results, f, indent=2, ensure_ascii=False)

    total_papers = sum(v.get("n_papers_found", 0) for v in results.values())
    empty = sum(1 for v in results.values() if v.get("n_papers_found", 0) == 0)
    print(f"\n[ok] Done. {len(results)} edges | {total_papers} total papers | {empty} still empty.")
    print(f"     Saved to: {output_path}")
    return results


# ---------------------------------------------------------------------------
# Helpers used by other modules
# ---------------------------------------------------------------------------

def load_literature_results(path: "str | Path") -> dict:
    path = Path(path)
    if not path.exists():
        return {}
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def get_edge_literature(results: dict, organ1: str, organ2: str) -> "dict | None":
    return results.get(f"{organ1}|{organ2}") or results.get(f"{organ2}|{organ1}")


def merge_with_edge_metadata(
    edge_metadata: dict[tuple[str, str], str],
    literature_results: dict,
    llm_descriptions: dict | None = None,
) -> dict[tuple[str, str], dict]:
    """Combine curated edge metadata with literature search results."""
    merged = {}

    all_pairs = {
        (o1, o2) if o1 < o2 else (o2, o1)
        for (o1, o2) in edge_metadata
    }

    llm_descriptions = llm_descriptions or {}

    for (o1, o2) in all_pairs:
        text   = edge_metadata.get((o1, o2), "")
        lit    = get_edge_literature(literature_results, o1, o2) or {}
        parsed = _parse_edge_text(text)
        llm_entry = (
            llm_descriptions.get(f"{o1}|{o2}")
            or llm_descriptions.get(f"{o2}|{o1}")
            or {}
        )
        ai_description = llm_entry.get("description", "")

        raw_kp = _filter_key_players(parsed.get("key_players_raw", []))
        lit_kp = lit.get("key_players", {})

        # Categorised lists come purely from literature (ranked by mention count).
        # raw_kp from the curated CSV is kept separate — it is NOT merged into
        # the categorised lists because it would appear in all three categories
        # with no mention-count evidence for this specific edge.
        merged_players = {
            "metabolites": lit_kp.get("metabolites", []),
            "hormones":    lit_kp.get("hormones", []),
            "proteins":    lit_kp.get("proteins", []),
        }
        # Mention counts from PubMed extraction (may be absent in old cached results)
        counts = {
            "metabolites": lit_kp.get("metabolites_counts", {}),
            "hormones":    lit_kp.get("hormones_counts", {}),
            "proteins":    lit_kp.get("proteins_counts", {}),
        }

        merged[(o1, o2)] = {
            "description":           text,
            "connection_type":       parsed.get("connection_type") or lit.get("connection_type", ""),
            "key_players_raw":       raw_kp,
            "key_players_merged":    merged_players,
            "key_players_counts":    counts,
            "notes":                 parsed.get("notes", ""),
            "sources":               parsed.get("sources", []),
            "ai_description":        ai_description,
            "pubmed": {
                "n_papers":     lit.get("n_papers_found", 0),
                "papers":       lit.get("papers", [])[:10],
                "query":        lit.get("pubmed_query", ""),
                "strategy":     lit.get("strategy_used", ""),
            },
        }
        merged[(o2, o1)] = merged[(o1, o2)]

    return merged


def _filter_key_players(items: list[str]) -> list[str]:
    """
    Remove anything that looks like a sentence or free-text note rather than
    an actual molecule/hormone/protein name. Keeps only short, name-like tokens.
    """
    clean = []
    for item in items:
        # Strip parenthetical expansions like "TH (T3, T4)" → keep outer token
        item = re.sub(r"\s*\(.*?\)", "", item).strip()
        if not item:
            continue
        # Drop if it looks like a sentence: ends with '.', '!', '?'
        if item[-1] in ".!?":
            continue
        # Drop if it contains a verb-like word (a sign of free text)
        if re.search(r"\b(is|are|was|were|has|have|the|that|which|this|via|through|by|from|with|into|and|or)\b",
                     item, re.IGNORECASE):
            continue
        # Drop if it's more than 6 words (too long to be a molecule name)
        if len(item.split()) > 6:
            continue
        clean.append(item)
    return clean


def _parse_edge_text(text: str) -> dict:
    result = {"connection_type": "", "key_players_raw": [], "notes": "", "sources": []}
    if not text:
        return result

    m = re.search(r"Type:\s*(.+?)(?:\n|$)", text, re.IGNORECASE)
    if m:
        result["connection_type"] = m.group(1).strip()

    m = re.search(r"Key Players:\s*(.+?)(?:\n|$)", text, re.IGNORECASE)
    if m:
        raw = m.group(1).strip()
        result["key_players_raw"] = _filter_key_players(
            [p.strip() for p in raw.split(",") if p.strip()]
        )

    m = re.search(r"Notes?:\s*(.*?)(?:Sources?:|$)", text, re.IGNORECASE | re.DOTALL)
    if m:
        result["notes"] = m.group(1).strip()

    m = re.search(r"Sources?:\s*(.*?)$", text, re.IGNORECASE | re.DOTALL)
    if m:
        result["sources"] = [
            s.strip().lstrip("- ").strip()
            for s in re.split(r"\n-|\n", m.group(1).strip())
            if s.strip().lstrip("- ").strip()
        ]
    return result


