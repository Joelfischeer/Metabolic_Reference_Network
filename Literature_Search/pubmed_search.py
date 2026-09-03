"""
Shared PubMed/NCBI E-utilities helpers and biomedical vocabulary used by all
three organ cross-talk pipelines (Edge_cosine_met_reference_network,
Edge_cosine_general_reference_network, reference_network_only_metabolic).

Each pipeline builds and issues its own PubMed queries (see each pipeline's
run_network.py / run_metabolic_lit_search.py); this module provides the
lower-level pieces they share: organ name aliases and MeSH terms
(ORGAN_ALIASES / ORGAN_MESH), abstract fetching (fetch_abstracts), key-player
extraction against curated hormone/metabolite/protein vocabularies
(extract_key_players), and synonym merging (_merge_synonyms).
"""

import json
import time
import re
from pathlib import Path
from datetime import datetime

import requests

# ---------------------------------------------------------------------------
# Organ name aliases — maps canonical name → list of PubMed search terms
# ---------------------------------------------------------------------------

ORGAN_ALIASES: dict[str, list[str]] = {
    "Adrenal Glands":   ["adrenal gland", "adrenal cortex", "adrenal medulla",
                         "adrenocortical"],
    "Bone Marrow":      ["bone marrow", "hematopoietic", "haematopoietic",
                         "hematopoiesis", "bone marrow niche"],
    "Brain":            ["brain", "cerebral", "hypothalamus", "hypothalamic", "pituitary", "pituitary gland",
                         "arcuate nucleus", "paraventricular nucleus", "nucleus tractus solitarius",
                         "area postrema", "hippocampus", "hippocampal", "amygdala", "brainstem", "brain stem",
                         "cerebellum", "cerebellar", "thalamus", "thalamic", "basal ganglia",
                         "prefrontal cortex", "cerebral cortex", "medulla oblongata", "pons",
                         "forebrain", "midbrain", "hindbrain", "central nervous system", "CNS",
                         "neuronal", "neuroendocrine"],
    "Colon":            ["colon", "large intestine", "colorectal",
                         "colonic", "gut microbiota", "microbiome"],
    "Heart":            ["heart", "cardiac", "myocardial", "cardiomyocyte",
                         "myocardium", "cardiac muscle"],
    "Kidney":           ["kidney", "renal", "nephron", "glomerular", "kidneys",
                         "tubular"],
    "Liver":            ["liver", "hepatic", "hepatocyte", "hepatocellular"],
    "Lung":             ["lung", "pulmonary", "alveolar", "bronchial"],
    "Muscle":           ["skeletal muscle", "myocyte"],
    "Pancreas":         ["pancreas", "pancreatic"],
    "Small Intestine":  ["small intestine", "duodenum", "jejunum", "ileum",
                         "intestinal", "gut", "enterocyte", "Peyer's patches"],
    "Spleen":           ["spleen", "splenic", "splenomegaly",
                         "lymphoid organ", "immune organ"],
    "Thyroid":          ["thyroid", "thyroid gland"],
    "WAT":              ["adipose tissue", "white adipose tissue", "white fat tissue",
                         "adipocyte", "adipogenesis", "visceral fat",
                         "subcutaneous fat", "WAT"],
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

# ---------------------------------------------------------------------------
# Synonym / alias groups
# Each entry maps a canonical display name → all vocabulary terms it covers.
# Terms not listed here keep their raw vocabulary form as the display name.
# ---------------------------------------------------------------------------

HORMONE_SYNONYMS: dict[str, set[str]] = {
    "epinephrine/adrenaline":        {"adrenaline", "epinephrine"},
    "norepinephrine/noradrenaline":  {"noradrenaline", "norepinephrine"},
    "T3 (triiodothyronine)":         {"t3", "triiodothyronine"},
    "T4 (thyroxine)":                {"t4", "thyroxine"},
    "TSH":                           {"tsh", "thyrotropin"},
    "thyroid hormone":               {"thyroid hormone", "thyroid hormones"},
    "IGF-1":                         {"igf-1", "igf1"},
    "vasopressin/ADH":               {"vasopressin", "adh"},
    "GLP-1":                         {"glp-1", "glp1"},
    "IL-6":                          {"il-6", "il6", "interleukin-6"},
    "TNF-α":                         {"tnf-alpha", "tnf", "tumor necrosis factor"},
    "IL-1β":                         {"il-1beta", "il1b", "interleukin-1"},
    "TGF-β":                         {"tgf-beta", "tgf-b"},
    "erythropoietin (EPO)":          {"erythropoietin", "epo"},
    "DHEA":                          {"dhea", "dhea-s"},
}

METABOLITE_SYNONYMS: dict[str, set[str]] = {
    "fatty acids":       {"fatty acid", "fatty acids"},
    "free fatty acids":  {"free fatty acid", "free fatty acids", "ffa", "nefa"},
    "triglycerides":     {"triglyceride", "triglycerides"},
    "ketone bodies":     {"ketone body", "ketone bodies"},
    "amino acids":       {"amino acid", "amino acids"},
    "bile acids":        {"bile acid", "bile acids"},
}

PROTEIN_SYNONYMS: dict[str, set[str]] = {
    "mTOR":                  {"mtor", "mtorc1", "mtorc2"},
    "IRS-1":                 {"irs-1", "irs1"},
    "IRS-2":                 {"irs-2", "irs2"},
    "PPARα":                 {"ppar-alpha", "ppara"},
    "PPARγ":                 {"ppar-gamma", "pparg"},
    "PGC-1α":                {"pgc-1alpha", "pgc1a", "pgc-1a", "pgc-1"},
    "SREBP":                 {"srebp", "srebp-1c"},
    "HIF-1α":                {"hif-1alpha", "hif1a"},
    "NF-κB":                 {"nf-kb", "nfkb"},
    "LPL":                   {"lpl", "lipoprotein lipase"},
    "HSL":                   {"hsl", "hormone-sensitive lipase"},
    "ATGL":                  {"atgl", "adipose triglyceride lipase"},
    "fatty acid synthase":   {"fas", "fatty acid synthase"},
    "hexokinase":            {"hexokinase", "hk1", "hk2"},
    "FATP":                  {"fatp", "fatp1", "fatp4"},
    "GLP-1 receptor":        {"glp-1 receptor", "glp1r"},
    "FOXO":                  {"foxo1", "foxo"},
    "COX-2":                 {"cox-2", "cox2"},
    "iNOS":                  {"nos2", "inos"},
    "ERK":                   {"erk1", "erk2"},
}


def _merge_synonyms(counts: dict[str, int],
                    synonym_groups: dict[str, set[str]]) -> dict[str, int]:
    """Sum counts for all synonym variants into their canonical display name."""
    term_to_canon: dict[str, str] = {
        term: canon
        for canon, terms in synonym_groups.items()
        for term in terms
    }
    merged: dict[str, int] = {}
    for term, n in counts.items():
        canon = term_to_canon.get(term, term)
        merged[canon] = merged.get(canon, 0) + n
    return merged


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
# Key player extraction
# ---------------------------------------------------------------------------

def extract_key_players(papers: list[dict]) -> dict:
    """
    Scan titles, abstracts, and MeSH keywords for known biomedical terms.

    Each term is counted at most once per paper regardless of how many times
    it appears in that paper's text.

    Returns a dict with two parallel structures per category:
      - "<cat>":        list of names sorted by descending paper count
      - "<cat>_counts": dict of {name: n_papers_containing_term}
    """
    paper_texts = [
        (
            p.get("title", "") + " " +
            p.get("abstract", "") + " " +
            " ".join(p.get("keywords", []))
        ).lower()
        for p in papers
    ]

    def find_terms(vocab: set[str],
                   synonyms: dict[str, set[str]]) -> tuple[list[str], dict[str, int]]:
        # Group raw vocab terms by canonical name FIRST, then count each
        # paper at most once per canonical term (across all its synonym
        # variants) — not once per raw variant. Counting per-variant and
        # summing afterward (the previous approach) double-counts any paper
        # that happens to use two synonyms of the same term (e.g. "il-6" and
        # "interleukin-6" in the same abstract), since the per-paper
        # membership information is already lost by the time raw counts are
        # merged. See threshold_utils.py-style discussion — synonym merging
        # must happen at match time, not after aggregation.
        term_to_canon = {
            term: canon for canon, terms in synonyms.items() for term in terms
        }
        canon_to_terms: dict[str, list[str]] = {}
        for term in vocab:
            canon_to_terms.setdefault(term_to_canon.get(term, term), []).append(term)

        counts: dict[str, int] = {}
        for canon, terms in canon_to_terms.items():
            patterns = [r"\b" + re.escape(t) + r"\b" for t in terms]
            n = sum(
                1 for txt in paper_texts
                if any(re.search(p, txt) for p in patterns)
            )
            if n > 0:
                counts[canon] = n
        ranked = [t for t, _ in sorted(counts.items(), key=lambda x: -x[1])]
        return ranked, counts

    metabolites, metabolites_counts = find_terms(METABOLITES, METABOLITE_SYNONYMS)
    hormones,    hormones_counts    = find_terms(HORMONES,    HORMONE_SYNONYMS)
    proteins,    proteins_counts    = find_terms(PROTEINS,    PROTEIN_SYNONYMS)

    return {
        "metabolites":        metabolites,
        "metabolites_counts": metabolites_counts,
        "hormones":           hormones,
        "hormones_counts":    hormones_counts,
        "proteins":           proteins,
        "proteins_counts":    proteins_counts,
    }


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


# ---------------------------------------------------------------------------
# Excel export of curated vocabulary lists
# ---------------------------------------------------------------------------

def export_vocabulary_to_excel(output_path: "str | Path" = None) -> Path:
    """
    Write the three curated key-player vocabulary lists to an Excel workbook.

    Each category gets its own sheet:
      - Hormones
      - Metabolites
      - Proteins

    Terms are sorted alphabetically within each sheet.
    """
    import pandas as pd

    if output_path is None:
        output_path = Path(__file__).parent.parent / "metabolic_data" / "key_player_vocabulary.xlsx"
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    sheets = {
        "Hormones":    sorted(HORMONES,    key=str.lower),
        "Metabolites": sorted(METABOLITES, key=str.lower),
        "Proteins":    sorted(PROTEINS,    key=str.lower),
    }

    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        for sheet_name, terms in sheets.items():
            df = pd.DataFrame({"Term": terms})
            df.index = range(1, len(df) + 1)
            df.to_excel(writer, sheet_name=sheet_name, index_label="No.")

            # Basic column width
            ws = writer.sheets[sheet_name]
            ws.column_dimensions["A"].width = 6
            ws.column_dimensions["B"].width = max(len(t) for t in terms) + 4

    print(f"[ok] Vocabulary exported: {output_path}  "
          f"({len(HORMONES)} hormones, {len(METABOLITES)} metabolites, {len(PROTEINS)} proteins)")
    return output_path

