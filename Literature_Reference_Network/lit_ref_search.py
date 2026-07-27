"""
Literature-derived reference network — organ-pair co-occurrence search.

For every organ-organ pair (91 pairs from 14 organs):
  1. One broad PubMed query retrieves papers mentioning both organs.
  2. Papers are post-filtered: kept only when the two organs co-occur in the
     SAME SENTENCE or as a COMPOUND / HYPHENATED term (e.g. "hepatorenal",
     "liver-gut axis").
  3. Each surviving paper is classified into one or more communication layers
     based on mechanism phrases found in its text.  Papers matching none of the
     five layer phrase sets are placed in the "undefined" layer.
  4. An edge is drawn in a layer when ≥ min_papers papers are classified there.

Results are cached incrementally in
  metabolic_data/lit_ref_results.json

Usage (standalone):
    uv run python Literature_Reference_Network/lit_ref_search.py
    uv run python Literature_Reference_Network/lit_ref_search.py --reset
    uv run python Literature_Reference_Network/lit_ref_search.py --min-papers 3
"""

import json
import re
import sys
import time
import argparse
from pathlib import Path
from datetime import datetime, timedelta
from itertools import combinations

HERE = Path(__file__).parent.parent
sys.path.insert(0, str(HERE))

DEFAULT_OUTPUT      = HERE / "metabolic_data" / "lit_ref_results.json"
DEFAULT_MIN_PAPERS  = 3
DEFAULT_MAX_RESULTS = 200
DEFAULT_YEARS_BACK  = 10
DELAY               = 0.4

# ── Layer metadata ─────────────────────────────────────────────────────────────

LAYER_COLORS = {
    "neural":     "#c084fc",
    "hormonal":   "#fb923c",
    "immune":     "#f87171",
    "metabolic":  "#34d399",
    "mechanical": "#60a5fa",
    "undefined":  "#94a3b8",
}

LAYER_LABELS = {
    "neural":     "Neural",
    "hormonal":   "Hormonal",
    "immune":     "Immune",
    "metabolic":  "Metabolic",
    "mechanical": "Mechanical",
    "undefined":  "Undefined",
}

ALL_LAYERS = list(LAYER_COLORS.keys())

# ── Compound organ stems ────────────────────────────────────────────────────────
# Used to detect merged compound terms like "hepatorenal", "cardiorenal",
# "entero-hepatic", "neuro-muscular", etc.

ORGAN_COMPOUND_STEMS: dict[str, list[str]] = {
    "Liver":          ["hepato", "hepat"],
    "Kidney":         ["nephro", "renal"],
    "Brain":          ["neuro", "cerebro", "encephalo", "hypothalamo"],
    "Heart":          ["cardio", "cardiac"],
    "Lung":           ["pulmo", "pneumo", "broncho"],
    "Pancreas":       ["pancreato"],
    "Colon":          ["colo", "colorect"],
    "Small Intestine":["entero", "intestin"],
    "Muscle":         ["myo", "musculo"],
    "WAT":            ["adipo", "lipid"],
    "Spleen":         ["spleno", "splen"],
    "Bone Marrow":    ["myelo"],
    "Adrenal Glands": ["adreno"],
    "Thyroid":        ["thyro"],
}

# ── Mechanism phrases for layer classification ─────────────────────────────────
# At least one phrase must appear in a paper's title+abstract for it to be
# classified into that layer.

LAYER_REQUIRED_PHRASES: dict[str, list[str]] = {
    "neural": [
        "autonomic innervation", "sympathetic innervation",
        "parasympathetic innervation", "vagal regulation",
        "vagal efferent", "vagal afferent", "neuroendocrine axis",
        "neural regulation", "autonomic regulation", "nerve-mediated",
        "neural crosstalk", "neuronal regulation",
    ],
    "hormonal": [
        "endocrine axis", "hormonal axis", "circulating hormone",
        "endocrine crosstalk", "hormonal crosstalk", "hormonal regulation",
        "endocrine regulation", "blood-borne signal",
        "hormonal communication", "endocrine communication",
        "hormonal signaling", "endocrine signaling",
    ],
    "immune": [
        "immune crosstalk", "inflammatory crosstalk", "cytokine-mediated",
        "immune-mediated", "inflammatory mediator", "immune regulation",
        "cytokine signaling", "immunomodulatory", "immune communication",
        "cytokine crosstalk", "immune axis",
    ],
    "metabolic": [
        "metabolic crosstalk", "inter-organ", "organ crosstalk",
        "organ communication", "substrate exchange", "metabolic communication",
        "metabolic axis", "metabolic interplay", "metabolic interaction",
        "metabolic relay",
    ],
    "mechanical": [
        "hemodynamic coupling", "blood flow-mediated", "flow-mediated",
        "pressure sensing", "stretch sensing", "pressure-mediated",
        "baroreceptor", "mechanotransduction", "mechanosensing",
        "vascular coupling", "mechanical coupling", "hemodynamic regulation",
    ],
}

# ── Layer-specific key-player vocabularies ─────────────────────────────────────

LAYER_VOCAB: dict[str, set] = {
    "metabolic": {
        "glucose", "fructose", "galactose", "glycogen", "lactate", "pyruvate",
        "acetyl-coa", "glucose-6-phosphate",
        "citrate", "succinate", "fumarate", "malate", "oxaloacetate",
        "alpha-ketoglutarate",
        "fatty acid", "free fatty acid", "triglyceride", "cholesterol",
        "phospholipid", "ceramide", "palmitate", "oleate",
        "hdl", "ldl", "vldl", "lipoprotein", "chylomicron", "ffa",
        "ketone body", "beta-hydroxybutyrate", "acetoacetate",
        "butyrate", "propionate", "acetate",
        "glutamine", "glutamate", "alanine", "leucine", "isoleucine", "valine",
        "branched-chain amino acid", "arginine", "glycine", "urea",
        "atp", "adp", "amp", "nadh", "nadph",
        "bile acid", "cholic acid", "deoxycholic acid",
        "bilirubin", "creatine", "uric acid",
    },
    "hormonal": {
        "insulin", "glucagon", "somatostatin", "amylin",
        "cortisol", "corticosterone", "aldosterone", "dhea",
        "adrenaline", "epinephrine", "noradrenaline", "norepinephrine",
        "glucocorticoid",
        "t3", "triiodothyronine", "t4", "thyroxine", "tsh", "thyroid hormone",
        "growth hormone", "igf-1", "acth", "prolactin", "vasopressin", "oxytocin",
        "leptin", "adiponectin", "resistin", "chemerin", "apelin",
        "ghrelin", "glp-1", "gip", "pyy", "cck", "secretin",
        "fgf21", "fgf19", "fetuin-a",
        "irisin", "fndc5",
        "testosterone", "estrogen", "estradiol", "progesterone",
        "angiotensin", "renin", "erythropoietin", "natriuretic peptide",
        "bnp", "anp",
        "osteocalcin", "fgf23", "klotho", "calcitriol", "vitamin d",
    },
    "immune": {
        "macrophage", "monocyte", "neutrophil", "lymphocyte",
        "t cell", "b cell", "nk cell", "dendritic cell", "mast cell",
        "eosinophil", "basophil", "treg", "th1", "th17",
        "kupffer cell", "microglia", "adipose tissue macrophage",
        "tnf-alpha", "tnf", "tumor necrosis factor",
        "il-1beta", "il-1", "interleukin-1",
        "il-6", "interleukin-6", "il-17", "interleukin-17",
        "ifn-gamma", "interferon",
        "il-10", "interleukin-10", "il-4", "il-13", "tgf-beta",
        "chemokine", "cxcl", "ccl", "mcp-1", "rantes",
        "nf-kb", "nlrp3", "inflammasome", "toll-like receptor", "tlr4",
        "complement", "reactive oxygen species", "ros", "oxidative stress",
        "lps", "lipopolysaccharide", "endotoxin",
    },
    "neural": {
        "acetylcholine", "norepinephrine", "dopamine", "serotonin",
        "gaba", "glutamate",
        "neuropeptide y", "npy", "vip", "substance p", "cgrp",
        "galanin", "neurotensin", "enkephalin",
        "sympathetic", "parasympathetic", "vagus nerve", "vagal",
        "autonomic nervous system",
        "adrenergic receptor", "muscarinic receptor", "nicotinic receptor",
        "catecholamine",
        "bdnf", "ngf", "neurotrophin", "glial cell",
    },
    "mechanical": {
        "shear stress", "blood flow", "hemodynamic", "wall stress",
        "pressure", "stretch", "strain", "compression", "tension",
        "baroreceptor", "mechanosensing", "mechanotransduction",
        "piezo1", "piezo2", "integrin",
        "vascular tone", "vasoconstriction", "vasodilation",
        "perfusion", "microcirculation", "endothelium",
        "nitric oxide", "no", "endothelin",
        "cardiac output", "preload", "afterload", "contractility",
    },
    "undefined": set(),  # uses combined vocab at extraction time
}


# ── Co-occurrence detection ────────────────────────────────────────────────────

def _organ_search_aliases(organ: str) -> list[str]:
    """Return short lowercase aliases suitable for sentence/text searching."""
    from Literature_Search.pubmed_search import ORGAN_ALIASES
    raw = ORGAN_ALIASES.get(organ, [organ])
    # Keep only single-word or short phrases (≤3 tokens) to avoid false positives
    return [a.lower() for a in raw if len(a.split()) <= 3]


def _paper_cooccurs(paper: dict, org1: str, org2: str) -> bool:
    """
    Return True if the two organs co-occur in this paper in any of:
      (a) a hyphenated compound term: "liver-gut", "hepato-renal"
      (b) a merged compound word: "hepatorenal", "cardiorenal"
      (c) the same sentence of the title or abstract
    """
    text = (paper.get("title", "") + ". " + paper.get("abstract", "")).lower()
    aliases1 = _organ_search_aliases(org1)
    aliases2 = _organ_search_aliases(org2)
    stems1   = [s.lower() for s in ORGAN_COMPOUND_STEMS.get(org1, [])]
    stems2   = [s.lower() for s in ORGAN_COMPOUND_STEMS.get(org2, [])]

    all1 = aliases1 + stems1
    all2 = aliases2 + stems2

    # (a) Hyphenated: "liver-gut", "gut-liver", "hepato-renal"
    for a in all1:
        for b in all2:
            if f"{a}-{b}" in text or f"{b}-{a}" in text:
                return True

    # (b) Merged compound word containing stems from both organs
    if stems1 and stems2:
        for word in re.findall(r'\b\w{6,}\b', text):
            has1 = any(s in word for s in stems1)
            has2 = any(s in word for s in stems2)
            if has1 and has2:
                return True

    # (c) Same sentence (split on ., !, ?, ;)
    sentences = re.split(r'(?<=[.!?;])\s+', text)
    for sent in sentences:
        has1 = any(a in sent for a in aliases1)
        has2 = any(a in sent for a in aliases2)
        if has1 and has2:
            return True

    return False


# ── Layer classification ───────────────────────────────────────────────────────

def _classify_paper(paper: dict) -> list[str]:
    """
    Classify a paper into one or more layers based on LAYER_REQUIRED_PHRASES.
    Returns ["undefined"] if no layer phrase is found.
    """
    text = (paper.get("title", "") + " " + paper.get("abstract", "")).lower()
    layers = [
        lyr for lyr, phrases in LAYER_REQUIRED_PHRASES.items()
        if any(ph in text for ph in phrases)
    ]
    return layers if layers else ["undefined"]


# ── Key player extraction ──────────────────────────────────────────────────────

def extract_layer_key_players(
    papers: list[dict], layer: str, top_n: int = 20
) -> tuple[list[str], dict[str, int]]:
    """
    Scan abstracts for vocabulary terms specific to a communication layer.
    For "undefined", uses the union of all layer vocabularies.
    Returns (terms_ranked_by_count, counts_dict), limited to top_n terms.
    """
    if layer == "undefined":
        vocab: set = set()
        for v in LAYER_VOCAB.values():
            vocab |= v
    else:
        vocab = LAYER_VOCAB.get(layer, set())

    counts: dict[str, int] = {}
    for p in papers:
        txt = (p.get("title", "") + " " + p.get("abstract", "")).lower()
        for term in vocab:
            if term in txt:
                counts[term] = counts.get(term, 0) + 1

    ranked = sorted(counts, key=lambda t: -counts[t])[:top_n]
    return ranked, {t: counts[t] for t in ranked}


# ── Filter cached papers by mechanism phrase ───────────────────────────────────

def filter_papers_by_layer(papers: list[dict], layer: str) -> list[dict]:
    """Keep only papers that contain at least one required phrase for the layer."""
    if layer == "undefined":
        return papers  # undefined has no required phrases
    phrases = LAYER_REQUIRED_PHRASES.get(layer, [])
    if not phrases:
        return papers
    return [
        p for p in papers
        if any(ph in (p.get("title", "") + " " + p.get("abstract", "")).lower()
               for ph in phrases)
    ]


# ── PubMed helpers ─────────────────────────────────────────────────────────────

def all_organ_pairs(organs: list[str]) -> list[tuple[str, str]]:
    return sorted(combinations(sorted(organs), 2))


def _search_pair(
    organ1: str, organ2: str,
    years_back: int, max_results: int, delay: float,
) -> list[dict]:
    """Broad PubMed search for papers mentioning both organs."""
    from Literature_Search.pubmed_search import (
        ORGAN_ALIASES, ORGAN_MESH, _ncbi_get, fetch_abstracts,
    )

    def organ_clause(organ: str) -> str:
        aliases = ORGAN_ALIASES.get(organ, [organ])
        mesh    = ORGAN_MESH.get(organ, f'"{organ}"[MeSH Terms]')
        alias_q = "(" + " OR ".join(f'"{a}"[Title/Abstract]' for a in aliases) + ")"
        return f"({mesh} OR {alias_q})"

    # Add a broad interaction filter to reduce completely unrelated papers
    interaction_filter = (
        '(axis[Title/Abstract] OR crosstalk[Title/Abstract] '
        'OR interaction[Title/Abstract] OR regulation[Title/Abstract] '
        'OR signaling[Title/Abstract] OR communication[Title/Abstract] '
        'OR coupling[Title/Abstract] OR cross-talk[Title/Abstract])'
    )

    query = (
        f"{organ_clause(organ1)} AND {organ_clause(organ2)}"
        f" AND {interaction_filter}"
    )

    min_date = (datetime.now() - timedelta(days=365 * years_back)).strftime("%Y/%m/%d")
    params = {
        "db":       "pubmed",
        "term":     query,
        "retmax":   max_results,
        "retmode":  "json",
        "mindate":  min_date,
        "datetype": "pdat",
    }
    resp = _ncbi_get("esearch.fcgi", params)
    if resp is None:
        return []
    pmids = resp.json().get("esearchresult", {}).get("idlist", [])
    if not pmids:
        return []
    time.sleep(delay)
    return fetch_abstracts(pmids, delay=delay)


# ── Main search runner ─────────────────────────────────────────────────────────

def run_lit_ref_search(
    organs: list[str],
    output_path: "str | Path" = DEFAULT_OUTPUT,
    min_papers: int  = DEFAULT_MIN_PAPERS,
    max_results: int = DEFAULT_MAX_RESULTS,
    years_back: int  = DEFAULT_YEARS_BACK,
    delay: float     = DELAY,
    resume: bool     = True,
    reset: bool      = False,
) -> dict:
    """
    Search all organ pairs with co-occurrence filtering and layer classification.
    Cache format: one entry per "organ1|organ2" key.
    """
    output_path = Path(output_path)
    results: dict = {}

    if reset and output_path.exists():
        output_path.unlink()
        print(f"[i] Cache deleted: {output_path.name}")
    elif resume and output_path.exists():
        with open(output_path, encoding="utf-8") as f:
            results = json.load(f)
        # Detect old format (keys contain two '|') and warn
        old_keys = [k for k in results if k.count("|") == 2]
        if old_keys:
            print(
                f"[!] Cache uses old per-layer format ({len(results)} entries). "
                "Run with --reset to rebuild with the new co-occurrence approach."
            )
            return results
        print(f"[i] Resuming: {len(results)} pairs already cached.")

    pairs  = all_organ_pairs(organs)
    total  = len(pairs)

    for idx, (o1, o2) in enumerate(pairs, 1):
        key = f"{o1}|{o2}"
        if key in results:
            d  = results[key]
            nc = d.get("n_papers_cooccur", 0)
            print(f"  [{idx:3d}/{total}] cached ({nc} co-occur): {o1} ↔ {o2}")
            continue

        print(f"  [{idx:3d}/{total}] {o1} ↔ {o2} … ", end="", flush=True)
        papers = _search_pair(o1, o2, years_back, max_results, delay)
        time.sleep(delay)

        # Post-filter: co-occurrence in same sentence or compound term
        cooccur = [p for p in papers if _paper_cooccurs(p, o1, o2)]

        # Classify each co-occurring paper into layers
        layer_papers: dict[str, list] = {lyr: [] for lyr in ALL_LAYERS}
        for p in cooccur:
            for lyr in _classify_paper(p):
                if lyr in layer_papers:
                    layer_papers[lyr].append(p)

        layer_counts = {lyr: len(ps) for lyr, ps in layer_papers.items()}
        active_layers = [lyr for lyr, n in layer_counts.items() if n >= 1]

        print(
            f"{len(papers)} fetched → {len(cooccur)} co-occur "
            f"| layers: {', '.join(f'{lyr}:{layer_counts[lyr]}' for lyr in active_layers) or 'none'}"
        )

        results[key] = {
            "organ1":           o1,
            "organ2":           o2,
            "n_papers_fetched": len(papers),
            "n_papers_cooccur": len(cooccur),
            "layer_papers":     {lyr: ps[:20] for lyr, ps in layer_papers.items()},
            "layer_counts":     layer_counts,
            "search_date":      datetime.now().isoformat(),
        }

        with open(output_path, "w", encoding="utf-8") as f:
            json.dump(results, f, indent=2, ensure_ascii=False)

    n_edges = sum(
        1 for v in results.values()
        if isinstance(v, dict) and v.get("n_papers_cooccur", 0) >= min_papers
    )
    print(
        f"\n[ok] Done. {n_edges}/{total} pairs have ≥{min_papers} co-occurring papers."
        f"\n     Saved to: {output_path}"
    )
    return results


# ── Standalone entry point ─────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Literature-derived organ network: co-occurrence + layer classification."
    )
    parser.add_argument("--reset",       action="store_true", help="Delete cache and restart.")
    parser.add_argument("--min-papers",  type=int, default=DEFAULT_MIN_PAPERS)
    parser.add_argument("--max-results", type=int, default=DEFAULT_MAX_RESULTS)
    parser.add_argument("--years-back",  type=int, default=DEFAULT_YEARS_BACK)
    args = parser.parse_args()

    from Data_Loader.load_data import load_node_metadata_from_csv
    organ_data    = HERE / "metabolic_data" / "organ_data.csv"
    node_metadata = load_node_metadata_from_csv(str(organ_data))
    organs        = list(node_metadata.keys())
    print(f"[i] {len(organs)} organs → {len(organs)*(len(organs)-1)//2} pairs.")

    run_lit_ref_search(
        organs      = organs,
        output_path = DEFAULT_OUTPUT,
        min_papers  = args.min_papers,
        max_results = args.max_results,
        years_back  = args.years_back,
        resume      = not args.reset,
        reset       = args.reset,
    )


if __name__ == "__main__":
    main()
