"""
run_network.py
==============
Condition-specific bootstrapped general edge network, using an
Otsuka–Ochiai (cosine) coefficient for edge strength instead of a raw
cross-mention count.

For each organ-organ pair the pipeline runs TWO queries in parallel:

  Query A (metabolic layer):
      (MeSH_A OR aliases_A) AND (MeSH_B OR aliases_B)
      AND METABOLIC_FILTER AND CONDITION_FILTER

  Query B (hormonal layer):
      (MeSH_A OR aliases_A) AND (MeSH_B OR aliases_B)
      AND HORMONAL_FILTER AND CONDITION_FILTER

  Papers from both queries are merged (union by PMID).  Each paper is tagged
  with the layer(s) it was found in — "metabolic", "hormonal", or both.

  3. Co-occurrence filter on the combined pool
  4. Bootstraps: N_BOOTSTRAP × SAMPLE_FRACTION random sample → mean cross-mention
     count per iteration. That raw mean is then normalized into an
     Otsuka–Ochiai coefficient:
         OO(A,B) = mean_cross_mentions / sqrt(n_found_A * n_found_B)
     A raw count is dominated by how much is published about an organ at
     all — this corrects for that, so edges reflect how specifically two
     organs are tied together rather than how well-studied either one is.
  5. Overview figure (heatmap + bar chart) and robust network — edges are
     kept when their mean Otsuka–Ochiai coefficient ≥ MIN_BOOTSTRAP_MEAN.

Run from the Metabolic_Reference_Network/ directory:
    uv run -m Edge_cosine_general_reference_network.run_network --condition healthy
    uv run -m Edge_cosine_general_reference_network.run_network --condition obese

    --reset        wipe search + bootstrap cache, start over
    --viz-only     skip search + bootstrap, rebuild HTML only
"""

import sys
import re
import json
import time
import random
import math
import argparse
import importlib.util
from pathlib import Path
from datetime import datetime, timedelta

import requests

HERE = Path(__file__).parent          # Edge_cosine_general_reference_network/
ROOT = HERE.parent                     # Metabolic_Reference_Network/
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from threshold_utils import Elbow, is_elbow, kneedle_elbow, resolve_min_bootstrap_mean
import run_comparison as _cmp  # sibling module — reused to embed the comparison view

NCBI_BASE  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

CONDITION_CONFIGS = {
    "healthy": HERE / "healthy" / "config_healthy.py",
    "obese":   HERE / "obese"   / "config_obese.py",
}

# Cohort connections CSV that defines which organs are in scope for each
# condition — only organs listed in its header row are searched/bootstrapped.
COHORT_CONNECTIONS_CSV = {
    "healthy": HERE / "healthy" / "healthy_cohort_connections.csv",
    "obese":   HERE / "obese"   / "obese_cohort_connections.csv",
}


# ── Config loader ─────────────────────────────────────────────────────────────

def _load_config(condition: str):
    path = CONDITION_CONFIGS[condition]
    spec = importlib.util.spec_from_file_location(f"_cfg_{condition}", path)
    mod  = importlib.util.module_from_spec(spec)
    mod.Elbow = Elbow  # lets the config file write `MIN_BOOTSTRAP_MEAN = Elbow`
    spec.loader.exec_module(mod)
    return mod


def load_cohort_organs(csv_path: Path) -> list[str]:
    """Return the organ names from the header row of a cohort connections CSV."""
    import csv as _csv
    with open(csv_path, encoding="utf-8") as f:
        header = next(_csv.reader(f))
    return [h.strip() for h in header[1:] if h.strip()]


# ── NCBI HTTP helper ──────────────────────────────────────────────────────────

def _ncbi_get(endpoint: str, params: dict, retries: int = 3, method: str = "GET"):
    params.setdefault("tool", "MetabolicReferenceNetwork")
    params.setdefault("email", "research@metabolic-network.org")
    for attempt in range(retries):
        try:
            if method == "POST":
                resp = requests.post(NCBI_BASE + endpoint, data=params, timeout=30)
            else:
                resp = requests.get(NCBI_BASE + endpoint, params=params, timeout=30)
            resp.raise_for_status()
            return resp
        except requests.RequestException as e:
            wait = 2 ** attempt
            print(f"    [!] HTTP error (attempt {attempt+1}/{retries}): {e} — retry in {wait}s")
            time.sleep(wait)
    return None


# ── Query builders ────────────────────────────────────────────────────────────

def _keyword_clause(keywords: list[str]) -> str:
    parts = []
    for kw in keywords:
        if " " in kw or "-" in kw:
            parts.append(f'"{kw}"[Title/Abstract]')
        else:
            parts.append(f'{kw}[Title/Abstract]')
    return "(" + " OR ".join(parts) + ")"


def _alias_clause(organ: str, organ_aliases: dict, organ_mesh: dict) -> str:
    mesh        = organ_mesh.get(organ, f'"{organ}"[MeSH Terms]')
    aliases     = organ_aliases.get(organ, [organ])
    alias_terms = " OR ".join(f'"{a}"[Title/Abstract]' for a in aliases)
    return f"({mesh} OR {alias_terms})"


def build_query(organ1: str, organ2: str,
                organ_aliases: dict, organ_mesh: dict,
                layer_filter: str, condition_filter: str,
                crosstalk_filter: str) -> str:
    """Build a PubMed query for one layer (metabolic or hormonal)."""
    c1 = _alias_clause(organ1, organ_aliases, organ_mesh)
    c2 = _alias_clause(organ2, organ_aliases, organ_mesh)
    return f"{c1} AND {c2} AND {layer_filter} AND {condition_filter} AND {crosstalk_filter}"


def build_per_organ_query(organ: str, organ_aliases: dict, organ_mesh: dict,
                          layer_filter: str, condition_filter: str,
                          crosstalk_filter: str) -> str:
    clause = _alias_clause(organ, organ_aliases, organ_mesh)
    return f"{clause} AND {layer_filter} AND {condition_filter} AND {crosstalk_filter}"


def _organ_aliases_set(organ: str, organ_aliases: dict) -> set[str]:
    """Lowercase alias set for cross-mention text search."""
    return {a.lower() for a in organ_aliases.get(organ, [organ])}


# ── Word-boundary / negation-safe term matching ────────────────────────────
# Plain substring checks (`"renal" in text`) false-positive on words that
# merely contain the term, e.g. "renal" inside "adrenal" — a paper about the
# adrenal gland would otherwise count as kidney evidence. \b...\b anchors the
# match to whole-word boundaries so this can't happen. A leading negative
# lookbehind additionally excludes negated mentions ("non-renal",
# "nonrenal", "non renal", "not renal") from counting as a positive mention
# of that organ/keyword — these explicitly say the opposite.

def _word_boundary_pattern(term: str):
    escaped = re.escape(term.lower())
    return re.compile(r"(?<!non-)(?<!non )(?<!non)(?<!not )\b" + escaped + r"\b")


def _compile_patterns(terms) -> list:
    return [_word_boundary_pattern(t) for t in terms]


def _any_match(patterns: list, text: str) -> bool:
    return any(p.search(text) for p in patterns)


def _filter_same_sentence_crosstalk(papers: list[dict], organ_patterns: list,
                                    crosstalk_patterns: list) -> list[dict]:
    """
    Keep only papers where this organ's own alias and at least one crosstalk
    keyword (network, axis, interplay, ...) appear together in the same
    sentence — not just somewhere in the same title/abstract. The PubMed
    query already requires the crosstalk filter to match somewhere in the
    document (cheap pre-filter, fewer papers fetched); this is the strict
    local check that actually enforces co-location.
    """
    kept = []
    for paper in papers:
        text = (paper.get("title", "") + ". " + paper.get("abstract", "")).lower()
        sentences = re.split(r"(?<=[.!?;])\s+", text)
        for sent in sentences:
            if _any_match(organ_patterns, sent) and _any_match(crosstalk_patterns, sent):
                kept.append(paper)
                break
    return kept


# ── PubMed search ─────────────────────────────────────────────────────────────

def search_pubmed(query: str, max_results: int, years_back: int) -> list[str]:
    min_date = (datetime.now() - timedelta(days=365 * years_back)).strftime("%Y/%m/%d")
    params = {
        "db": "pubmed", "term": query, "retmax": max_results,
        "retmode": "json", "mindate": min_date, "datetype": "pdat",
    }
    # POST avoids the URL-length ceiling GET requests hit once the
    # crosstalk/hormonal keyword clauses push the query past a few thousand
    # characters — NCBI returns HTTP 414 on GET in that case, which silently
    # degrades that organ's results to zero rather than erroring loudly.
    resp = _ncbi_get("esearch.fcgi", params, method="POST")
    if resp is None:
        return []
    return resp.json().get("esearchresult", {}).get("idlist", [])


def fetch_abstracts(pmids: list[str], delay: float = 0.4,
                    batch_size: int = 200) -> list[dict]:
    papers = []
    for i in range(0, len(pmids), batch_size):
        batch  = pmids[i:i + batch_size]
        params = {
            "db": "pubmed", "id": ",".join(batch),
            "rettype": "abstract", "retmode": "xml",
        }
        resp = _ncbi_get("efetch.fcgi", params)
        if resp is not None:
            papers.extend(_parse_xml_abstracts(resp.text))
        time.sleep(delay)
    return papers


def _parse_xml_abstracts(xml: str) -> list[dict]:
    import re
    articles = re.findall(r"<PubmedArticle>(.*?)</PubmedArticle>", xml, re.DOTALL)
    results  = []
    for art in articles:
        pmid      = _tag(art, "PMID") or ""
        title     = _tag(art, "ArticleTitle") or ""
        year      = _tag(art, "Year") or _tag(art, "MedlineDate") or ""
        abstracts = re.findall(r"<AbstractText[^>]*>(.*?)</AbstractText>", art, re.DOTALL)
        abstract  = " ".join(abstracts)
        title     = re.sub(r"<[^>]+>", "", title).strip()
        abstract  = re.sub(r"<[^>]+>", "", abstract).strip()
        if pmid and (title or abstract):
            results.append({"pmid": pmid, "title": title,
                            "abstract": abstract, "year": year})
    return results


def _tag(text: str, tag: str) -> str:
    import re
    m = re.search(rf"<{tag}[^>]*>(.*?)</{tag}>", text, re.DOTALL)
    return m.group(1).strip() if m else ""


# ── Co-occurrence filter ──────────────────────────────────────────────────────

def _filter_cooccur(papers: list[dict], organ1: str, organ2: str) -> list[dict]:
    from Literature_Reference_Network.lit_ref_search import _paper_cooccurs
    return [p for p in papers if _paper_cooccurs(p, organ1, organ2)]


# ── Key player vocabularies ───────────────────────────────────────────────────
# Proteins relevant to inter-organ metabolic communication.
# These are scanned in title + abstract (lowercased).

PROTEIN_VOCAB: set[str] = {
    # Glucose transporters
    "glut1", "glut2", "glut3", "glut4", "slc2a1", "slc2a2", "slc2a4",
    # Monocarboxylate transporters
    "mct1", "mct4", "slc16a1",
    # Fatty acid transport
    "cd36", "fatp1", "fatp4", "fabp",
    # Insulin / growth signalling
    "insulin receptor", "insr", "irs-1", "irs1", "irs2",
    "pi3k", "akt", "pdk1", "pkb",
    "mtor", "mtorc1", "mtorc2", "raptor", "rictor",
    "ampk", "lkb1",
    # Glucagon / cAMP
    "glucagon receptor", "gcgr", "pka", "creb",
    # Nuclear receptors / transcription factors
    "ppar-alpha", "ppara", "ppar-gamma", "pparg",
    "pgc-1alpha", "pgc-1α", "pgc1a", "pgc-1",
    "srebp", "srebp-1c", "chrebp",
    "foxo1", "foxo3",
    "lxr", "rxr", "hnf4a",
    # Lipid metabolism enzymes
    "cpt1", "cpt-1", "cpt2",
    "lpl", "lipoprotein lipase",
    "fasn", "fatty acid synthase",
    "acc", "acetyl-coa carboxylase",
    "atgl", "hsl",
    # Glucose metabolism enzymes
    "hexokinase", "hk2",
    "pyruvate kinase", "pkm2",
    "pyruvate dehydrogenase", "pdh", "pdha",
    "g6pase", "pepck", "pck1",
    "pfkfb3", "pfk",
    # Mitochondrial
    "uncoupling protein", "ucp1", "ucp2", "ucp3",
    # Fibroblast growth factors
    "fgf21", "fgf19", "fgf15",
    # Other metabolic regulators
    "sirt1", "sirt3",
    "glp-1 receptor", "glp1r",
    "adiponectin receptor", "adipor1", "adipor2",
    "leptin receptor", "lepr",
    "glucokinase", "gck",
}


def _canon_map(synonym_groups: dict[str, set[str]]) -> dict[str, str]:
    """raw term -> canonical display name, from a pubmed_search.py-style
    {canonical: {raw variants}} synonym dict. Terms not covered by any group
    are left to fall back to themselves at lookup time."""
    return {term: canon for canon, terms in synonym_groups.items() for term in terms}


def _precompute_paper_terms(papers: list[dict]) -> dict[str, dict]:
    """
    Pre-scan each co-occurring paper for hormone, metabolite, and protein
    terms. Synonym variants (e.g. "il-6"/"interleukin-6", "mtor"/"mtorc1") are
    merged to a single canonical name per paper — reusing
    Literature_Search.pubmed_search's HORMONE_SYNONYMS/METABOLITE_SYNONYMS/
    PROTEIN_SYNONYMS — so a paper using two synonyms of the same term counts
    as one key-player mention there, not two.
    Returns {pmid: {category: [matched_canonical_terms]}} for fast lookup in bootstrap.
    """
    from Literature_Reference_Network.lit_ref_search import LAYER_VOCAB
    from Literature_Search.pubmed_search import (
        HORMONE_SYNONYMS, METABOLITE_SYNONYMS, PROTEIN_SYNONYMS,
    )

    hormone_vocab    = LAYER_VOCAB["hormonal"]
    metabolite_vocab = LAYER_VOCAB["metabolic"]

    hormone_canon    = _canon_map(HORMONE_SYNONYMS)
    metabolite_canon = _canon_map(METABOLITE_SYNONYMS)
    protein_canon    = _canon_map(PROTEIN_SYNONYMS)

    def _matched_canonical(vocab, canon_map, text):
        # dict.fromkeys dedupes while preserving first-seen order — two
        # synonym variants matching the same paper collapse to one entry.
        return list(dict.fromkeys(canon_map.get(t, t) for t in vocab if t in text))

    result = {}
    for paper in papers:
        pmid = paper.get("pmid", "")
        if not pmid:
            continue
        text = (paper.get("title", "") + " " + paper.get("abstract", "")).lower()
        result[pmid] = {
            "hormones":    _matched_canonical(hormone_vocab, hormone_canon, text),
            "metabolites": _matched_canonical(metabolite_vocab, metabolite_canon, text),
            "proteins":    _matched_canonical(PROTEIN_VOCAB, protein_canon, text),
        }
    return result


def _kp_bootstrap_to_ranked(kp_accum: dict, n_iters: int, top_n: int = 15) -> dict:
    """
    Convert accumulated {cat: {term: {sum, appear_count}}} into
    {cat: [{term, mean, freq}]} sorted by frequency, top_n each.
    """
    result = {}
    for cat, terms in kp_accum.items():
        scored = []
        for term, acc in terms.items():
            mean_count = acc["sum"] / n_iters
            freq       = acc["appear_count"] / n_iters
            scored.append({"term": term, "mean": round(mean_count, 4),
                           "freq": round(freq, 4)})
        scored.sort(key=lambda x: -x["freq"])
        result[cat] = scored[:top_n]
    return result


def _format_kp_labels(kp_ranked: dict) -> dict:
    """Convert bootstrap-ranked key players to label strings for the viz."""
    return {
        cat: [f"{e['term']} ({e['freq']*100:.0f}%)" for e in entries]
        for cat, entries in kp_ranked.items()
    }


# ── Search loop ───────────────────────────────────────────────────────────────

def run_search(organs: list[str], output_path: Path,
               organ_aliases: dict, organ_mesh: dict,
               metabolic_filter: str, hormonal_filter: str,
               condition_filter: str, crosstalk_filter: str,
               crosstalk_keywords: list[str], cfg,
               resume: bool = True) -> dict:
    """
    Per-organ PubMed search: runs a metabolic query AND a hormonal query per
    organ, merges the paper pools (union by PMID), and tags each paper with its
    source layer(s).  Returns {organ: {query_m, query_h, pmids, papers, ...}}.

    The PubMed queries' crosstalk clause only requires a crosstalk keyword to
    appear somewhere in the document — a cheap pre-filter. After fetching,
    papers are further filtered down to those where the organ's own alias
    and a crosstalk keyword appear together in the SAME SENTENCE
    (word-boundary/negation-safe — see _filter_same_sentence_crosstalk).
    """
    results: dict = {}
    if resume and output_path.exists():
        with open(output_path, encoding="utf-8") as f:
            results = json.load(f)
        print(f"[i] Resuming: {len(results)}/{len(organs)} organs cached.")

    crosstalk_patterns = _compile_patterns(crosstalk_keywords)

    total = len(organs)
    for idx, organ in enumerate(organs):
        if organ in results:
            r = results[organ]
            print(f"  [{idx+1}/{total}] Cached {organ}: "
                  f"{r.get('n_found', 0)} papers "
                  f"(M={r.get('n_found_metabolic', 0)}, H={r.get('n_found_hormonal', 0)})")
            continue

        query_m = build_per_organ_query(organ, organ_aliases, organ_mesh,
                                        metabolic_filter, condition_filter,
                                        crosstalk_filter)
        query_h = build_per_organ_query(organ, organ_aliases, organ_mesh,
                                        hormonal_filter, condition_filter,
                                        crosstalk_filter)
        print(f"\n  [{idx+1}/{total}] {organ}")
        time.sleep(cfg.DELAY)

        pmids_m = search_pubmed(query_m, cfg.MAX_PAPERS, cfg.YEARS_BACK)
        time.sleep(cfg.DELAY)
        pmids_h = search_pubmed(query_h, cfg.MAX_PAPERS, cfg.YEARS_BACK)

        set_m, set_h   = set(pmids_m), set(pmids_h)
        combined_pmids = list(set_m | set_h)
        n_raw          = len(combined_pmids)
        print(f"    metabolic={len(pmids_m)}  hormonal={len(pmids_h)}  "
              f"combined={n_raw}", end="", flush=True)

        if n_raw == 0:
            print()
            papers = []
        else:
            print(" — fetching abstracts...", flush=True)
            papers_raw = fetch_abstracts(combined_pmids, delay=cfg.DELAY)
            for paper in papers_raw:
                pmid   = paper.get("pmid", "")
                layers = []
                if pmid in set_m:
                    layers.append("metabolic")
                if pmid in set_h:
                    layers.append("hormonal")
                paper["layers"] = layers
            organ_patterns = _compile_patterns(organ_aliases.get(organ, [organ]))
            papers = _filter_same_sentence_crosstalk(papers_raw, organ_patterns, crosstalk_patterns)
            print(f"    => {len(papers_raw)} abstracts retrieved, "
                  f"{len(papers)} kept (organ + crosstalk term in same sentence)")

        results[organ] = {
            "organ":               organ,
            "query_metabolic":     query_m,
            "query_hormonal":      query_h,
            "pmids":               [p.get("pmid", "") for p in papers],
            "papers":              papers,
            "n_found":             len(papers),
            "n_found_raw":         n_raw,
            "n_found_metabolic":   len(pmids_m),
            "n_found_hormonal":    len(pmids_h),
            "search_date":         datetime.now().isoformat(),
        }

        with open(output_path, "w", encoding="utf-8") as f:
            json.dump(results, f, indent=2, ensure_ascii=False)

    print(f"\n[ok] {len(results)} organs | "
          f"{sum(v.get('n_found', 0) for v in results.values())} total papers")
    return results


# ── Bootstrap ─────────────────────────────────────────────────────────────────

def run_bootstrap(search_results: dict, output_path: Path, cfg,
                  organ_aliases: dict,
                  resume: bool = True) -> dict:
    """
    Per iteration: sample SAMPLE_FRACTION of each organ's pool and count how
    often each other organ is mentioned.  Layer tags on papers drive the stacked
    bar chart in the overview.  The mean cross-mention count is then normalized
    into an Otsuka–Ochiai coefficient — mean / sqrt(n_found_o1 * n_found_o2) —
    so an organ with a much larger literature pool doesn't dominate every edge
    just by volume.  mean Otsuka–Ochiai coefficient >= MIN_BOOTSTRAP_MEAN
    → robust edge.
    """
    from itertools import combinations

    boot: dict = {}
    if resume and output_path.exists():
        with open(output_path, encoding="utf-8") as f:
            boot = json.load(f)
        print(f"[i] Resuming bootstrap: {len(boot)} pairs cached.")

    organs  = sorted(search_results.keys())
    pairs   = list(combinations(sorted(organs), 2))
    n_pairs = len(pairs)
    all_cached = boot and all(
        f"{o1}|{o2}" in boot or f"{o2}|{o1}" in boot for o1, o2 in pairs
    )
    if all_cached and resume:
        print(f"[i] All {n_pairs} pairs already cached.")
        return boot

    alias_sets = {organ: _compile_patterns(_organ_aliases_set(organ, organ_aliases))
                  for organ in organs}

    # Pre-scan each organ's pool: for each paper, which other organs are
    # mentioned in the SAME SENTENCE as this organ's own alias? A cross-mention
    # requires both organs' aliases to co-occur within one sentence of the
    # title/abstract — a whole-document match no longer counts (e.g. an organ
    # named once in the intro and an unrelated organ named once in the
    # discussion no longer creates an edge).
    print("[i] Pre-scanning organ paper pools for cross-organ mentions (same-sentence)...")
    organ_paper_mentions: dict[str, dict[str, set]] = {}
    for organ, data in search_results.items():
        own_aliases = alias_sets[organ]
        pm: dict[str, set] = {}
        for paper in data.get("papers", []):
            pmid = paper.get("pmid", "")
            if not pmid:
                continue
            # ". " between title and abstract guarantees a sentence boundary
            # there, so a title-only alias can never falsely co-locate with
            # an abstract-only alias just because they were concatenated.
            text = (paper.get("title", "") + ". " + paper.get("abstract", "")).lower()
            sentences = re.split(r"(?<=[.!?;])\s+", text)
            mentioned: set = set()
            for sent in sentences:
                if not _any_match(own_aliases, sent):
                    continue
                for other, patterns in alias_sets.items():
                    if other == organ or other in mentioned:
                        continue
                    if _any_match(patterns, sent):
                        mentioned.add(other)
            pm[pmid] = mentioned
        organ_paper_mentions[organ] = pm
        print(f"  {organ}: {len(data.get('papers', []))} papers scanned")

    # Pre-compute key-player terms for all papers
    all_paper_terms: dict[str, dict] = {}
    for data in search_results.values():
        all_paper_terms.update(_precompute_paper_terms(data.get("papers", [])))

    # pmid → paper and pmid → layer set for viz and stacked bar
    pmid_to_paper:  dict[str, dict]     = {}
    pmid_to_layers: dict[str, set[str]] = {}
    for data in search_results.values():
        for paper in data.get("papers", []):
            pmid = paper.get("pmid", "")
            if pmid:
                pmid_to_paper[pmid] = paper
                for lyr in paper.get("layers", []):
                    pmid_to_layers.setdefault(pmid, set()).add(lyr)

    # Per-pair accumulators
    pair_keys    = [f"{o1}|{o2}" for o1, o2 in pairs]
    pair_counts  = {k: [] for k in pair_keys}
    pair_kp      = {k: {cat: {} for cat in ("hormones", "metabolites", "proteins")}
                    for k in pair_keys}
    pair_cooccur = {k: set() for k in pair_keys}

    print(f"[i] Running {cfg.N_BOOTSTRAP} bootstrap iterations "
          f"({int(cfg.SAMPLE_FRACTION * 100)}% per organ)...")
    for it in range(cfg.N_BOOTSTRAP):
        organ_samples: dict[str, list] = {}
        for organ, data in search_results.items():
            pmids    = data.get("pmids", [])
            n_sample = max(1, round(len(pmids) * cfg.SAMPLE_FRACTION))
            organ_samples[organ] = random.sample(pmids, min(n_sample, len(pmids)))

        for (o1, o2), key in zip(pairs, pair_keys):
            count     = 0
            hit_pmids: set = set()

            pm1 = organ_paper_mentions.get(o1, {})
            for pmid in organ_samples.get(o1, []):
                if o2 in pm1.get(pmid, set()):
                    count += 1
                    hit_pmids.add(pmid)

            pm2 = organ_paper_mentions.get(o2, {})
            for pmid in organ_samples.get(o2, []):
                if o1 in pm2.get(pmid, set()):
                    count += 1
                    hit_pmids.add(pmid)

            pair_counts[key].append(count)
            pair_cooccur[key] |= hit_pmids

            iter_kp: dict = {cat: {} for cat in pair_kp[key]}
            for pmid in hit_pmids:
                for cat, terms in all_paper_terms.get(pmid, {}).items():
                    for term in terms:
                        iter_kp[cat][term] = iter_kp[cat].get(term, 0) + 1
            for cat, tc in iter_kp.items():
                for term, cnt in tc.items():
                    if term not in pair_kp[key][cat]:
                        pair_kp[key][cat][term] = {"sum": 0, "appear_count": 0}
                    pair_kp[key][cat][term]["sum"]          += cnt
                    pair_kp[key][cat][term]["appear_count"] += 1

        if (it + 1) % 10 == 0:
            print(f"  [{it+1}/{cfg.N_BOOTSTRAP}] iterations done")

    print(f"\n[i] Finalising {n_pairs} pairs...")
    for (o1, o2), key in zip(pairs, pair_keys):
        if key in boot and resume:
            continue
        counts   = pair_counts[key]
        raw_mean = sum(counts) / len(counts)
        raw_std  = math.sqrt(sum((c - raw_mean) ** 2 for c in counts) / len(counts))
        cooccur  = pair_cooccur[key]
        papers   = [pmid_to_paper[p] for p in cooccur if p in pmid_to_paper]
        kp       = _kp_bootstrap_to_ranked(pair_kp[key], cfg.N_BOOTSTRAP)

        # Layer breakdown for stacked bar chart
        n_both = sum(1 for p in cooccur
                     if "metabolic" in pmid_to_layers.get(p, set())
                     and "hormonal"  in pmid_to_layers.get(p, set()))
        n_m    = sum(1 for p in cooccur if "metabolic" in pmid_to_layers.get(p, set()))
        n_h    = sum(1 for p in cooccur if "hormonal"  in pmid_to_layers.get(p, set()))

        n1 = search_results.get(o1, {}).get("n_found", 0)
        n2 = search_results.get(o2, {}).get("n_found", 0)

        # Otsuka–Ochiai coefficient: normalizes the raw cross-mention count by
        # the geometric mean of each organ's full literature pool size, so an
        # organ that's simply heavily studied overall doesn't dominate edges
        # by volume alone.  OO = mean_cross_mentions / sqrt(n_found_o1 * n_found_o2).
        denom   = math.sqrt(n1 * n2)
        oo_mean = (raw_mean / denom) if denom > 0 else 0.0
        oo_std  = (raw_std  / denom) if denom > 0 else 0.0
        cv      = (oo_std / oo_mean) if oo_mean > 0 else 0.0

        boot[key] = {
            "organ1":                o1,
            "organ2":                o2,
            "n_found":               n1 + n2,
            "n_found_o1":            n1,
            "n_found_o2":            n2,
            "n_cooccur_total":       len(cooccur),
            "n_cooccur_metabolic":   n_m,
            "n_cooccur_hormonal":    n_h,
            "n_cooccur_both":        n_both,
            "sample_size":           (max(1, round(n1 * cfg.SAMPLE_FRACTION)) +
                                      max(1, round(n2 * cfg.SAMPLE_FRACTION))),
            "raw_counts":            counts,
            "raw_mean":              round(raw_mean, 3),
            "raw_std":               round(raw_std, 3),
            "mean":                  round(oo_mean, 6),
            "std":                   round(oo_std, 6),
            "cv":                    round(cv, 4),
            "papers":                papers,
            "key_players_bootstrap": kp,
            "bootstrap_date":        datetime.now().isoformat(),
        }
        flag  = "✓" if (not is_elbow(cfg.MIN_BOOTSTRAP_MEAN) and oo_mean >= cfg.MIN_BOOTSTRAP_MEAN) else "·"
        n_kph = len(kp.get("hormones", []))
        n_kpm = len(kp.get("metabolites", []))
        n_kpp = len(kp.get("proteins", []))
        print(f"  {flag} OO={oo_mean:.5f} ± {oo_std:.5f}  cv={cv:.3f}  "
              f"(raw mean={raw_mean:.1f})  kp: {n_kph}H/{n_kpm}M/{n_kpp}P  ({o1} <-> {o2})")

    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(boot, f, indent=2, ensure_ascii=False)

    robust = sum(1 for v in boot.values() if v["mean"] >= cfg.MIN_BOOTSTRAP_MEAN)
    print(f"\n[ok] Bootstrap done. {robust} robust pairs "
          f"(mean Otsuka–Ochiai coefficient ≥ {cfg.MIN_BOOTSTRAP_MEAN}).")
    return boot


# ── Overview HTML (heatmap + bar chart) ──────────────────────────────────────
# Elbow detection (kneedle_elbow) lives in threshold_utils.py, shared with the
# MIN_BOOTSTRAP_MEAN = Elbow config resolution.

def _build_overview_html_string(boot: dict, cfg) -> str:
    """Build the standalone bootstrap-overview HTML (heatmap + bar chart +
    elbow plot) as a string, embedded into robust_network_{condition}.html
    behind the "Bootstrap Overview" button — no separate file is written."""
    organs = sorted({v["organ1"] for v in boot.values()} |
                    {v["organ2"] for v in boot.values()})
    n = len(organs)

    lookup: dict = {}
    for v in boot.values():
        key = (min(v["organ1"], v["organ2"]), max(v["organ1"], v["organ2"]))
        lookup[key] = v

    matrix = []
    for oa in organs:
        row = []
        for ob in organs:
            if oa == ob:
                row.append({"mean": -1, "std": 0, "n_found": 0, "n_cooccur": 0})
                continue
            key = (min(oa, ob), max(oa, ob))
            v   = lookup.get(key, {})
            row.append({
                "mean":     v.get("mean", 0),
                "std":      v.get("std", 0),
                "n_found":  v.get("n_found", 0),
                "n_cooccur": v.get("n_cooccur_total", 0),
            })
        matrix.append(row)

    bar_pairs = sorted(
        [v for v in boot.values() if v["mean"] > 0],
        key=lambda x: -x["mean"]
    )

    sorted_means = [p["mean"] for p in bar_pairs]
    elbow_mean   = round(kneedle_elbow(sorted_means), 6) if sorted_means else 0.0
    elbow_idx    = next((i for i, m in enumerate(sorted_means) if m <= elbow_mean), len(sorted_means) - 1)

    data_js = {
        "organs": organs, "matrix": matrix,
        "barPairs": [
            {"label": f"{v['organ1']} — {v['organ2']}",
             "mean": v["mean"], "std": v["std"],
             "n_found": v["n_found"], "n_cooccur": v["n_cooccur_total"],
             "n_cooccur_metabolic": v.get("n_cooccur_metabolic", 0),
             "n_cooccur_hormonal":  v.get("n_cooccur_hormonal", 0),
             "n_cooccur_both":      v.get("n_cooccur_both", 0)}
            for v in bar_pairs
        ],
        "sortedMeans": sorted_means,
        "elbowMean":   elbow_mean,
        "elbowIdx":    elbow_idx,
        "threshold": cfg.MIN_BOOTSTRAP_MEAN,
        "n_bootstrap": cfg.N_BOOTSTRAP,
        "sample_fraction": cfg.SAMPLE_FRACTION,
        "title": cfg.VIZ_TITLE_OVERVIEW,
        "condition": cfg.VIZ_LABEL,
    }

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>{cfg.VIZ_TITLE_OVERVIEW}</title>
<style>
  *, *::before, *::after {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
    background: #0f172a; color: #e2e8f0; min-height: 100vh; padding: 32px 24px;
  }}
  h1 {{ font-size: 1.5rem; font-weight: 700; color: #f1f5f9; margin-bottom: 6px; }}
  .subtitle {{ color: #94a3b8; font-size: 0.88rem; margin-bottom: 28px; }}
  .card {{
    background: #1e293b; border: 1px solid #334155; border-radius: 12px;
    padding: 24px; margin-bottom: 28px;
  }}
  .card-title {{ font-size: 1rem; font-weight: 600; color: #f1f5f9; margin-bottom: 12px; }}
  .card-hint  {{ font-size: 0.8rem; color: #64748b; margin-bottom: 14px; }}
  .params {{ display: flex; gap: 20px; flex-wrap: wrap; margin-bottom: 28px; }}
  .param  {{ background: #0f172a; border-radius: 8px; padding: 12px 18px; text-align: center; }}
  .param-val {{ font-size: 1.35rem; font-weight: 700; color: #38bdf8; }}
  .param-lbl {{ font-size: 0.73rem; color: #94a3b8; margin-top: 2px; }}
  canvas {{ display: block; }}
  .legend {{ display:flex; align-items:center; gap:8px; margin-top:10px;
    font-size:0.76rem; color:#94a3b8; }}
  .legend-grad {{ width:140px; height:10px; border-radius:4px;
    background: linear-gradient(to right,#1e293b,#0ea5e9,#0369a1); }}
  .threshold-note {{ font-size:0.78rem; color:#fbbf24; margin-top:8px; }}
  .elbow-note {{ font-size:0.78rem; color:#a78bfa; margin-top:4px; }}
  .legend-line {{ display:inline-block; width:28px; height:0; border-top:2px dashed currentColor; vertical-align:middle; margin-right:6px; }}
</style>
</head>
<body>
<h1 id="pg-title"></h1>
<p class="subtitle" id="pg-sub"></p>
<div class="params" id="params-row"></div>

<div class="card">
  <div class="card-title">Connection Strength Heatmap</div>
  <div class="card-hint">Mean Otsuka–Ochiai coefficient per 10 % bootstrap sample
    (cross-mention count normalized by each organ's literature pool size).
    Hover for details. Gold outline = robust (mean ≥ threshold).</div>
  <div style="overflow-x:auto"><canvas id="heatmap"></canvas></div>
  <div class="legend">
    <span>0</span><div class="legend-grad"></div>
    <span id="leg-max"></span>
    <span style="margin-left:14px">mean Otsuka–Ochiai coefficient</span>
  </div>
  <div class="threshold-note" id="thresh-note"></div>
</div>

<div class="card">
  <div class="card-title">Connection Strength Ranking</div>
  <div class="card-hint">Sorted by mean Otsuka–Ochiai coefficient. Error bars = ± 1 SD. Bars show metabolic-only (teal), both layers (violet), and hormonal-only (pink) paper contributions per edge.</div>
  <div style="display:flex;gap:24px;flex-wrap:wrap;margin-bottom:10px;font-size:0.78rem;color:#94a3b8">
    <span><span style="display:inline-block;width:14px;height:10px;background:#0ea5e9;border-radius:2px;vertical-align:middle;margin-right:4px"></span>Metabolic only</span>
    <span><span style="display:inline-block;width:14px;height:10px;background:#a78bfa;border-radius:2px;vertical-align:middle;margin-right:4px"></span>Both layers</span>
    <span><span style="display:inline-block;width:14px;height:10px;background:#f472b6;border-radius:2px;vertical-align:middle;margin-right:4px"></span>Hormonal only</span>
    <span><span class="legend-line" style="color:#fbbf24"></span>Current threshold (<span id="lbl-threshold"></span>)</span>
    <span><span class="legend-line" style="color:#a78bfa"></span>Elbow suggestion (<span id="lbl-elbow"></span>)</span>
  </div>
  <div style="overflow-x:auto"><canvas id="barchart"></canvas></div>
</div>

<div class="card">
  <div class="card-title">Elbow Plot — Suggested Threshold</div>
  <div class="card-hint" id="elbow-hint"></div>
  <div style="overflow-x:auto"><canvas id="elbowchart"></canvas></div>
  <div class="elbow-note" id="elbow-note"></div>
</div>

<script>
const D = {json.dumps(data_js, ensure_ascii=False)};

document.getElementById('pg-title').textContent = D.title;
document.getElementById('pg-sub').textContent =
  D.condition + ' condition | ' + D.n_bootstrap + ' bootstrap iterations × ' +
  (D.sample_fraction*100).toFixed(0) + '% random sample per pair';
document.getElementById('thresh-note').textContent =
  'Gold outline = robust (mean ≥ ' + D.threshold + ')';

const robustN = D.barPairs.filter(p => p.mean >= D.threshold).length;
document.getElementById('lbl-threshold').textContent = D.threshold;
document.getElementById('lbl-elbow').textContent     = D.elbowMean;
[
  [D.n_bootstrap, 'Bootstrap iterations'],
  [(D.sample_fraction*100).toFixed(0)+'%', 'Sample fraction'],
  [D.threshold, 'Current threshold'],
  [D.elbowMean, 'Elbow suggestion'],
  [robustN, 'Robust pairs'],
].forEach(([v,l]) => {{
  const d = document.createElement('div');
  d.className='param';
  const isElbow = l === 'Elbow suggestion';
  d.innerHTML=`<div class="param-val" style="${{isElbow?'color:#a78bfa':''}}">${{v}}</div><div class="param-lbl">${{l}}</div>`;
  document.getElementById('params-row').appendChild(d);
}});

// ── Heatmap ──────────────────────────────────────────────────────────────────
const ORG=D.organs, N=ORG.length, CELL=52, LW=130, LH=110;
let maxM=0;
D.matrix.forEach(r=>r.forEach(c=>{{if(c.mean>maxM)maxM=c.mean;}}));
document.getElementById('leg-max').textContent=maxM.toFixed(4);

const hm=document.getElementById('heatmap');
hm.width=LW+N*CELL+8; hm.height=LH+N*CELL+8;
const hx=hm.getContext('2d');
hx.fillStyle='#0f172a'; hx.fillRect(0,0,hm.width,hm.height);

function cellColor(mean,max){{
  if(mean<0) return '#0f172a';
  if(mean===0) return '#1e293b';
  const t=Math.min(mean/max,1);
  return `rgb(${{Math.round(30+t*(3-30))}},${{Math.round(41+t*(105-41))}},${{Math.round(59+t*(161-59))}})`;
}}

for(let i=0;i<N;i++) for(let j=0;j<N;j++){{
  const c=D.matrix[i][j], x=LW+j*CELL, y=LH+i*CELL;
  hx.fillStyle=cellColor(c.mean,maxM);
  hx.fillRect(x+1,y+1,CELL-2,CELL-2);
  if(c.mean>=D.threshold){{
    hx.strokeStyle='#fbbf24'; hx.lineWidth=2;
    hx.strokeRect(x+2,y+2,CELL-4,CELL-4);
  }}
  if(c.mean>0){{
    hx.fillStyle=c.mean>=D.threshold?'#fef9c3':'#e2e8f0';
    hx.font='bold 9px system-ui'; hx.textAlign='center'; hx.textBaseline='middle';
    hx.fillText(c.mean.toFixed(3),x+CELL/2,y+CELL/2);
  }}
}}
hx.fillStyle='#94a3b8'; hx.font='11px system-ui'; hx.textAlign='right'; hx.textBaseline='middle';
for(let i=0;i<N;i++) hx.fillText(ORG[i],LW-6,LH+i*CELL+CELL/2);
hx.textAlign='left';
for(let j=0;j<N;j++){{
  const x=LW+j*CELL+CELL/2, y=LH-6;
  hx.save(); hx.translate(x,y); hx.rotate(-Math.PI/4);
  hx.fillText(ORG[j],0,0); hx.restore();
}}

const tip=document.createElement('div');
tip.style.cssText='position:fixed;pointer-events:none;display:none;background:#1e293b;'+
  'border:1px solid #475569;border-radius:6px;padding:8px 12px;font-size:0.78rem;'+
  'color:#e2e8f0;z-index:999;line-height:1.6';
document.body.appendChild(tip);

hm.addEventListener('mousemove',e=>{{
  const r=hm.getBoundingClientRect();
  const j=Math.floor((e.clientX-r.left-LW)/CELL);
  const i=Math.floor((e.clientY-r.top-LH)/CELL);
  if(i>=0&&i<N&&j>=0&&j<N&&i!==j){{
    const c=D.matrix[i][j];
    tip.innerHTML=`<strong>${{ORG[i]}} ↔ ${{ORG[j]}}</strong><br>`+
      `Otsuka–Ochiai coefficient: <b>${{c.mean.toFixed(5)}}</b> ± ${{c.std.toFixed(5)}}<br>`+
      `Total found: ${{c.n_found}}<br>Co-occurring (total): ${{c.n_cooccur}}`;
    tip.style.display='block';
    tip.style.left=(e.clientX+14)+'px'; tip.style.top=(e.clientY-20)+'px';
  }} else tip.style.display='none';
}});
hm.addEventListener('mouseleave',()=>{{tip.style.display='none';}});

// ── Bar chart ─────────────────────────────────────────────────────────────────
const BP=D.barPairs, BH=24, BG=5, LFT=220, TP=40, RPD=60;
const bc=document.getElementById('barchart');
bc.width=LFT+560+RPD;
bc.height=TP+Math.max(BP.length,1)*(BH+BG)+30;
const bx=bc.getContext('2d');
bx.fillStyle='#0f172a'; bx.fillRect(0,0,bc.width,bc.height);

if(BP.length===0){{
  bx.fillStyle='#64748b'; bx.font='13px system-ui'; bx.textAlign='center';
  bx.fillText('No pairs with mean > 0',bc.width/2,bc.height/2);
}} else {{
  const maxB=BP[0].mean+BP[0].std, BW=bc.width-LFT-RPD;
  const sc=v=>v/(maxB*1.05)*BW;
  bx.strokeStyle='#334155'; bx.lineWidth=1;
  bx.beginPath(); bx.moveTo(LFT,TP-10); bx.lineTo(LFT,TP+BP.length*(BH+BG)); bx.stroke();

  // Elbow line (violet, drawn behind threshold line)
  const ex2=LFT+sc(D.elbowMean);
  bx.strokeStyle='#a78bfa'; bx.lineWidth=1.5; bx.setLineDash([3,3]);
  bx.beginPath(); bx.moveTo(ex2,TP-10); bx.lineTo(ex2,TP+BP.length*(BH+BG)); bx.stroke();
  bx.setLineDash([]); bx.fillStyle='#a78bfa'; bx.font='10px system-ui';
  bx.textAlign='center'; bx.fillText('elbow',ex2,TP-14);

  // Threshold line (gold)
  const tx=LFT+sc(D.threshold);
  bx.strokeStyle='#fbbf24'; bx.lineWidth=1.5; bx.setLineDash([5,4]);
  bx.beginPath(); bx.moveTo(tx,TP-10); bx.lineTo(tx,TP+BP.length*(BH+BG)); bx.stroke();
  bx.setLineDash([]); bx.fillStyle='#fbbf24'; bx.font='10px system-ui';
  bx.textAlign='center'; bx.fillText('threshold',tx,TP-4);

  BP.forEach((p,i)=>{{
    const y=TP+i*(BH+BG), robust=p.mean>=D.threshold;
    bx.fillStyle=robust?'#fef9c3':'#94a3b8';
    bx.font=(robust?'bold ':'')+'11px system-ui';
    bx.textAlign='right'; bx.textBaseline='middle';
    bx.fillText(p.label,LFT-8,y+BH/2);

    // Stacked bar: metabolic-only (teal) | both (violet) | hormonal-only (pink)
    const total=p.n_cooccur||1;
    const nb=p.n_cooccur_both||0;
    const nMOnly=Math.max(0,(p.n_cooccur_metabolic||0)-nb);
    const nHOnly=Math.max(0,(p.n_cooccur_hormonal||0)-nb);
    const barW=sc(p.mean);
    const mW=barW*(nMOnly/total), bW=barW*(nb/total), hW=barW*(nHOnly/total);
    const restW=Math.max(0,barW-mW-bW-hW);

    // metabolic-only (teal)
    bx.fillStyle=robust?'rgba(14,165,233,0.9)':'rgba(14,165,233,0.42)';
    bx.fillRect(LFT,y,mW,BH);
    // both layers (violet)
    bx.fillStyle=robust?'rgba(167,139,250,0.9)':'rgba(167,139,250,0.42)';
    bx.fillRect(LFT+mW,y,bW,BH);
    // hormonal-only (pink)
    bx.fillStyle=robust?'rgba(244,114,182,0.9)':'rgba(244,114,182,0.42)';
    bx.fillRect(LFT+mW+bW,y,hW,BH);
    // unclassified remainder (dim)
    if(restW>0){{
      bx.fillStyle=robust?'rgba(148,163,184,0.25)':'rgba(148,163,184,0.12)';
      bx.fillRect(LFT+mW+bW+hW,y,restW,BH);
    }}

    if(p.std>0){{
      const cx=LFT+sc(p.mean), ex=sc(p.std);
      bx.strokeStyle=robust?'#fbbf24':'#64748b'; bx.lineWidth=1.5;
      bx.beginPath();
      bx.moveTo(cx-ex,y+BH/2); bx.lineTo(cx+ex,y+BH/2);
      bx.moveTo(cx-ex,y+5); bx.lineTo(cx-ex,y+BH-5);
      bx.moveTo(cx+ex,y+5); bx.lineTo(cx+ex,y+BH-5);
      bx.stroke();
    }}
    bx.fillStyle='#f1f5f9'; bx.font='10px system-ui'; bx.textAlign='left';
    bx.fillText(p.mean.toFixed(4),LFT+sc(p.mean)+sc(p.std)+4,y+BH/2);
  }});
}}

// ── Elbow plot ────────────────────────────────────────────────────────────────
const SM=D.sortedMeans, EI=D.elbowIdx, EM=D.elbowMean;
const EP_PAD={{l:60,r:40,t:30,b:50}};
const ec=document.getElementById('elbowchart');
ec.width=700; ec.height=260;
const ex=ec.getContext('2d');
ex.fillStyle='#0f172a'; ex.fillRect(0,0,ec.width,ec.height);

document.getElementById('elbow-hint').textContent =
  'Each point = one organ pair, ranked by bootstrap mean (descending). ' +
  'The violet marker shows the suggested elbow — detected on a smoothed, ' +
  'log-scaled curve so the threshold captures a meaningful cluster of edges ' +
  'rather than only the very strongest pair.';
document.getElementById('elbow-note').textContent =
  'Elbow suggestion: ' + EM + '  |  Current threshold: ' + D.threshold +
  (EM !== D.threshold ? '  ← consider updating MIN_BOOTSTRAP_MEAN in config' : '  ✓ threshold matches elbow');

if(SM.length > 1){{
  const PW=ec.width-EP_PAD.l-EP_PAD.r, PH=ec.height-EP_PAD.t-EP_PAD.b;
  const maxSM=SM[0], minSM=SM[SM.length-1];
  const px=i=>(EP_PAD.l + i/(SM.length-1)*PW);
  const py=v=>(EP_PAD.t + PH - (v-minSM)/(maxSM-minSM||1)*PH);

  // Grid lines
  ex.strokeStyle='#1e293b'; ex.lineWidth=1;
  [0,0.25,0.5,0.75,1].forEach(t=>{{
    const yy=EP_PAD.t+PH*(1-t);
    ex.beginPath(); ex.moveTo(EP_PAD.l,yy); ex.lineTo(EP_PAD.l+PW,yy); ex.stroke();
    ex.fillStyle='#475569'; ex.font='10px system-ui'; ex.textAlign='right';
    ex.fillText((minSM+(maxSM-minSM)*t).toFixed(4),EP_PAD.l-6,yy+4);
  }});

  // Diagonal reference line (normalised space projected back)
  ex.strokeStyle='#334155'; ex.lineWidth=1; ex.setLineDash([4,4]);
  ex.beginPath(); ex.moveTo(px(0),py(maxSM)); ex.lineTo(px(SM.length-1),py(minSM)); ex.stroke();
  ex.setLineDash([]);

  // Elbow vertical
  ex.strokeStyle='#a78bfa'; ex.lineWidth=1.5; ex.setLineDash([3,3]);
  const epx=px(EI);
  ex.beginPath(); ex.moveTo(epx,EP_PAD.t); ex.lineTo(epx,EP_PAD.t+PH); ex.stroke();
  ex.setLineDash([]);

  // Threshold vertical
  const threshIdx=SM.findIndex(m=>m<=D.threshold);
  if(threshIdx>=0){{
    const tpx=px(threshIdx);
    ex.strokeStyle='#fbbf24'; ex.lineWidth=1.5; ex.setLineDash([5,4]);
    ex.beginPath(); ex.moveTo(tpx,EP_PAD.t); ex.lineTo(tpx,EP_PAD.t+PH); ex.stroke();
    ex.setLineDash([]);
    ex.fillStyle='#fbbf24'; ex.font='10px system-ui'; ex.textAlign='center';
    ex.fillText('threshold',tpx,EP_PAD.t-6);
  }}

  // Curve
  ex.strokeStyle='#0ea5e9'; ex.lineWidth=2;
  ex.beginPath();
  SM.forEach((m,i)=>{{ i===0 ? ex.moveTo(px(i),py(m)) : ex.lineTo(px(i),py(m)); }});
  ex.stroke();

  // Dots
  SM.forEach((m,i)=>{{
    const isElbow=(i===EI);
    ex.fillStyle=isElbow?'#a78bfa':'rgba(14,165,233,0.5)';
    ex.beginPath(); ex.arc(px(i),py(m),isElbow?5:3,0,Math.PI*2); ex.fill();
  }});

  // Elbow label
  ex.fillStyle='#a78bfa'; ex.font='bold 11px system-ui'; ex.textAlign='center';
  ex.fillText('elbow: '+EM, epx, EP_PAD.t-6);

  // Axes labels
  ex.fillStyle='#64748b'; ex.font='11px system-ui'; ex.textAlign='center';
  ex.fillText('Organ pairs (rank)', EP_PAD.l+PW/2, ec.height-8);
  ex.save(); ex.translate(14, EP_PAD.t+PH/2); ex.rotate(-Math.PI/2);
  ex.fillText('Otsuka–Ochiai coefficient', 0, 0); ex.restore();
}}

bc.addEventListener('mousemove',e=>{{
  if(!BP.length) return;
  const r=bc.getBoundingClientRect();
  const i=Math.floor((e.clientY-r.top-TP)/(BH+BG));
  if(i>=0&&i<BP.length){{
    const p=BP[i];
    const nb2=p.n_cooccur_both||0;
    const nMOnly2=Math.max(0,(p.n_cooccur_metabolic||0)-nb2);
    const nHOnly2=Math.max(0,(p.n_cooccur_hormonal||0)-nb2);
    tip.innerHTML=`<strong>${{p.label}}</strong><br>`+
      `Otsuka–Ochiai coefficient: <b>${{p.mean.toFixed(5)}}</b> ± ${{p.std.toFixed(5)}}<br>`+
      `Total found: ${{p.n_found}}<br>Co-occurring (total): ${{p.n_cooccur}}<br>`+
      `<span style="color:#38bdf8">Metabolic only: ${{nMOnly2}}</span><br>`+
      `<span style="color:#c4b5fd">Both layers: ${{nb2}}</span><br>`+
      `<span style="color:#f472b6">Hormonal only: ${{nHOnly2}}</span><br>`+
      (p.mean>=D.threshold
        ?'<span style="color:#fbbf24">✓ Robust</span>'
        :'<span style="color:#94a3b8">Below threshold</span>');
    tip.style.display='block';
    tip.style.left=(e.clientX+14)+'px'; tip.style.top=(e.clientY-20)+'px';
  }} else tip.style.display='none';
}});
bc.addEventListener('mouseleave',()=>{{tip.style.display='none';}});
</script>
</body>
</html>
"""
    return html


# ── Robust network visualization ──────────────────────────────────────────────
# The bootstrap overview is embedded directly into robust_network_{condition}.html
# (inline content behind its top-bar button, see
# export_network_to_cytoscape_dashboard's `html`-mode extra_topbar_buttons) — no
# separate bootstrap_overview_*.html file is written. The comparison-vs-reference
# view is no longer a separate page/overlay either: comparison_data (from
# run_comparison.compute_comparison_data) is used to recolor this same graph's
# edges and swap their descriptions in place via the comparison toggle button
# (see comparison_toggle_label in export_network_to_cytoscape_dashboard).


def build_robust_network(boot: dict, search_results: dict, output_path: Path, cfg,
                         organs: list[str], overview_html_str: str,
                         comparison_data: dict | None):
    try:
        import networkx as nx
        from Visualisation.networkBuilderUtils import export_network_to_cytoscape_dashboard
        from Literature_Search.llm_connection_type import (
            load_connection_type_classifications, get_connection_types,
        )

        # Optional LLM-generated summaries (run_llm_descriptions.py) — no-op
        # if that script hasn't been run for this condition yet. Organ-level
        # text is grounded in search_results (per-organ query, same as edge
        # definition); pair-level text and connection type are grounded in
        # bootstrap "papers" (already same-sentence cross-mention evidence —
        # see run_bootstrap()).
        out_dir = output_path.parent
        organ_llm: dict = {}
        pair_llm: dict = {}
        organ_llm_path = out_dir / f"llm_organ_descriptions_{cfg.CONDITION_NAME}.json"
        pair_llm_path  = out_dir / f"llm_pair_descriptions_{cfg.CONDITION_NAME}.json"
        type_path      = out_dir / f"connection_types_{cfg.CONDITION_NAME}.json"
        if organ_llm_path.exists():
            with open(organ_llm_path, encoding="utf-8") as f:
                organ_llm = json.load(f)
        if pair_llm_path.exists():
            with open(pair_llm_path, encoding="utf-8") as f:
                pair_llm = json.load(f)
        llm_types = load_connection_type_classifications(type_path)

        G = nx.Graph()
        for organ in organs:
            G.add_node(organ, llm_description=organ_llm.get(organ, {}).get("description", ""))

        categories = comparison_data["categories"] if comparison_data else {}
        ref_edges  = comparison_data["ref_edges"]   if comparison_data else {}
        lit_label  = f"{cfg.VIZ_LABEL} Literature-Based Network"

        # Every organ pair with bootstrap data becomes a graph edge — not
        # just ones clearing MIN_BOOTSTRAP_MEAN at build time — so the
        # dashboard's threshold slider can reveal/hide edges live instead of
        # a single cutoff being baked in permanently.
        n_robust = 0
        for edge_key, b in boot.items():
            o1, o2  = b["organ1"], b["organ2"]
            papers  = b.get("papers", [])
            if b["mean"] >= cfg.MIN_BOOTSTRAP_MEAN:
                n_robust += 1

            kp_ranked = b.get("key_players_bootstrap", {})

            def _kp_lists_and_counts(cat: str):
                entries = kp_ranked.get(cat, [])
                terms  = [e["term"] for e in entries]
                counts = {e["term"]: int(round(e["freq"] * 100)) for e in entries}
                return terms, counts

            h_terms, h_counts = _kp_lists_and_counts("hormones")
            m_terms, m_counts = _kp_lists_and_counts("metabolites")
            p_terms, p_counts = _kp_lists_and_counts("proteins")

            n_cooccur_m = b.get("n_cooccur_metabolic", 0)
            n_cooccur_h = b.get("n_cooccur_hormonal", 0)
            layer_parts = []
            if n_cooccur_m:
                layer_parts.append(f"metabolic: {n_cooccur_m}")
            if n_cooccur_h:
                layer_parts.append(f"hormonal: {n_cooccur_h}")
            layer_str = "  |  Layers — " + ", ".join(layer_parts) if layer_parts else ""

            qm1 = search_results.get(o1, {}).get("query_metabolic", "")
            qh1 = search_results.get(o1, {}).get("query_hormonal", "")
            description = (
                f"{cfg.VIZ_LABEL} general connection — {o1} ↔ {o2}\n"
                f"Otsuka–Ochiai coefficient: {b['mean']:.5f} ± {b['std']:.5f}  "
                f"(cv={b['cv']:.3f},  {cfg.N_BOOTSTRAP} iterations × "
                f"{int(cfg.SAMPLE_FRACTION*100)}% per organ,  "
                f"raw cross-mentions/iter: {b.get('raw_mean', 0):.1f} ± {b.get('raw_std', 0):.1f})\n"
                f"Papers in pools: {b.get('n_found_o1',0)} ({o1}) + "
                f"{b.get('n_found_o2',0)} ({o2})  |  "
                f"Cross-mentioning: {b['n_cooccur_total']}"
                f"{layer_str}"
            )

            pair     = (min(o1, o2), max(o1, o2))
            ref_data = ref_edges.get(pair)

            # LLM-generated summary (run_llm_descriptions.py), grounded in
            # this pair's same-sentence cross-mention evidence papers. Empty
            # string (no-op) if that script hasn't been run yet. When present,
            # its citation numbers [1]-[5] index into its OWN selected_papers
            # list (llm_pair_descriptions_*.json), not the full bootstrap
            # evidence list — "papers" below must match whichever text is
            # actually shown (ai_description takes priority in the sidebar),
            # or citation links resolve to the wrong paper / nothing at all.
            llm_entry           = pair_llm.get(edge_key, {})
            llm_description     = llm_entry.get("description", "")
            citation_papers     = llm_entry.get("papers") if llm_description else None

            # LLM-classified connection type(s) (run_llm_descriptions.py),
            # 1-3 categories from cfg.CONNECTION_TYPES via 3-run majority
            # vote — see Literature_Search/llm_connection_type.py. Empty
            # list (no badge) if that script hasn't been run yet.
            type_labels       = get_connection_types(llm_types, o1, o2, cfg.CONNECTION_TYPES)
            conn_type         = type_labels[0] if type_labels else ""
            conn_type_others  = type_labels[1:]

            G.add_edge(o1, o2)
            G.edges[o1, o2]['color']       = "#0ea5e9"
            G.edges[o1, o2]['description'] = description
            G.edges[o1, o2]['merged_data'] = {
                "n_papers_found":               b["n_found"],
                "n_papers_cooccur":             b["n_cooccur_total"],
                "n_papers_cooccur_metabolic":   n_cooccur_m,
                "n_papers_cooccur_hormonal":    n_cooccur_h,
                "papers":                       citation_papers if citation_papers else papers,
                "key_players_hormones":         h_terms,
                "key_players_metabolites":      m_terms,
                "key_players_proteins":         p_terms,
                "key_players_counts_hormones":    h_counts,
                "key_players_counts_metabolites": m_counts,
                "key_players_counts_proteins":    p_counts,
                "key_players_bootstrap":        True,
                "connection_type":        conn_type,
                "connection_type_others": conn_type_others,
                "pubmed_query":          qm1,
                "pubmed_query_hormonal": qh1,
                "ai_description": llm_description,
                # Adaptive-threshold / comparison-toggle fields — see
                # threshold_control / comparison_toggle_label in
                # export_network_to_cytoscape_dashboard's docstring. No-ops
                # (edge just never appears "only in reference") if
                # comparison_data isn't available yet.
                "bootstrap_mean": b["mean"],
                "is_ref_edge":    ref_data is not None,
                "ref_n_papers":   ref_data.get("n_papers_found", 0) if ref_data else 0,
            }

        # Elbow suggestion recomputed here (not just taken from
        # cfg.MIN_BOOTSTRAP_MEAN) so the slider's reset target is always the
        # true elbow even if the config pins a literal number instead of
        # Elbow.
        sorted_means = sorted((v["mean"] for v in boot.values() if v["mean"] > 0), reverse=True)
        elbow_mean   = round(kneedle_elbow(sorted_means), 6) if sorted_means else 0.0

        comparison_toggle_label = (
            f"Comparison vs {cfg.VIZ_LABEL} Reference" if comparison_data else None
        )
        comparison_tabs = (
            _cmp._info_tabs(categories, cfg, comparison_data["cohort_csv_name"])
            if comparison_data else []
        )

        export_network_to_cytoscape_dashboard(
            graph=G,
            filename=str(output_path),
            include_legend=False,
            title=cfg.VIZ_TITLE_ROBUST,
            info_panel_tabs=_robust_info_tabs(cfg, organs, search_results, boot) + comparison_tabs,
            extra_topbar_buttons=[
                {
                    "icon": "📊",
                    "label": "Bootstrap Overview",
                    "html": overview_html_str,
                },
            ],
            comparison_toggle_label=comparison_toggle_label,
            comparison_lit_label=lit_label,
            comparison_ref_label="Reference Metabolic Network",
            threshold_control={
                "default": cfg.MIN_BOOTSTRAP_MEAN,
                "elbow":   elbow_mean,
                "label":   "Otsuka–Ochiai threshold",
            },
            total_possible_edges=len(boot),
        )
        print(f"[ok] Robust network saved: {output_path}  "
              f"({n_robust}/{len(boot)} robust edges at threshold {cfg.MIN_BOOTSTRAP_MEAN}"
              + (f", {len(ref_edges)} reference edges" if comparison_data else "") + ")")
    except Exception as e:
        import traceback
        print(f"[!] Robust network visualization failed: {e}")
        traceback.print_exc()


def _literature_stats_section(organs: list[str], search_results: dict, boot: dict) -> str:
    """
    "Literature Statistics" block for the bottom of the Search Filters tab:
    per-organ and per-connection paper counts, plus a canvas-drawn Sankey
    diagram (ribbon width = paper count for that connection), in the style
    of Fig. 1 in "Organ cross-talk: molecular mechanisms, biological
    functions, and therapeutic interventions for diseases" (Che et al.,
    2026) — "the thickness of the connecting chords represents the number
    of articles related to each organ pair."
    """
    from Visualisation.networkBuilderUtils import ORGAN_COLORS, DEFAULT_NODE_COLOR

    organ_papers = {o: search_results.get(o, {}).get("n_found", 0) for o in organs}

    links = []
    for b in boot.values():
        o1, o2 = b.get("organ1"), b.get("organ2")
        n = b.get("n_cooccur_total", 0)
        if o1 in organ_papers and o2 in organ_papers and n > 0:
            # Consistent left/right assignment: the organ with the larger
            # total paper pool is the "source" (left column).
            if organ_papers[o1] < organ_papers[o2]:
                o1, o2 = o2, o1
            links.append({"source": o1, "target": o2, "value": n})

    organ_rows = "".join(
        f'<div style="display:flex;justify-content:space-between;padding:3px 0;'
        f'border-bottom:1px solid #1e293b;font-size:0.78rem">'
        f'<span style="color:#e2e8f0">{o}</span>'
        f'<span style="color:#64748b">{organ_papers[o]:,} papers</span></div>'
        for o in sorted(organs, key=lambda o: -organ_papers[o])
    )
    link_rows = "".join(
        f'<div style="display:flex;justify-content:space-between;padding:3px 0;'
        f'border-bottom:1px solid #1e293b;font-size:0.78rem">'
        f'<span style="color:#e2e8f0">{l["source"]} ↔ {l["target"]}</span>'
        f'<span style="color:#64748b">{l["value"]:,} papers</span></div>'
        for l in sorted(links, key=lambda l: -l["value"])
    )

    sankey_organs = [o for o in sorted(organs, key=lambda o: -organ_papers[o])
                     if any(l["source"] == o or l["target"] == o for l in links)]
    sankey_data = {
        "organs":      sankey_organs,
        "organColors": {o: ORGAN_COLORS.get(o, DEFAULT_NODE_COLOR) for o in sankey_organs},
        "links":       links,
    }
    sankey_json = json.dumps(sankey_data, ensure_ascii=False)

    return f"""
    <div class="info-h2">Literature Statistics</div>
    <div class="info-h2">Organ Cross-Talk Sankey</div>
    <p class="info-p" style="color:#94a3b8">
      Ribbon thickness = number of same-sentence co-occurring papers for that
      organ pair.
    </p>
    <div style="overflow-x:auto">
      <canvas id="lit-sankey" style="display:block"></canvas>
    </div>
    <script>
    (function() {{
      const D = {sankey_json};
      const canvas = document.getElementById('lit-sankey');
      if (!D.links.length) {{ canvas.style.display = 'none'; return; }}
      const ctx = canvas.getContext('2d');

      // Circular chord diagram: organs arranged around a ring (arc length
      // proportional to total evidence touching that organ), links drawn as
      // ribbons curving through the center — ribbon width = paper count.
      const organs = D.organs, links = D.links;
      const SIZE = 640, CX = SIZE / 2, CY = SIZE / 2, R = 210, ARC_W = 14, LABEL_GAP = 12;
      const GAP_DEG = organs.length > 1 ? 1.4 : 0;

      canvas.width = SIZE; canvas.height = SIZE;
      ctx.clearRect(0, 0, SIZE, SIZE);

      const totals = {{}};
      organs.forEach(o => totals[o] = 0);
      links.forEach(l => {{ totals[l.source] += l.value; totals[l.target] += l.value; }});
      const grand = organs.reduce((s, o) => s + totals[o], 0) || 1;
      const usableDeg = 360 - GAP_DEG * organs.length;

      const toRad = d => d * Math.PI / 180;
      const pointOnCircle = (deg, radius) => {{
        const rad = toRad(deg);
        return [CX + radius * Math.cos(rad), CY + radius * Math.sin(rad)];
      }};
      const colorFor = o => D.organColors[o] || '#2563eb';

      // Assign each organ an angular range around the ring, starting at top.
      let angle = -90;
      const organRange = {{}};
      organs.forEach(o => {{
        const span = (totals[o] / grand) * usableDeg;
        organRange[o] = {{ start: angle, end: angle + span }};
        angle += span + GAP_DEG;
      }});

      // Sub-divide each organ's arc per link (deterministic order) so each
      // ribbon gets an exact angular range on both ends.
      const sortedLinks = [...links].sort((a, b) => (a.source + a.target).localeCompare(b.source + b.target));
      const cursor = {{}};
      organs.forEach(o => cursor[o] = organRange[o].start);

      const ribbonPaths = [];
      sortedLinks.forEach(l => {{
        const span = (l.value / grand) * usableDeg;
        const a0 = cursor[l.source]; cursor[l.source] += span;
        const b0 = cursor[l.target]; cursor[l.target] += span;
        const a1 = a0 + span, b1 = b0 + span;

        const [x0, y0]   = pointOnCircle(a0, R);
        const [x0b, y0b] = pointOnCircle(a1, R);
        const [x1, y1]   = pointOnCircle(b0, R);
        const [x1b, y1b] = pointOnCircle(b1, R);

        const path = new Path2D();
        path.moveTo(x0, y0);
        path.quadraticCurveTo(CX, CY, x1, y1);
        path.lineTo(x1b, y1b);
        path.quadraticCurveTo(CX, CY, x0b, y0b);
        path.closePath();

        const grad = ctx.createLinearGradient(x0, y0, x1, y1);
        grad.addColorStop(0, colorFor(l.source) + '99');
        grad.addColorStop(1, colorFor(l.target) + '55');
        ctx.fillStyle = grad;
        ctx.fill(path);

        ribbonPaths.push({{ path, link: l }});
      }});

      // Organ arcs + radial labels
      organs.forEach(o => {{
        const {{ start, end }} = organRange[o];
        if (end <= start) return;
        ctx.beginPath();
        ctx.arc(CX, CY, R, toRad(start), toRad(end));
        ctx.lineWidth = ARC_W;
        ctx.strokeStyle = colorFor(o);
        ctx.stroke();

        const mid = (start + end) / 2;
        const [lx, ly] = pointOnCircle(mid, R + ARC_W / 2 + LABEL_GAP);
        const flip = mid > 90 && mid < 270;
        ctx.save();
        ctx.translate(lx, ly);
        ctx.rotate(toRad(mid) + (flip ? Math.PI : 0));
        ctx.fillStyle = '#e2e8f0';
        ctx.font = '11px system-ui';
        ctx.textAlign = flip ? 'right' : 'left';
        ctx.textBaseline = 'middle';
        ctx.fillText(o, 0, 0);
        ctx.restore();
      }});

      const tip = document.createElement('div');
      tip.style.cssText = 'position:fixed;pointer-events:none;display:none;background:#1e293b;' +
        'border:1px solid #475569;border-radius:6px;padding:6px 10px;font-size:0.76rem;' +
        'color:#e2e8f0;z-index:999';
      document.body.appendChild(tip);

      canvas.addEventListener('mousemove', e => {{
        const r = canvas.getBoundingClientRect();
        const x = (e.clientX - r.left) * (canvas.width / r.width);
        const y = (e.clientY - r.top) * (canvas.height / r.height);
        let hit = null;
        for (const rp of ribbonPaths) {{
          if (ctx.isPointInPath(rp.path, x, y)) {{ hit = rp.link; break; }}
        }}
        if (hit) {{
          tip.innerHTML = `<strong>${{hit.source}} ↔ ${{hit.target}}</strong><br>${{hit.value.toLocaleString()}} co-occurring papers`;
          tip.style.display = 'block';
          tip.style.left = (e.clientX + 14) + 'px';
          tip.style.top = (e.clientY - 10) + 'px';
        }} else {{
          tip.style.display = 'none';
        }}
      }});
      canvas.addEventListener('mouseleave', () => {{ tip.style.display = 'none'; }});
    }})();
    </script>
    <p class="info-p">

    </p>
    <div style="display:flex;gap:24px;flex-wrap:wrap;margin-bottom:18px">
      <div style="flex:1;min-width:220px">
        <div style="font-size:0.74rem;color:#94a3b8;margin-bottom:6px;text-transform:uppercase;letter-spacing:0.03em">Papers per organ</div>
        {organ_rows}
      </div>
      <div style="flex:1;min-width:220px">
        <div style="font-size:0.74rem;color:#94a3b8;margin-bottom:6px;text-transform:uppercase;letter-spacing:0.03em">Papers per connection</div>
        {link_rows if link_rows else '<p style="color:#64748b;font-size:0.8rem">No cross-mention evidence yet.</p>'}
      </div>
    </div>
"""


def _connection_types_section(cfg) -> str:
    """
    "Connection Types" tab: explains the LLM classification step and lists
    the current CONNECTION_TYPES categories (kept identical to
    reference_network_only_metabolic/config.py so a label means the same
    thing across all five dashboards).
    """
    ct_rows = "".join(
        f'<tr><td style="padding:3px 8px 3px 0;color:#94a3b8;white-space:nowrap">'
        f'<strong style="color:#e2e8f0">{v["label"]}</strong></td>'
        f'<td style="padding:3px 0;color:#94a3b8">{v["description"]}</td></tr>'
        for v in cfg.CONNECTION_TYPES.values()
    )
    return f"""
    <div class="info-h2">LLM Classification</div>
    <p class="info-p">
      For each robust organ pair an LLM (<code>llama3.2</code>, via Ollama) reads up to
      <strong>{cfg.LLM_MAX_PAPERS} of this pair's same-sentence cross-mention papers</strong>
      (the same bootstrap evidence used for the AI summary) and assigns
      <strong>1 to 3 connection types</strong> from the categories below, ranked
      from most to least feasible — a pair only gets a second or third type if it is
      independently well-supported by the papers.
    </p>
    <p class="info-p">
      To reduce single-sample noise, the classification is run
      <strong>3 times independently</strong> per pair; only types that appear in at
      least 2 of the 3 runs are kept (majority vote). If the 3 runs disagree on
      everything, the single most-voted type is kept, so every classified pair
      always ends up with at least one type.
    </p>
    <p class="info-p" style="color:#94a3b8">
      Generated by <code>run_llm_descriptions.py</code> — run it (or add
      <code>--reset</code> to reclassify after editing CONNECTION_TYPES in this
      condition's config) then rebuild with <code>--viz-only</code>.
    </p>
    <div class="info-h2">Connection Type Categories</div>
    <table style="width:100%;border-collapse:collapse;font-size:12px">
      {ct_rows}
    </table>
    <p class="info-p" style="margin-top:10px">
      All assigned types are shown as coloured badges on each edge in the network.
    </p>
"""


def _robust_info_tabs(cfg, organs: list[str], search_results: dict, boot: dict) -> list[dict]:
    metabolic_kws = "\n".join(f"  {kw}" for kw in cfg.METABOLIC_KEYWORDS)
    hormonal_kws  = "\n".join(f"  {kw}" for kw in cfg.HORMONAL_KEYWORDS)
    condition_kws = "\n".join(f"  {kw}" for kw in cfg.CONDITION_KEYWORDS)
    crosstalk_kws = "\n".join(f"  {kw}" for kw in cfg.CROSSTALK_KEYWORDS)
    return [
        {
            "id": "bootstrap",
            "label": "Bootstrap Analysis",
            "content": f"""
    <div class="info-h2">Condition: {cfg.VIZ_LABEL}</div>
    <p class="info-p">
      Organ pairs with mean Otsuka–Ochiai coefficient ≥ <strong>{cfg.MIN_BOOTSTRAP_MEAN}</strong>
      per iteration in the <strong>{cfg.VIZ_LABEL.lower()}</strong> literature context.
      Papers are fetched per organ (two queries per organ: metabolic + hormonal layer)
      and cross-organ mentions are counted during bootstrapping.
    </p>
    <div class="info-h2">Per-Organ Query (run twice per organ)</div>
    <div class="info-code">(MeSH_ORGAN OR aliases_ORGAN)
AND [METABOLIC_FILTER | HORMONAL_FILTER]
AND CONDITION_FILTER ({cfg.VIZ_LABEL})</div>
    <p class="info-p">
      In each bootstrap iteration {int(cfg.SAMPLE_FRACTION*100)}% of each organ's paper
        pool is sampled and cross-organ mentions are counted. A paper only counts as a
        cross-mention if the two organs' names (or aliases) appear together in the
        <strong>same sentence</strong> of its title/abstract. The resulting mean count is then normalized into an
        <strong>Otsuka–Ochiai coefficient</strong>:
    </p>
    <div style="display:flex;align-items:center;gap:10px;margin:14px 0;padding:12px 16px;
                background:#0f172a;border-radius:8px;border:1px solid #334155;font-size:1.05rem">
      <span style="font-style:italic;color:#38bdf8">OO</span><span style="color:#e2e8f0">(A, B)&nbsp;=</span>
      <span style="display:inline-flex;flex-direction:column;align-items:center;line-height:1.3;color:#e2e8f0">
        <span style="padding:0 6px 4px;border-bottom:1.5px solid #94a3b8">mean cross-mentions</span>
        <span style="padding:4px 6px 0">√(n<sub>A</sub> × n<sub>B</sub>)</span>
      </span>
    </div>
    <p class="info-p">
      so an organ with a much larger literature pool doesn't dominate edges just
      by volume.
    </p>
    <div class="info-h2">Bootstrap Parameters</div>
    <div class="info-stat-grid">
      <div class="info-stat"><div class="info-stat-val">{cfg.N_BOOTSTRAP}</div><div class="info-stat-lbl">Iterations</div></div>
      <div class="info-stat"><div class="info-stat-val">{int(cfg.SAMPLE_FRACTION*100)}%</div><div class="info-stat-lbl">Sample fraction</div></div>
      <div class="info-stat"><div class="info-stat-val">{cfg.MAX_PAPERS:,}</div><div class="info-stat-lbl">Max papers/organ/layer</div></div>
      <div class="info-stat"><div class="info-stat-val">{cfg.MIN_BOOTSTRAP_MEAN}</div><div class="info-stat-lbl">Robust threshold (OO coeff.)</div></div>
      <div class="info-stat"><div class="info-stat-val" id="stat-val-edges">—</div><div class="info-stat-lbl">Robust edges</div></div>
    </div>
    <div class="info-h2">Key Player Scores</div>
    <p class="info-p">
      Each edge shows <strong>metabolites</strong>, <strong>hormones</strong>, and
      <strong>proteins</strong> with their bootstrap frequency — the fraction of
      the {cfg.N_BOOTSTRAP} iterations in which that molecule was mentioned in
      the sampled co-occurring papers. A score of 90% means the molecule appeared
      in 9 out of 10 random samples, indicating consistent literature support.
    </p>
    <p class="info-p">
      <strong>CV (coefficient of variation)</strong>: std / mean of the Otsuka–Ochiai
      coefficient across bootstrap samples. Lower CV = evidence is spread across many
      papers rather than a few outlier samples.
    </p>
    <p class="info-p">
      See <strong>bootstrap_overview_{cfg.CONDITION_NAME}.html</strong> for the
      heatmap and full strength ranking.
    </p>
""",
        },
        {
            "id": "filters",
            "label": "Search Filters",
            "content": f"""
    <div class="info-h2">How papers are selected</div>
    <p class="info-p">
      Two PubMed queries per organ — one for the <strong>metabolic layer</strong>
      and one for the <strong>hormonal layer</strong> — each ANDing the organ terms
      with its layer filter, the condition filter, and the crosstalk filter
      below, fetching up to {cfg.MAX_PAPERS:,} papers each. Papers from both
      queries are merged per organ; each paper is tagged with the layer(s) it
      was retrieved from. During bootstrapping, papers from one organ's pool
      that mention another organ in title/abstract count as cross-organ
      evidence and drive the stacked bar chart layer breakdown. The mean cross-mention
      count is then normalized into an Otsuka–Ochiai coefficient using each organ's
      total literature pool size (n_found) so heavily-studied organs don't dominate
      purely by volume.
    </p>
    <p class="info-p">
      Organ and crosstalk terms are matched as <strong>whole words</strong>,
      not substrings — "renal" no longer matches inside "adrenal", so a
      paper about the adrenal gland can't be mistaken for kidney evidence.
      A match is also discarded if it's <strong>negated</strong> — "non-renal",
      "nonrenal", "non renal", and "not renal" don't count as a positive
      mention of that organ or crosstalk term.
    </p>
    <div class="info-h2">METABOLIC_FILTER</div>
    <p class="info-p">
      At least one of these keywords must appear in title or abstract.
    </p>
    <div class="info-code">{metabolic_kws}</div>
    <div class="info-h2">HORMONAL_FILTER</div>
    <p class="info-p">
      At least one of these keywords must appear in title or abstract.
    </p>
    <div class="info-code">{hormonal_kws}</div>
    <div class="info-h2">CONDITION_FILTER — {cfg.VIZ_LABEL}</div>
    <p class="info-p">
      Required by both queries — restricts results to the
      <strong>{cfg.VIZ_LABEL.lower()}</strong> physiological context.
    </p>
    <div class="info-code">{condition_kws}</div>
    <div class="info-h2">CROSSTALK_FILTER</div>
    <p class="info-p">
      Required by both queries too — at least one of these keywords must also
      appear somewhere in the document (a coarse pre-filter), restricting
      results to papers actually discussing inter-organ relationships, not
      just papers that happen to mention an organ and a metabolic/hormonal
      term without relating them to anything else. After fetching, papers are
      filtered further: kept only if the organ's own name/alias and one of
      these crosstalk terms appear together in the <strong>same sentence</strong>
      of the title/abstract, not merely somewhere in the same document.
    </p>
    <div class="info-code">{crosstalk_kws}</div>
""",
        },
        {
            "id": "literature_stats",
            "label": "Literature Statistics",
            "content": _literature_stats_section(organs, search_results, boot),
        },
        {
            "id": "connection_types",
            "label": "Connection Types",
            "content": _connection_types_section(cfg),
        },
    ]


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--condition", required=True,
                        choices=list(CONDITION_CONFIGS),
                        help="Which condition config to use.")
    parser.add_argument("--reset",    action="store_true",
                        help="Delete search + bootstrap cache and restart.")
    parser.add_argument("--viz-only", action="store_true",
                        help="Skip search + bootstrap; rebuild HTML only.")
    args = parser.parse_args()

    cfg = _load_config(args.condition)

    from Literature_Search.pubmed_search import ORGAN_ALIASES, ORGAN_MESH

    cohort_csv = COHORT_CONNECTIONS_CSV[cfg.CONDITION_NAME]
    if not cohort_csv.exists():
        print(f"[!] Cohort connections file not found: {cohort_csv}")
        sys.exit(1)
    cohort_organs = load_cohort_organs(cohort_csv)
    organs  = sorted(o for o in cohort_organs if o in ORGAN_ALIASES)
    skipped = sorted(set(cohort_organs) - set(ORGAN_ALIASES))
    if skipped:
        print(f"[!] Cohort organs not in ORGAN_ALIASES, skipping: {skipped}")

    metabolic_filter = _keyword_clause(cfg.METABOLIC_KEYWORDS)
    hormonal_filter  = _keyword_clause(cfg.HORMONAL_KEYWORDS)
    condition_filter = _keyword_clause(cfg.CONDITION_KEYWORDS)
    crosstalk_filter = _keyword_clause(cfg.CROSSTALK_KEYWORDS)

    out_dir       = HERE / cfg.CONDITION_NAME
    out_dir.mkdir(parents=True, exist_ok=True)
    search_json   = out_dir / f"search_results_{cfg.CONDITION_NAME}.json"
    boot_json     = out_dir / f"bootstrap_results_{cfg.CONDITION_NAME}.json"
    robust_html   = out_dir / f"robust_network_{cfg.CONDITION_NAME}.html"

    print(f"\n[i] Condition  : {cfg.VIZ_LABEL}")
    print(f"[i] {len(organs)} organs (from {cohort_csv.name})")
    print(f"[i] Max papers : {cfg.MAX_PAPERS:,}/organ/layer | Years back: {cfg.YEARS_BACK}")
    print(f"[i] Layers     : metabolic ({len(cfg.METABOLIC_KEYWORDS)} kw) + "
          f"hormonal ({len(cfg.HORMONAL_KEYWORDS)} kw)")
    print(f"[i] Bootstrap  : {cfg.N_BOOTSTRAP} × {int(cfg.SAMPLE_FRACTION*100)}% per organ")
    print(f"[i] Threshold  : mean Otsuka–Ochiai coefficient ≥ {cfg.MIN_BOOTSTRAP_MEAN}\n")

    if args.reset:
        for p in (search_json, boot_json):
            if p.exists():
                p.unlink()
                print(f"[i] Deleted cache: {p.name}")

    # ── Step 1: Per-organ PubMed search ──────────────────────────────────────
    if not args.viz_only:
        search_results = run_search(
            organs, search_json, ORGAN_ALIASES, ORGAN_MESH,
            metabolic_filter, hormonal_filter, condition_filter,
            crosstalk_filter, cfg.CROSSTALK_KEYWORDS, cfg,
            resume=not args.reset,
        )
    else:
        if not search_json.exists():
            print(f"[!] No search results at {search_json}. Run without --viz-only first.")
            sys.exit(1)
        with open(search_json, encoding="utf-8") as f:
            search_results = json.load(f)
        print(f"[i] Loaded {len(search_results)} organs from cache.")

    # Restrict to the current cohort's organs — cached search results may
    # still hold entries for organs outside this condition's cohort CSV
    # (e.g. from before this filter existed, or from a broader prior run).
    n_before = len(search_results)
    search_results = {o: v for o, v in search_results.items() if o in organs}
    if len(search_results) != n_before:
        print(f"[i] Filtered search cache to cohort organs: "
              f"{n_before} -> {len(search_results)}")

    # ── Step 2: Bootstrap ─────────────────────────────────────────────────────
    if not args.viz_only:
        print(f"\n[i] Running bootstrap ({cfg.N_BOOTSTRAP} iterations × "
              f"{int(cfg.SAMPLE_FRACTION*100)}% per organ)...")
        boot_results = run_bootstrap(search_results, boot_json, cfg, ORGAN_ALIASES,
                                     resume=not args.reset)
    else:
        if not boot_json.exists():
            print(f"[!] No bootstrap results at {boot_json}. Run without --viz-only first.")
            sys.exit(1)
        with open(boot_json, encoding="utf-8") as f:
            boot_results = json.load(f)
        print(f"[i] Loaded {len(boot_results)} cached bootstrap results.")

    # Restrict to pairs where both organs are in the current cohort — the
    # bootstrap cache may still hold pairs outside this condition's cohort.
    n_before = len(boot_results)
    boot_results = {k: v for k, v in boot_results.items()
                    if v.get("organ1") in organs and v.get("organ2") in organs}
    if len(boot_results) != n_before:
        print(f"[i] Filtered bootstrap cache to cohort organ pairs: "
              f"{n_before} -> {len(boot_results)}")

    # Resolve MIN_BOOTSTRAP_MEAN = Elbow (if used) now that the cohort's
    # bootstrap means are known — every downstream use of cfg.MIN_BOOTSTRAP_MEAN
    # (overview HTML, robust network, printed summaries) sees the real number.
    if is_elbow(cfg.MIN_BOOTSTRAP_MEAN):
        means = [v.get("mean", 0) for v in boot_results.values()]
        cfg.MIN_BOOTSTRAP_MEAN = resolve_min_bootstrap_mean(cfg.MIN_BOOTSTRAP_MEAN, means)
        print(f"[i] Resolved 'Elbow' threshold -> mean ≥ {cfg.MIN_BOOTSTRAP_MEAN}")

    # ── Step 3: Visualizations ────────────────────────────────────────────────
    # Overview (heatmap/bar chart) and the comparison-vs-reference view are
    # both built as strings and embedded into robust_network_*.html itself —
    # no separate bootstrap_overview_*.html or comparison_*_vs_reference.html
    # files are written. One file per condition.
    print(f"\n[i] Building overview HTML...")
    overview_html_str = _build_overview_html_string(boot_results, cfg)

    print(f"[i] Loading comparison data...")
    comparison_data = _cmp.compute_comparison_data(cfg.CONDITION_NAME, cfg)
    if comparison_data is None:
        print(f"[i] Note: reference or bootstrap data not ready yet — the "
              f"Comparison toggle will be unavailable until you re-run "
              f"this with --viz-only after both are available.")

    print(f"[i] Building robust network...")
    build_robust_network(boot_results, search_results, robust_html, cfg, organs,
                         overview_html_str, comparison_data)

    print(f"\n[ok] All outputs in {cfg.CONDITION_NAME}/")


if __name__ == "__main__":
    main()
