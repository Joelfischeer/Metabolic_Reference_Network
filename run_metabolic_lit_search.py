"""
run_metabolic_lit_search.py
===========================
Metabolic-focused PubMed search for all organ-organ pairs.

Query per pair (single strategy, no cascade):
    (MeSH_A OR aliases_A) AND (MeSH_B OR aliases_B)
    AND <METABOLIC_FILTER>

where METABOLIC_FILTER is the curated metabolic keyword clause below.

Results saved to:  metabolic_data/metabolic_literature_results.json
Visualization at:  metabolic_data/metabolic_literature_network.html

Run:
    python run_metabolic_lit_search.py              # first run or resume
    python run_metabolic_lit_search.py --reset      # wipe cache and restart
    python run_metabolic_lit_search.py --force-empty  # re-search 0-paper edges
    python run_metabolic_lit_search.py --viz-only   # skip search, rebuild viz
"""

import sys
import json
import time
import re
import argparse
from pathlib import Path
from datetime import datetime, timedelta

import requests

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))

# Windows console/redirected-output encoding defaults to cp1252, which can't
# encode the arrows/ellipses used in progress prints below — force UTF-8 so
# the run doesn't crash mid-way through (e.g. when stdout is piped to a file).
for _stream in (sys.stdout, sys.stderr):
    if hasattr(_stream, "reconfigure"):
        _stream.reconfigure(encoding="utf-8", errors="replace")

from Literature_Search.pubmed_search import (
    ORGAN_ALIASES,
    ORGAN_MESH,
    fetch_abstracts,
    extract_key_players,
    _merge_synonyms,
    HORMONE_SYNONYMS,
    METABOLITE_SYNONYMS,
    PROTEIN_SYNONYMS,
)

# ── Load configuration from the output folder ─────────────────────────────────
# Edit:  reference_network_only_metabolic/config.py
OUTPUT_DIR = HERE / "reference_network_only_metabolic"
sys.path.insert(0, str(OUTPUT_DIR))
import config as _cfg

MAX_PAPERS         = _cfg.MAX_PAPERS
YEARS_BACK         = _cfg.YEARS_BACK
DELAY              = _cfg.DELAY
METABOLIC_KEYWORDS = _cfg.METABOLIC_KEYWORDS
CROSSTALK_KEYWORDS = _cfg.CROSSTALK_KEYWORDS
CONNECTION_TYPES   = _cfg.CONNECTION_TYPES
LLM_MODEL          = _cfg.LLM_MODEL
LLM_MAX_PAPERS     = _cfg.LLM_MAX_PAPERS
VIZ_TITLE          = _cfg.VIZ_TITLE

# Condition keyword filters -- see config.py for the exact lists and why
# these are applied locally (post-extraction) instead of as a PubMed
# query-time filter like Edge_cosine_met/general_reference_network.
CONDITION_KEYWORDS = {
    "healthy": _cfg.CONDITION_KEYWORDS_HEALTHY,
    "obese":   _cfg.CONDITION_KEYWORDS_OBESE,
}

OUTPUT_JSON         = OUTPUT_DIR / "metabolic_literature_results.json"
OUTPUT_HTML         = OUTPUT_DIR / "metabolic_literature_network.html"
OUTPUT_LLM_TYPES    = OUTPUT_DIR / "metabolic_connection_types.json"
# Search scope: every pair across the full 14-organ set (91 pairs). The
# healthy cohort CSV is no longer a search-time filter -- it's read
# separately in build_viz() to flag which pairs/organs belong to the
# "Healthy" display toggle in the dashboard.
EDGE_FILTER_CSV     = OUTPUT_DIR / "all_organ_connections.csv"
HEALTHY_FILTER_CSV  = OUTPUT_DIR / "healthy_cohort_connections.csv"
OBESE_FILTER_CSV    = OUTPUT_DIR / "obese_cohort_connections.csv"
LLM_DESCRIPTIONS    = OUTPUT_DIR / "metabolic_llm_descriptions.json"
ORGAN_DESCRIPTIONS  = HERE / "metabolic_data" / "organ_descriptions.json"

# Per-condition cohort CSV + output files -- kept fully separate from the
# "all" outputs above (OUTPUT_JSON/OUTPUT_LLM_TYPES/LLM_DESCRIPTIONS), which
# stay untouched by --condition runs.
CONDITION_COHORT_CSV = {
    "healthy": HEALTHY_FILTER_CSV,
    "obese":   OBESE_FILTER_CSV,
}


def condition_output_paths(condition: str) -> dict[str, Path]:
    return {
        "results":     OUTPUT_DIR / f"metabolic_literature_results_{condition}.json",
        "llm_types":   OUTPUT_DIR / f"metabolic_connection_types_{condition}.json",
        "llm_descs":   OUTPUT_DIR / f"metabolic_llm_descriptions_{condition}.json",
    }

def load_edge_filter(csv_path: Path) -> list[tuple[str, str]]:
    """
    Parse a 0/1 adjacency matrix CSV.  Returns sorted canonical pairs where
    the cell value is 1.  Both upper and lower triangle are accepted; diagonal
    and empty cells are ignored.
    """
    import csv as _csv
    pairs: set[tuple[str, str]] = set()
    with open(csv_path, encoding="utf-8") as f:
        reader = _csv.reader(f)
        header = next(reader)
        col_organs = header[1:]          # first column is the row-label
        for row in reader:
            if not row:
                continue
            row_organ = row[0].strip()
            for col_idx, val in enumerate(row[1:], start=0):
                val = val.strip()
                if val == "1" and col_idx < len(col_organs):
                    col_organ = col_organs[col_idx].strip()
                    if row_organ and col_organ and row_organ != col_organ:
                        pairs.add((min(row_organ, col_organ),
                                   max(row_organ, col_organ)))
    return sorted(pairs)


def _build_metabolic_filter(keywords: list[str]) -> str:
    """Convert the keyword list into a PubMed boolean clause."""
    parts = []
    for kw in keywords:
        # Multi-word phrases need quoting; single words do not
        if " " in kw or "-" in kw:
            parts.append(f'"{kw}"[Title/Abstract]')
        else:
            parts.append(f'{kw}[Title/Abstract]')
    return "(" + " OR ".join(parts) + ")"

METABOLIC_FILTER = _build_metabolic_filter(METABOLIC_KEYWORDS)
CROSSTALK_FILTER = _build_metabolic_filter(CROSSTALK_KEYWORDS)

NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


# ── HTTP helper ───────────────────────────────────────────────────────────────

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
            print(f"    [!] HTTP error (attempt {attempt+1}/{retries}): {e} — retrying in {wait}s")
            time.sleep(wait)
    return None


# ── Query builder ─────────────────────────────────────────────────────────────

def _alias_clause(organ: str) -> str:
    """Build (MeSH OR alias1 OR alias2 ...) clause for one organ."""
    mesh    = ORGAN_MESH.get(organ, f'"{organ}"[MeSH Terms]')
    aliases = ORGAN_ALIASES.get(organ, [organ])
    alias_terms = " OR ".join(f'"{a}"[Title/Abstract]' for a in aliases)
    return f"({mesh} OR {alias_terms})"


def build_query(organ1: str, organ2: str) -> str:
    """Single metabolic-focused query for a pair of organs."""
    c1 = _alias_clause(organ1)
    c2 = _alias_clause(organ2)
    return f"{c1} AND {c2} AND {METABOLIC_FILTER} AND {CROSSTALK_FILTER}"


# ── Word-boundary / negation-safe term matching ────────────────────────────
# Plain substring checks (`"renal" in text`) false-positive on words that
# merely contain the term, e.g. "renal" inside "adrenal" — a paper about the
# adrenal gland would otherwise count as kidney evidence. \b...\b anchors the
# match to whole-word boundaries so this can't happen. A leading negative
# lookbehind additionally excludes negated mentions ("non-renal", "nonrenal",
# "non renal", "not renal") from counting as a positive mention of that
# organ/keyword — these explicitly say the opposite.

def _word_boundary_pattern(term: str):
    escaped = re.escape(term.lower())
    return re.compile(r"(?<!non-)(?<!non )(?<!non)(?<!not )\b" + escaped + r"\b")


def _compile_patterns(terms) -> list:
    return [_word_boundary_pattern(t) for t in terms]


def _any_match(patterns: list, text: str) -> bool:
    return any(p.search(text) for p in patterns)


# ── Multi-organ same-sentence disambiguation ────────────────────────────────
# A qualifying sentence (organ1 + organ2 + crosstalk keyword co-present) can
# still be a false positive for the pair (organ1, organ2) specifically when a
# THIRD organ is also named in that sentence, e.g. "the pancreas modulates
# ... the heart and colon" -- Pancreas is the true agent, Heart/Colon are
# just co-listed targets, not a Heart<->Colon relationship. These helpers add
# a purely positional/lexical check (no LLM, no parsing) that generalizes to
# any number of organs mentioned in one sentence:
#
#  1. _list_groups() groups organ mentions that are joined ONLY by
#     enumeration glue (",", "and", "or", "/", "&", "as well as", articles)
#     with nothing else between them -- "heart, colon, and kidney" becomes
#     one group regardless of how many members it has. Two organs in the
#     same group are siblings in a list, not claimed to interact with each
#     other.
#  2. For organs in different groups, a pair is rejected if some OTHER organ
#     (not a fellow list-member of either side) is mentioned before the
#     later of the two -- the classic "X regulates ... A and B" pattern --
#     while organs mentioned only AFTER both (e.g. a trailing contrastive
#     "...unlike the kidney") are not held against the pair.
ALL_ORGAN_PATTERNS = {o: _compile_patterns(aliases) for o, aliases in ORGAN_ALIASES.items()}

_LIST_GLUE_RE = re.compile(r"^\s*(,|;|&|/)?\s*(and|or|as well as|&|/)?\s*(the|a|an)?\s*$")


def _find_organ_mentions(sent: str, all_organ_patterns: dict) -> list[tuple[int, int, str]]:
    """All (start, end, organ_name) mentions in sent, sorted by position."""
    hits = [(m.start(), m.end(), organ)
            for organ, pats in all_organ_patterns.items()
            for p in pats for m in p.finditer(sent)]
    hits.sort(key=lambda t: t[0])
    return hits


def _list_groups(mentions: list[tuple[int, int, str]], sent: str) -> dict[str, int]:
    """Assign each organ an enumeration-list group id. Consecutive mentions
    (in sentence order) separated only by list glue ("," "and" "/" ...) get
    the same id; anything else (a verb, a clause) starts a new group."""
    groups: dict[str, int] = {}
    group_id = 0
    prev_end = None
    for start, end, organ in mentions:
        if prev_end is not None and not _LIST_GLUE_RE.match(sent[prev_end:start]):
            group_id += 1
        groups.setdefault(organ, group_id)
        prev_end = end
    return groups


def _passes_multi_organ_check(sent: str, organ1: str, organ2: str,
                               all_organ_patterns: dict) -> bool:
    """False if the sentence's phrasing indicates organ1/organ2 are not
    directly linked to each other (see module comment above)."""
    mentions = _find_organ_mentions(sent, all_organ_patterns)
    if len(mentions) <= 2:
        return True  # nothing else mentioned -- no ambiguity to resolve

    groups = _list_groups(mentions, sent)
    if groups.get(organ1) == groups.get(organ2):
        return False  # co-listed siblings ("...in the heart and colon")

    first_pos: dict[str, int] = {}
    for start, _end, organ in mentions:
        first_pos.setdefault(organ, start)
    later = max(first_pos[organ1], first_pos[organ2])
    own_groups = {groups[organ1], groups[organ2]}

    for start, _end, organ in mentions:
        if organ in (organ1, organ2) or start >= later:
            continue
        if groups.get(organ) in own_groups:
            continue  # fellow list-member of organ1/organ2, not a third party
        return False  # a genuine third organ precedes the later of the two
    return True


def _filter_same_sentence_crosstalk(papers: list[dict], organ1_patterns: list,
                                    organ2_patterns: list, crosstalk_patterns: list,
                                    organ1_name: str | None = None,
                                    organ2_name: str | None = None,
                                    all_organ_patterns: dict | None = None) -> list[dict]:
    """
    Keep only papers where both organs' aliases AND at least one crosstalk
    keyword (network, axis, interplay, ...) appear together in the same
    sentence — not just somewhere across the same title/abstract. The
    PubMed query already requires both organs and a crosstalk keyword to
    match somewhere in the document (cheap pre-filter, fewer papers
    fetched); this is the strict local check that actually enforces
    three-way co-location.

    When organ1_name/organ2_name/all_organ_patterns are given, a qualifying
    sentence is further required to pass _passes_multi_organ_check() --
    guarding against sentences that mention a third organ and are really
    describing that third organ's relationship to organ1/organ2, not a
    direct organ1<->organ2 link (see module comment above).
    """
    kept = []
    for paper in papers:
        text = (paper.get("title", "") + ". " + paper.get("abstract", "")).lower()
        sentences = re.split(r"(?<=[.!?;])\s+", text)
        for sent in sentences:
            if not (_any_match(organ1_patterns, sent)
                    and _any_match(organ2_patterns, sent)
                    and _any_match(crosstalk_patterns, sent)):
                continue
            if (all_organ_patterns is not None and organ1_name and organ2_name
                    and not _passes_multi_organ_check(sent, organ1_name, organ2_name,
                                                       all_organ_patterns)):
                continue
            kept.append(paper)
            break
    return kept


def _filter_by_condition(papers: list[dict], condition_patterns: list) -> list[dict]:
    """
    Keep only papers whose title+abstract contains at least one condition
    keyword (e.g. "healthy", "obesity") anywhere in the text.

    Applied locally to papers that have ALREADY been fetched and passed
    through _filter_same_sentence_crosstalk -- no new PubMed query. This
    mirrors the Edge_cosine CONDITION_FILTER's matching semantics (whole
    title/abstract, not restricted to one sentence) but as a post-hoc
    filter on already-cached papers instead of a query-time clause.
    """
    kept = []
    for paper in papers:
        text = (paper.get("title", "") + ". " + paper.get("abstract", "")).lower()
        if _any_match(condition_patterns, text):
            kept.append(paper)
    return kept


def build_condition_results(condition: str, base_results: dict, pairs: list[tuple[str, str]]) -> dict:
    """
    Build a condition-specific results dict from the already-searched "all"
    results: for each pair in `pairs`, filter its already-extracted papers
    down to those also matching one of that condition's keywords, then
    recompute key players from the smaller, condition-filtered paper set.

    Pairs with no cached entry in base_results, or that have zero papers
    after the condition filter, still get an entry (n_papers_found=0) so
    downstream steps (key player counts, LLM classification/description)
    see a consistent pair set matching `pairs`.
    """
    condition_patterns = _compile_patterns(CONDITION_KEYWORDS[condition])

    base_by_pair: dict[tuple[str, str], dict] = {}
    for entry in base_results.values():
        o1, o2 = entry.get("organ1", ""), entry.get("organ2", "")
        if o1 and o2:
            base_by_pair[(min(o1, o2), max(o1, o2))] = entry

    results: dict = {}
    for (o1, o2) in pairs:
        base_entry = base_by_pair.get((o1, o2))
        if base_entry is None:
            print(f"  [!] no cached 'all' entry for {o1} <-> {o2} -- run the base search first")
            continue

        all_papers = base_entry.get("papers", [])
        cond_papers = _filter_by_condition(all_papers, condition_patterns)
        key_players = extract_key_players(cond_papers)

        print(f"  {o1} <-> {o2}: {len(all_papers)} papers -> {len(cond_papers)} after "
              f"'{condition}' filter")

        results[f"{o1}|{o2}"] = {
            "organ1":                 o1,
            "organ2":                 o2,
            "pubmed_query":           base_entry.get("pubmed_query", ""),
            "strategy_used":          f"{base_entry.get('strategy_used', '')}+condition:{condition}",
            "n_papers_found":         len(cond_papers),
            "n_papers_found_raw":     len(all_papers),
            "papers":                 cond_papers,
            "key_players":            key_players,
            "search_date":            datetime.now().isoformat(),
        }

    return results


# ── PubMed search ─────────────────────────────────────────────────────────────

def search_pubmed(query: str, max_results: int, years_back: int) -> list[str]:
    """Return up to max_results PMIDs for query."""
    min_date = (datetime.now() - timedelta(days=365 * years_back)).strftime("%Y/%m/%d")
    params = {
        "db":      "pubmed",
        "term":    query,
        "retmax":  max_results,
        "retmode": "json",
        "mindate": min_date,
        "datetype": "pdat",
    }
    # POST avoids the URL-length ceiling GET requests hit once the pair's
    # combined organ+metabolic+crosstalk clauses push the query past a few
    # thousand characters — NCBI returns HTTP 414 on GET in that case, which
    # silently degrades that pair's results to zero rather than erroring loudly.
    resp = _ncbi_get("esearch.fcgi", params, method="POST")
    if resp is None:
        return []
    data = resp.json()
    return data.get("esearchresult", {}).get("idlist", [])


# ── Main search loop ──────────────────────────────────────────────────────────

def run_search(organ_pairs, output_path, resume=True, force_empty=False):
    output_path = Path(output_path)
    results: dict = {}

    if resume and output_path.exists():
        with open(output_path, encoding="utf-8") as f:
            results = json.load(f)
        n_cached = len(results)
        n_empty  = sum(1 for v in results.values() if v.get("n_papers_found", 0) == 0)
        print(f"[i] Resuming: {n_cached} cached ({n_empty} empty).")
        if force_empty:
            empty_keys = [k for k, v in results.items() if v.get("n_papers_found", 0) == 0]
            for k in empty_keys:
                del results[k]
            print(f"[i] Cleared {len(empty_keys)} empty entries for re-search.")

    crosstalk_patterns = _compile_patterns(CROSSTALK_KEYWORDS)

    total = len(organ_pairs)
    for idx, (o1, o2) in enumerate(organ_pairs):
        edge_key = f"{o1}|{o2}"
        sym_key  = f"{o2}|{o1}"

        if edge_key in results or sym_key in results:
            n = (results.get(edge_key) or results.get(sym_key, {})).get("n_papers_found", 0)
            print(f"  [{idx+1}/{total}] Cached ({n} papers): {o1} <-> {o2}")
            continue

        query = build_query(o1, o2)
        print(f"\n  [{idx+1}/{total}] {o1} <-> {o2}")
        print(f"    Query: {query[:120]}...")
        time.sleep(DELAY)

        pmids = search_pubmed(query, MAX_PAPERS, YEARS_BACK)
        print(f"    Found {len(pmids)} PMIDs — fetching abstracts...", flush=True)
        time.sleep(DELAY)

        papers_raw = fetch_abstracts(pmids, delay=DELAY)
        organ1_patterns = _compile_patterns(ORGAN_ALIASES.get(o1, [o1]))
        organ2_patterns = _compile_patterns(ORGAN_ALIASES.get(o2, [o2]))
        papers = _filter_same_sentence_crosstalk(papers_raw, organ1_patterns,
                                                 organ2_patterns, crosstalk_patterns,
                                                 organ1_name=o1, organ2_name=o2,
                                                 all_organ_patterns=ALL_ORGAN_PATTERNS)
        key_players = extract_key_players(papers)

        n = len(papers)
        print(
            f"    => {len(papers_raw)} papers retrieved, {n} kept "
            f"(both organs + crosstalk term in same sentence)"
            f" | hormones: {len(key_players['hormones'])}"
            f" | metabolites: {len(key_players['metabolites'])}"
            f" | proteins: {len(key_players['proteins'])}"
        )

        results[edge_key] = {
            "organ1":          o1,
            "organ2":          o2,
            "pubmed_query":    query,
            "strategy_used":   "metabolic-single-query",
            "n_papers_found":  n,
            "n_papers_found_raw": len(papers_raw),
            "papers":          papers,
            "key_players":     key_players,
            "search_date":     datetime.now().isoformat(),
        }

        with open(output_path, "w", encoding="utf-8") as f:
            json.dump(results, f, indent=2, ensure_ascii=False)

    total_papers = sum(v.get("n_papers_found", 0) for v in results.values())
    empty = sum(1 for v in results.values() if v.get("n_papers_found", 0) == 0)
    print(f"\n[ok] Done. {len(results)} edges | {total_papers} total papers | {empty} empty.")
    print(f"     Saved: {output_path}")
    return results


# ── Visualization ─────────────────────────────────────────────────────────────

def _literature_stats_section(organs: list[str], results_by_pair: dict, allowed_pairs: set,
                               variant_id: str = "all") -> str:
    """
    "Literature Statistics" content for one variant (all/healthy/obese):
    per-organ and per-connection paper counts, plus a canvas-drawn circular
    Sankey diagram. Ribbon width = same-sentence co-occurring paper count
    for that organ pair (further condition-filtered for the healthy/obese
    variants -- see condition_variants in build_viz()).

    Unlike the per-organ cosine pipelines, this pipeline searches per PAIR,
    not per organ, so there's no single fetch a "papers per organ" count
    could be read directly from. It's derived here as the union of PMIDs
    across every pair touching that organ (a paper supporting two of an
    organ's connections is only counted once).

    variant_id must be unique across the three calls in one page (used to
    namespace the canvas element id -- otherwise three "lit-sankey" ids in
    one DOM would collide and only one canvas would ever get its diagram).
    """
    from Visualisation.networkBuilderUtils import ORGAN_COLORS, DEFAULT_NODE_COLOR

    organ_pmids: dict[str, set] = {o: set() for o in organs}
    for pair in allowed_pairs:
        o1, o2 = pair
        data = results_by_pair.get(pair, {})
        if data.get("n_papers_found", 0) <= 0:
            continue
        pmids = {p.get("pmid") for p in data.get("papers", []) if p.get("pmid")}
        organ_pmids.setdefault(o1, set()).update(pmids)
        organ_pmids.setdefault(o2, set()).update(pmids)
    organ_papers = {o: len(organ_pmids.get(o, set())) for o in organs}

    links = []
    for pair in allowed_pairs:
        o1, o2 = pair
        data = results_by_pair.get(pair, {})
        n = data.get("n_papers_found", 0)
        if n <= 0:
            continue
        # Consistent left/right assignment: the organ with the larger total
        # paper count is the "source" (ribbon origin).
        if organ_papers.get(o1, 0) < organ_papers.get(o2, 0):
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

    canvas_id = f"lit-sankey-{variant_id}"
    return f"""
    <div class="info-h2">Organ Cross-Talk Sankey</div>
    <p class="info-p" style="color:#94a3b8">
      Ribbon thickness = number of same-sentence co-occurring papers for that
      organ pair.
    </p>
    <div style="overflow-x:auto">
      <canvas id="{canvas_id}" style="display:block"></canvas>
    </div>
    <script>
    (function() {{
      const D = {sankey_json};
      const canvas = document.getElementById('{canvas_id}');
      if (!D.links.length) {{ canvas.style.display = 'none'; return; }}
      const ctx = canvas.getContext('2d');

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

      let angle = -90;
      const organRange = {{}};
      organs.forEach(o => {{
        const span = (totals[o] / grand) * usableDeg;
        organRange[o] = {{ start: angle, end: angle + span }};
        angle += span + GAP_DEG;
      }});

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


def _literature_stats_tab_content(organs: list[str], results_by_pair: dict, allowed_pairs: set,
                                   condition_data: dict[str, tuple | None] | None) -> str:
    """
    Wraps three _literature_stats_section() variants (All / Healthy /
    Obese) with a small selector so the Literature Statistics tab can show
    condition-filtered stats side by side with the full search -- same
    layout and Sankey diagram for each, just computed from different data
    (see condition_data in _build_info_tabs()).
    """
    condition_data = condition_data or {}

    variant_sections = {
        "all": _literature_stats_section(organs, results_by_pair, allowed_pairs, variant_id="all"),
    }
    for cond in ("healthy", "obese"):
        data = condition_data.get(cond)
        if data is None:
            variant_sections[cond] = (
                f'<p class="info-p" style="color:#64748b">'
                f'The {cond} condition pipeline has not been run yet '
                f'(<code>uv run python run_metabolic_lit_search.py --condition {cond}</code>) '
                f'-- no condition-filtered literature statistics are available.</p>'
            )
        else:
            cond_organs, cond_results_by_pair, cond_pairs = data
            variant_sections[cond] = _literature_stats_section(
                cond_organs, cond_results_by_pair, cond_pairs, variant_id=cond
            )

    labels = {"all": "All", "healthy": "Healthy", "obese": "Obese"}
    selector_buttons = "".join(
        f'<button class="lit-stats-selector-btn{" active" if key == "all" else ""}" '
        f'data-variant="{key}" onclick="showLitStatsVariant(\'{key}\')">{labels[key]}</button>'
        for key in ("all", "healthy", "obese")
    )
    variant_panels = "".join(
        f'<div class="lit-stats-variant" id="lit-stats-variant-{key}" '
        f'style="display:{"block" if key == "all" else "none"}">{content}</div>'
        for key, content in variant_sections.items()
    )

    return f"""
    <div class="info-h2">Literature Statistics</div>
    <p class="info-p" style="color:#94a3b8">
      Choose which paper set these statistics are computed from: the full
      condition-agnostic search, or one of the two condition-filtered views
      (see "Condition Filter" under Literature Search).
    </p>
    <div style="display:flex;gap:6px;margin-bottom:14px">
      {selector_buttons}
    </div>
    {variant_panels}
    <style>
      .lit-stats-selector-btn {{
        background:#1e293b;border:1px solid #334155;color:#94a3b8;
        border-radius:6px;padding:5px 12px;font-size:0.76rem;cursor:pointer;
      }}
      .lit-stats-selector-btn:hover {{ background:#334155;color:#e2e8f0; }}
      .lit-stats-selector-btn.active {{ background:#6366f1;border-color:#6366f1;color:#fff; }}
    </style>
    <script>
    function showLitStatsVariant(variant) {{
      document.querySelectorAll('.lit-stats-variant').forEach(el => {{
        el.style.display = (el.id === 'lit-stats-variant-' + variant) ? 'block' : 'none';
      }});
      document.querySelectorAll('.lit-stats-selector-btn').forEach(btn => {{
        btn.classList.toggle('active', btn.dataset.variant === variant);
      }});
    }}
    </script>
"""


def _build_info_tabs(organs: list[str], results_by_pair: dict, allowed_pairs: set,
                      condition_data: dict[str, tuple | None] | None = None) -> list[dict]:
    """
    Build info panel content dynamically from the current configuration.

    condition_data (optional): {"healthy": (organs, results_by_pair, pairs) | None,
    "obese": (...) | None} -- when present, the Literature Statistics tab
    gets an All/Healthy/Obese selector using this data for the latter two
    (see build_condition_results() / run_metabolic_lit_search.py
    --condition). A condition mapped to None (pipeline not yet run for it)
    is shown as unavailable rather than omitted, so the selector's presence
    doesn't silently depend on which steps happen to have been run.
    """
    condition_data = condition_data or {}
    kw_block = " OR ".join(
        f'"{kw}"' if (" " in kw or "-" in kw) else kw
        for kw in METABOLIC_KEYWORDS
    )
    crosstalk_kw_block = " OR ".join(
        f'"{kw}"' if (" " in kw or "-" in kw) else kw
        for kw in CROSSTALK_KEYWORDS
    )
    healthy_kw_block = ", ".join(f'"{kw}"' for kw in CONDITION_KEYWORDS["healthy"])
    obese_kw_block = ", ".join(f'"{kw}"' for kw in CONDITION_KEYWORDS["obese"])
    from Literature_Search.llm_connection_type import VOTE_RUNS as vote_runs
    vote_majority = vote_runs // 2 + 1
    ct_rows = "".join(
        f'<tr><td style="padding:3px 8px 3px 0;color:#94a3b8;white-space:nowrap">'
        f'<strong style="color:#e2e8f0">{v["label"]}</strong></td>'
        f'<td style="padding:3px 0;color:#94a3b8">{v["description"]}</td></tr>'
        for v in CONNECTION_TYPES.values()
    )
    return [
        {
            "id": "search",
            "label": "Literature Search",
            "content": f"""
        <div class="info-h2">Overview</div>
        <p class="info-p">
          This is the <strong>glucose metabolism-based reference network</strong>.
          Organ–organ connections (both the Healthy and Obese predefined edge sets
          selectable via the topbar toggle) are based on partial correlation networks
          from
          <a href="https://doi.org/10.1016/j.medj.2025.100881" target="_blank"
             style="color:#818cf8;text-decoration:underline"
             title="Geist, B. K. et al. The metabolic organ connectome: A novel approach to measure allostatic load during health-to-disease transition. Med 6, 100881 (2025)">Geist et al. (2025)</a>
          -- the Healthy cohort is n=241, BMI &lt;24.4 kg/m².
          For each included pair the pipeline issues
          <strong>one PubMed query</strong> requiring both organs, at least one
          metabolic keyword, and at least one crosstalk keyword to all
          co-occur somewhere in the document — a coarse pre-filter. After
          fetching, papers are filtered further: kept only if both organs'
          names/aliases <strong>and</strong> a crosstalk keyword appear
          together in the <strong>same sentence</strong>, not merely
          somewhere in the same title/abstract.
        </p>
        <div class="info-h2">Search Parameters</div>
        <div class="info-stat-grid">
          <div class="info-stat"><div class="info-stat-val">{YEARS_BACK} yrs</div><div class="info-stat-lbl">Look-back window</div></div>
          <div class="info-stat"><div class="info-stat-val">{MAX_PAPERS:,}</div><div class="info-stat-lbl">Max papers per pair</div></div>
          <div class="info-stat"><div class="info-stat-val">single</div><div class="info-stat-lbl">Query strategy</div></div>
          <div class="info-stat"><div class="info-stat-val" id="stat-val-edges">—</div><div class="info-stat-lbl">Edges with ≥1 paper</div></div>
        </div>
        <div class="info-h2">Query Construction</div>
        <p class="info-p">Each organ contributes its <strong>MeSH term</strong> and all
        <strong>text aliases</strong> searched in Title/Abstract. Both organ clauses must
        co-occur with at least one metabolic keyword and at least one crosstalk keyword:</p>
        <div class="info-code">(MeSH_A OR aliases_A[Title/Abstract])
AND (MeSH_B OR aliases_B[Title/Abstract])
AND ({kw_block})
AND ({crosstalk_kw_block})</div>
        <p class="info-p">
          Organ and crosstalk terms are matched as <strong>whole words</strong>,
          not substrings — "renal" no longer matches inside "adrenal", so a
          paper about the adrenal gland can't be mistaken for kidney evidence.
          A match is also discarded if it's <strong>negated</strong> —
          "non-renal", "nonrenal", "non renal", and "not renal" don't count
          as a positive mention of that organ or crosstalk term.
        </p>
        <div class="info-h2">Same-Sentence Filter</div>
        <p class="info-p">
          The query above is a cheap pre-filter — it only guarantees both
          organs and a crosstalk keyword appear <em>somewhere</em> in the
          document. After fetching, a paper is kept only if there's a single
          sentence in its title/abstract containing <strong>organ A's
          name/alias, organ B's name/alias, and a crosstalk keyword all
          together</strong>. A paper discussing organ A in one paragraph and
          organ B in a completely separate one no longer counts as evidence
          for that pair, even though it passed the PubMed query.
        </p>
        <div class="info-h2">Condition Filter (Healthy / Obese)</div>
        <p class="info-p">
          The Healthy and Obese views shown via the topbar toggle are <strong>not</strong>
          separate PubMed searches — there is only ever one shared, condition-agnostic
          search per pair (above). Instead, a condition-specific view is built by taking
          that pair's already same-sentence-filtered papers and keeping only the ones
          whose title/abstract <strong>also</strong> contains at least one of that
          condition's keywords, matched as a whole word anywhere in the text (not
          restricted to the same sentence as the organ/crosstalk match):
        </p>
        <div class="info-code">Healthy keywords: {healthy_kw_block}
Obese keywords:   {obese_kw_block}</div>
        <p class="info-p">
          This condition-filtered paper set is what feeds the key player counts,
          LLM connection-type classification, and LLM edge descriptions shown when
          Healthy or Obese is selected — each condition can therefore show a different
          paper count, different key players, and even a different connection type or
          description for the same organ pair than the All view. A pair with no papers
          left after this filter is dropped from that condition's view (but stays
          available under All).
        </p>
        <div class="info-h2">Key Player Extraction</div>
        <p class="info-p">All abstracts are scanned against three curated vocabulary lists:</p>
        <ul style="margin:0 0 10px 16px;padding:0;color:#cbd5e1">
          <li style="margin-bottom:4px"><strong style="color:#fdba74">Hormones</strong> — insulin, glucagon, cortisol, leptin, GLP-1, IGF-1, FGF21, irisin, …</li>
          <li style="margin-bottom:4px"><strong style="color:#5eead4">Metabolites</strong> — glucose, lactate, fatty acids, bile acids, amino acids, ATP, …</li>
          <li style="margin-bottom:4px"><strong style="color:#a5b4fc">Proteins</strong> — GLUT1–5, AMPK, mTOR, PPARs, PGC-1α, CPT1, LPL, …</li>
        </ul>
        <p class="info-p">Key players are <strong>ranked by mention count</strong> and shown for each edge.</p>
""",
        },
        {
            "id": "literature_stats",
            "label": "Literature Statistics",
            "content": _literature_stats_tab_content(organs, results_by_pair, allowed_pairs, condition_data),
        },
        {
            "id": "llm",
            "label": "Connection Types",
            "content": f"""
        <div class="info-h2">LLM Classification</div>
        <p class="info-p">
          For each organ pair an LLM (<code>{LLM_MODEL}</code>, via Ollama) reads up to
          <strong>{LLM_MAX_PAPERS} paper abstracts</strong> and assigns
          <strong>1 to 3 connection types</strong> from the categories below, ranked
          from most to least feasible — a pair only gets a second or third type if it is
          independently well-supported by the papers.
        </p>
        <p class="info-p">
          To reduce single-sample noise, the classification is run
          <strong>{vote_runs} times independently</strong> per pair; only types that appear
          in at least {vote_majority} of the {vote_runs} runs are kept (majority vote). If the
          runs disagree on everything, the single most-voted type is kept, so every pair
          always ends up with at least one type.
        </p>
        <div class="info-h2">Connection Type Categories</div>
        <table style="width:100%;border-collapse:collapse;font-size:12px">
          {ct_rows}
        </table>
        <p class="info-p" style="margin-top:10px">
          All assigned types are shown as coloured badges on each edge in the network.
        </p>
""",
        },
        {
            "id": "descriptions",
            "label": "Descriptions",
            "content": f"""
        <div class="info-h2">Edge Descriptions</div>
        <p class="info-p">
          Each organ–organ connection receives a short, cited scientific summary generated
          by <code>{LLM_MODEL}</code> in a two-step process:
        </p>
        <ol style="margin:0 0 12px 16px;padding:0;color:#cbd5e1;font-size:13px;line-height:1.7">
          <li style="margin-bottom:6px">
            <strong style="color:#e2e8f0">Sample</strong>: up to
            <strong>25 papers</strong> are randomly sampled from the full PubMed result
            set for that pair (if fewer papers exist, all are used).
          </li>
          <li style="margin-bottom:6px">
            <strong style="color:#e2e8f0">Select</strong>: the LLM reads the 25 paper
            titles and abstracts and identifies the <strong>5 most relevant</strong> ones
            for the specific metabolic or hormonal interaction between the two organs.
          </li>
          <li>
            <strong style="color:#e2e8f0">Summarise</strong>: the LLM writes a
            3–4 sentence cited summary based solely on those 5 papers.
            Citations in the text correspond to paper numbers <strong>[1]–[5]</strong>
            shown in the edge panel.
          </li>
        </ol>
        <p class="info-p">
          This ensures the description always cites a curated, relevant subset rather than
          generic highly-cited papers that may not address the specific connection.
        </p>
        <div class="info-h2">Organ Descriptions</div>
        <p class="info-p">
          Each organ node receives a 5-sentence metabolic overview generated by
          <code>{LLM_MODEL}</code>:
        </p>
        <ol style="margin:0 0 12px 16px;padding:0;color:#cbd5e1;font-size:13px;line-height:1.7">
          <li style="margin-bottom:6px">
            <strong style="color:#e2e8f0">Search</strong>: PubMed is queried for the
            organ combined with metabolic keywords, retrieving up to <strong>50 papers</strong>
            from the last <strong>5 years</strong>.
          </li>
          <li style="margin-bottom:6px">
            <strong style="color:#e2e8f0">Summarise</strong>: the LLM writes exactly
            <strong>5 sentences</strong> covering the organ's primary metabolic substrates,
            energy production pathways, and key hormonal regulation signals.
          </li>
        </ol>
        <p class="info-p">
          Organ descriptions appear in the node panel when you click on an organ in the network.
        </p>
""",
        },
    ]


def _index_results_by_pair(results: dict) -> dict[tuple[str, str], dict]:
    by_pair: dict[tuple[str, str], dict] = {}
    for v in results.values():
        o1, o2 = v.get("organ1", ""), v.get("organ2", "")
        if o1 and o2:
            by_pair[(min(o1, o2), max(o1, o2))] = v
    return by_pair


def _compute_edge_display_data(o1: str, o2: str, edge_data: dict,
                                llm_types: dict, llm_descs: dict,
                                get_connection_types, CONNECTION_TYPES) -> dict:
    """
    Build the display-ready dict for one edge (key players, connection
    type, LLM description, papers) from a given set of literature/LLM
    sources. Reused for both the "all" baseline and each condition variant
    (healthy/obese) -- same shape, different source data.
    """
    papers  = edge_data.get("papers", [])
    _kp_raw = edge_data.get("key_players", {})

    def _remerge(raw_counts, synonyms):
        merged = _merge_synonyms(raw_counts, synonyms)
        ranked = [t for t, _ in sorted(merged.items(), key=lambda x: -x[1])]
        return ranked, merged

    _h_list, _h_counts = _remerge(_kp_raw.get("hormones_counts", {}), HORMONE_SYNONYMS)
    _m_list, _m_counts = _remerge(_kp_raw.get("metabolites_counts", {}), METABOLITE_SYNONYMS)
    _p_list, _p_counts = _remerge(_kp_raw.get("proteins_counts", {}), PROTEIN_SYNONYMS)

    type_labels = get_connection_types(llm_types, o1, o2, CONNECTION_TYPES)
    conn_type        = type_labels[0] if type_labels else ""
    conn_type_others = type_labels[1:]

    llm_entry = (
        llm_descs.get(f"{o1}|{o2}")
        or llm_descs.get(f"{o2}|{o1}")
        or {}
    )
    # When an LLM summary exists, its citation numbers [1]-[5] index into
    # its OWN selected_papers list, not the full evidence list -- "papers"
    # below must match whichever text is actually shown, or citation links
    # resolve to the wrong paper / nothing.
    llm_description = llm_entry.get("description", "")
    citation_papers = llm_entry.get("papers") if llm_description else None

    return {
        "pubmed_query":              edge_data.get("pubmed_query", ""),
        "n_papers_found":            edge_data.get("n_papers_found", 0),
        "papers":                    citation_papers if citation_papers else papers,
        "connection_type":        conn_type,
        "connection_type_others": conn_type_others,
        "search_date":               edge_data.get("search_date", ""),
        "ai_description":            llm_description,
        "key_players_hormones":           _h_list,
        "key_players_metabolites":        _m_list,
        "key_players_proteins":           _p_list,
        "key_players_counts_hormones":    _h_counts,
        "key_players_counts_metabolites": _m_counts,
        "key_players_counts_proteins":    _p_counts,
    }


def build_viz(results: dict, output_path: Path):
    """Build the interactive HTML network visualization."""
    try:
        import networkx as nx
        from Visualisation.networkBuilderUtils import export_network_to_cytoscape_dashboard
        from Literature_Search.llm_connection_type import (
            load_connection_type_classifications, get_connection_types,
        )
        from Literature_Search.llm_descriptions import load_llm_descriptions
        from Literature_Search.organ_descriptions import load_organ_descriptions

        allowed_pairs  = set(load_edge_filter(EDGE_FILTER_CSV))
        allowed_organs = {o for pair in allowed_pairs for o in pair}

        # Healthy cohort subset used only to flag which nodes/edges the
        # dashboard's "Healthy" display toggle should show -- it no longer
        # restricts what was searched (see EDGE_FILTER_CSV above).
        healthy_pairs  = set(load_edge_filter(HEALTHY_FILTER_CSV)) if HEALTHY_FILTER_CSV.exists() else set()
        healthy_organs = {o for pair in healthy_pairs for o in pair}
        obese_pairs    = set(load_edge_filter(OBESE_FILTER_CSV)) if OBESE_FILTER_CSV.exists() else set()
        obese_organs   = {o for pair in obese_pairs for o in pair}

        llm_types   = load_connection_type_classifications(OUTPUT_LLM_TYPES)
        llm_descs   = load_llm_descriptions(LLM_DESCRIPTIONS)
        organ_descs = load_organ_descriptions(ORGAN_DESCRIPTIONS)

        # Index results by canonical pair key for fast lookup
        results_by_pair = _index_results_by_pair(results)

        # Condition-specific sources (healthy/obese), built by
        # run_metabolic_lit_search.py --condition healthy|obese. Missing
        # files just mean no condition_variants get attached for that
        # condition (the toggle still shows/hides pairs, just without
        # condition-specific key players/connection type/description).
        condition_sources: dict[str, tuple | None] = {}
        for cond in ("healthy", "obese"):
            paths = condition_output_paths(cond)
            if not paths["results"].exists():
                condition_sources[cond] = None
                continue
            with open(paths["results"], encoding="utf-8") as f:
                cond_results = json.load(f)
            cond_results_by_pair = _index_results_by_pair(cond_results)
            cond_llm_types = (load_connection_type_classifications(paths["llm_types"])
                               if paths["llm_types"].exists() else {})
            cond_llm_descs = (load_llm_descriptions(paths["llm_descs"])
                               if paths["llm_descs"].exists() else {})
            condition_sources[cond] = (cond_results_by_pair, cond_llm_types, cond_llm_descs)

        G = nx.Graph()
        for organ in sorted(allowed_organs):
            organ_entry = organ_descs.get(organ, {})
            G.add_node(organ,
                       llm_description=organ_entry.get("description", ""),
                       llm_papers=organ_entry.get("papers", []),
                       is_healthy=organ in healthy_organs,
                       is_obese=organ in obese_organs)

        for (o1, o2) in allowed_pairs:
            G.add_edge(o1, o2)
            is_h = (o1, o2) in healthy_pairs
            is_o = (o1, o2) in obese_pairs
            G.edges[o1, o2]['is_healthy'] = is_h
            G.edges[o1, o2]['is_obese']   = is_o

            edge_data = results_by_pair.get((o1, o2), {})
            G.edges[o1, o2]['merged_data'] = _compute_edge_display_data(
                o1, o2, edge_data, llm_types, llm_descs,
                get_connection_types, CONNECTION_TYPES,
            )
            G.edges[o1, o2]['color'] = "#64748b"

            # Condition variants: only computed for pairs that actually
            # belong to that condition's cohort, and only if that
            # condition's pipeline has been run at least once.
            variants = {}
            for cond, is_member in (("healthy", is_h), ("obese", is_o)):
                if not is_member or condition_sources[cond] is None:
                    continue
                cr_by_pair, c_types, c_descs = condition_sources[cond]
                cond_edge_data = cr_by_pair.get((o1, o2), {})
                variants[cond] = _compute_edge_display_data(
                    o1, o2, cond_edge_data, c_types, c_descs,
                    get_connection_types, CONNECTION_TYPES,
                )
            if variants:
                G.edges[o1, o2]['condition_variants'] = variants

        info_condition_data = {}
        for cond, cond_organs, cond_pairs in (("healthy", healthy_organs, healthy_pairs),
                                               ("obese", obese_organs, obese_pairs)):
            if condition_sources[cond] is None:
                info_condition_data[cond] = None
            else:
                cr_by_pair, _c_types, _c_descs = condition_sources[cond]
                info_condition_data[cond] = (sorted(cond_organs), cr_by_pair, cond_pairs)

        export_network_to_cytoscape_dashboard(
            graph=G,
            filename=str(output_path),
            include_legend=False,
            title=VIZ_TITLE,
            info_panel_tabs=_build_info_tabs(sorted(allowed_organs), results_by_pair, allowed_pairs,
                                             condition_data=info_condition_data),
            healthy_toggle_label="Healthy",
            obese_toggle_label="Obese",
            comparison_upload_label="Upload comparison network",
        )
        print(f"[ok] Visualization saved: {output_path}")
    except Exception as e:
        import traceback
        print(f"[!] Visualization failed: {e}")
        traceback.print_exc()


# ── Condition pipeline (--condition healthy|obese) ─────────────────────────────

def run_condition_pipeline(condition: str, args) -> None:
    """
    Build a condition-filtered results file from the already-cached "all"
    search results (no new PubMed calls), then run LLM connection-type
    classification and LLM connection descriptions against it.

    Entirely separate from the normal flow below: does not touch
    OUTPUT_JSON/OUTPUT_LLM_TYPES/LLM_DESCRIPTIONS and does not rebuild the
    dashboard -- the "all" case stays exactly as it is today.
    """
    if not OUTPUT_JSON.exists():
        print(f"[!] No base results file found at {OUTPUT_JSON}.")
        print(f"    Run without --condition first to build the base 'all' search.")
        sys.exit(1)
    with open(OUTPUT_JSON, encoding="utf-8") as f:
        base_results = json.load(f)
    print(f"[i] Loaded {len(base_results)} cached 'all' edges.")

    cohort_csv = CONDITION_COHORT_CSV[condition]
    if not cohort_csv.exists():
        print(f"[!] Cohort CSV not found: {cohort_csv}")
        sys.exit(1)
    cond_pairs = load_edge_filter(cohort_csv)
    print(f"[i] Condition: {condition} | Cohort: {cohort_csv.name} → {len(cond_pairs)} pairs")
    print(f"[i] Condition keywords: {', '.join(CONDITION_KEYWORDS[condition])}\n")

    out_paths = condition_output_paths(condition)

    print("[i] Building condition-filtered results (filtering already-fetched papers, "
          "no new PubMed queries)...")
    cond_results = build_condition_results(condition, base_results, cond_pairs)
    with open(out_paths["results"], "w", encoding="utf-8") as f:
        json.dump(cond_results, f, indent=2, ensure_ascii=False)
    total_papers = sum(v.get("n_papers_found", 0) for v in cond_results.values())
    print(f"\n[ok] {len(cond_results)} edges | {total_papers} papers after condition filter.")
    print(f"     Saved: {out_paths['results']}")

    if not args.skip_llm_type:
        print(f"\n[i] Running LLM connection-type classification ({condition})...")
        from Literature_Search.llm_connection_type import generate_connection_type_classifications
        generate_connection_type_classifications(
            organ_pairs         = cond_pairs,
            literature_results  = cond_results,
            connection_types    = CONNECTION_TYPES,
            output_path         = out_paths["llm_types"],
            model               = LLM_MODEL,
            max_papers          = LLM_MAX_PAPERS,
            resume              = not args.reset_llm_type,
            reset               = args.reset_llm_type,
        )
    else:
        print("[i] Skipping LLM type classification (using cache).")

    print(f"\n[i] Running LLM connection descriptions ({condition})...")
    from Literature_Search.llm_descriptions import generate_llm_descriptions
    generate_llm_descriptions(
        organ_pairs         = cond_pairs,
        literature_results  = cond_results,
        output_path         = out_paths["llm_descs"],
        model               = LLM_MODEL,
        resume              = not args.reset,
        reset               = args.reset,
    )

    print(f"\n[ok] Condition pipeline done for '{condition}'.")
    print("     Note: the dashboard is not yet wired to these condition-specific files -- "
          "it still shows the same 'all' data regardless of the Healthy/Obese toggle.")


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--reset",          action="store_true",
                        help="Delete search cache and start from scratch.")
    parser.add_argument("--force-empty",    action="store_true",
                        help="Re-search edges that previously returned 0 papers.")
    parser.add_argument("--skip-llm-type",  action="store_true",
                        help="Skip LLM connection-type classification; use cached results.")
    parser.add_argument("--reset-llm-type", action="store_true",
                        help="Delete LLM type cache and reclassify all pairs.")
    parser.add_argument("--viz-only",       action="store_true",
                        help="Skip search and LLM; just rebuild the visualization.")
    parser.add_argument("--condition", choices=["healthy", "obese"],
                        help="Build a condition-filtered results file from the already-cached "
                             "'all' search (no new PubMed queries) plus its own LLM "
                             "connection-type classification and descriptions. Separate "
                             "pipeline -- does not touch the 'all' case or the dashboard.")
    args = parser.parse_args()

    if args.condition:
        run_condition_pipeline(args.condition, args)
        return

    if args.reset and OUTPUT_JSON.exists():
        OUTPUT_JSON.unlink()
        print(f"[i] Search cache deleted: {OUTPUT_JSON}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    if not EDGE_FILTER_CSV.exists():
        print(f"[!] Edge filter not found: {EDGE_FILTER_CSV}")
        sys.exit(1)
    pairs = load_edge_filter(EDGE_FILTER_CSV)
    print(f"[i] Edge filter: {EDGE_FILTER_CSV.name} → {len(pairs)} pairs selected")
    print(f"[i] Max {MAX_PAPERS} papers each | {YEARS_BACK} years back\n")

    # ── Step 1: PubMed search ─────────────────────────────────────────────────
    if not args.viz_only:
        results = run_search(pairs, OUTPUT_JSON,
                             resume=not args.reset,
                             force_empty=args.force_empty)
    else:
        if not OUTPUT_JSON.exists():
            print(f"[!] No results file found at {OUTPUT_JSON}. Run without --viz-only first.")
            sys.exit(1)
        with open(OUTPUT_JSON, encoding="utf-8") as f:
            results = json.load(f)
        print(f"[i] Loaded {len(results)} cached edges.")

    # ── Step 2: LLM connection-type classification ────────────────────────────
    if not args.viz_only and not args.skip_llm_type:
        print("\n[i] Running LLM connection-type classification...")
        from Literature_Search.llm_connection_type import generate_connection_type_classifications
        generate_connection_type_classifications(
            organ_pairs       = pairs,
            literature_results = results,
            connection_types  = CONNECTION_TYPES,
            output_path       = OUTPUT_LLM_TYPES,
            model             = LLM_MODEL,
            max_papers        = LLM_MAX_PAPERS,
            resume            = not args.reset_llm_type,
            reset             = args.reset_llm_type,
        )
    elif args.skip_llm_type or args.viz_only:
        print("[i] Skipping LLM type classification (using cache).")

    # ── Step 3: Visualization ─────────────────────────────────────────────────
    print("\n[i] Building visualization...")
    build_viz(results, OUTPUT_HTML)


if __name__ == "__main__":
    main()
