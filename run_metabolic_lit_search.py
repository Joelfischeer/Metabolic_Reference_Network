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

def _lit_sankey_html(variant_id: str, title: str, sankey_json: str) -> str:
    """Canvas-drawn circular Sankey card for one variant."""
    canvas_id = f"lit-sankey-{variant_id}"
    return f"""
    <div class="ls-card">
    <div class="ls-card-title">{title}</div>
    <div style="text-align:center">
      <canvas id="{canvas_id}" style="display:block"></canvas>
    </div>
    <script>
    (function() {{
      const D = {sankey_json};
      const canvas = document.getElementById('{canvas_id}');
      if (!D.links.length) {{ canvas.style.display = 'none'; return; }}
      const ctx = canvas.getContext('2d');

      const organs = D.organs, links = D.links;
      const SIZE = 640, CX = SIZE / 2, CY = SIZE / 2, R = 200, ARC_W = 14, LABEL_GAP = 10;
      const GAP_DEG = organs.length > 1 ? 1.4 : 0;

      const SCALE = 2;
      canvas.width = SIZE * SCALE; canvas.height = SIZE * SCALE;
      canvas.style.width = '100%'; canvas.style.maxWidth = '460px'; canvas.style.height = 'auto'; canvas.style.margin = '0 auto';
      ctx.scale(SCALE, SCALE);
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
        ctx.font = '14px system-ui';
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
    </div>
"""


def _lit_lists_html(label: str, organ_rows: str, link_rows: str) -> str:
    """Papers per organ / per connection lists for one variant."""
    return f"""
    <div style="display:flex;gap:18px;flex-wrap:wrap">
      <div style="flex:1;min-width:190px">
        <div class="ls-list-title">Papers per organ</div>
        {organ_rows}
      </div>
      <div style="flex:1;min-width:210px">
        <div class="ls-list-title">Papers per connection</div>
        <div class="ls-scroll">
        {link_rows if link_rows else '<p style="color:#64748b;font-size:0.8rem">No cross-mention evidence yet.</p>'}
        </div>
      </div>
    </div>
"""


def _literature_stats_section(organs: list[str], results_by_pair: dict, allowed_pairs: set) -> dict:
    """
    Literature-statistics data for one variant (all/healthy/obese): the
    circular-Sankey data plus per-organ and per-connection paper counts.
    Ribbon width = same-sentence co-occurring paper count for that organ pair
    (further condition-filtered for the healthy/obese variants -- see
    condition_variants in build_viz()).

    This pipeline searches per PAIR, not per organ, so "papers per organ" is
    derived as the union of PMIDs across every pair touching that organ (a
    paper supporting two of an organ's connections is only counted once).

    Returns {"sankey_json": str, "organ_rows": html, "link_rows": html}.
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
        n = results_by_pair.get(pair, {}).get("n_papers_found", 0)
        if n <= 0:
            continue
        # The organ with the larger total paper count is the ribbon "source".
        if organ_papers.get(o1, 0) < organ_papers.get(o2, 0):
            o1, o2 = o2, o1
        links.append({"source": o1, "target": o2, "value": n})

    row = ('<div style="display:flex;justify-content:space-between;padding:3px 0;'
           'border-bottom:1px solid #1e293b;font-size:0.78rem">'
           '<span style="color:#e2e8f0">{a}</span>'
           '<span style="color:#64748b">{n:,} papers</span></div>')
    organ_rows = "".join(row.format(a=o, n=organ_papers[o])
                         for o in sorted(organs, key=lambda o: -organ_papers[o]))
    link_rows = "".join(row.format(a=f'{l["source"]} ↔ {l["target"]}', n=l["value"])
                        for l in sorted(links, key=lambda l: -l["value"]))

    sankey_organs = [o for o in sorted(organs, key=lambda o: -organ_papers[o])
                     if any(l["source"] == o or l["target"] == o for l in links)]
    sankey_json = json.dumps({
        "organs":      sankey_organs,
        "organColors": {o: ORGAN_COLORS.get(o, DEFAULT_NODE_COLOR) for o in sankey_organs},
        "links":       links,
    }, ensure_ascii=False)
    return {"sankey_json": sankey_json, "organ_rows": organ_rows, "link_rows": link_rows}


def _lit_stats_tab_wrap(variants: dict) -> str:
    """
    Two-column Literature Statistics layout. Left: condition selector with
    papers per organ / per connection for the selected condition. Right: the
    Sankey diagrams for All, Healthy and Obese, stacked and always visible.

    variants: key -> {"sankey_json", "organ_rows", "link_rows"}, or
    {"message": html} for a condition whose pipeline has not been run.
    """
    labels = {"all": "All", "healthy": "Healthy", "obese": "Obese"}
    keys = [k for k in ("all", "healthy", "obese") if k in variants]

    buttons = "".join(
        f'<button class="lit-stats-selector-btn{" active" if k == "all" else ""}" '
        f'data-variant="{k}" onclick="showLitStatsVariant(\'{k}\')">{labels[k]}</button>'
        for k in keys)

    lists, cards = "", ""
    for k in keys:
        v = variants[k]
        display = "block" if k == "all" else "none"
        if "message" in v:
            lists += (f'<div class="lit-stats-variant" id="lit-stats-variant-{k}" '
                      f'style="display:{display}">{v["message"]}</div>')
            cards += (f'<div class="ls-card"><div class="ls-card-title">{labels[k]}</div>'
                      f'{v["message"]}</div>')
        else:
            lists += (f'<div class="lit-stats-variant" id="lit-stats-variant-{k}" '
                      f'style="display:{display}">'
                      f'{_lit_lists_html(labels[k], v["organ_rows"], v["link_rows"])}</div>')
            cards += _lit_sankey_html(k, labels[k], v["sankey_json"])

    return f"""
    <style>
      .ls-split {{ display:grid; grid-template-columns:minmax(280px,5fr) minmax(300px,6fr);
                   gap:28px; align-items:start; }}
      @media (max-width:900px) {{ .ls-split {{ grid-template-columns:1fr; }} }}
      .ls-left {{ position:sticky; top:0; }}
      .ls-right {{ display:flex; flex-direction:column; gap:12px; min-width:0; }}
      .ls-card {{ background:#0f172a; border:1px solid #334155; border-radius:8px;
                  padding:8px 10px 10px; }}
      .ls-card-title {{ font-size:13px; font-weight:700; color:#e2e8f0; margin-bottom:4px; }}
      .ls-list-title {{ font-size:0.74rem; color:#94a3b8; margin-bottom:6px;
                        text-transform:uppercase; letter-spacing:0.03em; }}
      .ls-scroll {{ max-height:calc(100vh - 300px); min-height:160px; overflow-y:auto; padding-right:6px; }}
      .lit-stats-selector-btn {{
        background:#1e293b;border:1px solid #334155;color:#94a3b8;
        border-radius:6px;padding:5px 12px;font-size:0.76rem;cursor:pointer;
      }}
      .lit-stats-selector-btn:hover {{ background:#334155;color:#e2e8f0; }}
      .lit-stats-selector-btn.active {{ background:#6366f1;border-color:#6366f1;color:#fff; }}
    </style>
    <div class="ls-split">
      <div class="ls-left">
        <div class="info-h2">Literature Statistics</div>
        <p class="info-p" style="color:#94a3b8">
          Choose which paper set the counts are computed from: the full
          condition-agnostic search, or one of the two condition-filtered views
          (see "Condition Filter" under Literature Search). The Sankey diagrams
          on the right always show all three; ribbon thickness = number of
          same-sentence co-occurring papers for that organ pair.
        </p>
        <div style="display:flex;gap:6px;margin-bottom:14px">{buttons}</div>
        {lists}
      </div>
      <div class="ls-right">{cards}</div>
    </div>
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


def _literature_stats_tab_content(organs: list[str], results_by_pair: dict, allowed_pairs: set,
                                   condition_data: dict[str, tuple | None] | None) -> str:
    """
    Literature Statistics tab: All / Healthy / Obese variants computed from
    different paper sets (see condition_data in _build_info_tabs()), laid out
    by _lit_stats_tab_wrap().
    """
    condition_data = condition_data or {}

    variants = {"all": _literature_stats_section(organs, results_by_pair, allowed_pairs)}
    for cond in ("healthy", "obese"):
        data = condition_data.get(cond)
        if data is None:
            variants[cond] = {"message": (
                f'<p class="info-p" style="color:#64748b">'
                f'The {cond} condition pipeline has not been run yet '
                f'(<code>uv run python run_metabolic_lit_search.py --condition {cond}</code>) '
                f'-- no condition-filtered literature statistics are available.</p>')}
        else:
            cond_organs, cond_results_by_pair, cond_pairs = data
            variants[cond] = _literature_stats_section(cond_organs, cond_results_by_pair, cond_pairs)
    return _lit_stats_tab_wrap(variants)


def _venn_svg(n_only_a: int, n_only_b: int, n_both: int,
              label_a: str, color_a: str,
              label_b: str, color_b: str,
              n_outer: int | None = None,
              outer_label: str = "All",
              subtitle_a: str = "", subtitle_b: str = "") -> str:
    """
    Euler diagram: outer rounded-rect (All), two overlapping circles (A and B).
    All counts are integers; positions are fixed (non-proportional) for clarity.
    """
    W, H = 420, 210
    cx_a, cx_b, cy, r = 155, 265, 108, 82

    # Region label x positions
    lx_a   = cx_a - r // 2 - 12   # left of left circle
    lx_ab  = (cx_a + cx_b) // 2   # centre of overlap
    lx_b   = cx_b + r // 2 + 12   # right of right circle
    ly_num = cy + 5
    ly_lbl = cy + 22

    outer_text = ""
    if n_outer is not None:
        outer_text = (
            f'<text x="22" y="{H - 14}" fill="#475569" font-size="10" font-family="inherit">'
            f'{outer_label} only: {n_outer:,}</text>'
        )

    total_a = n_only_a + n_both
    total_b = n_only_b + n_both

    sub_a = subtitle_a or f"total {total_a:,}"
    sub_b = subtitle_b or f"total {total_b:,}"

    return (
        f'<svg viewBox="0 0 {W} {H}" xmlns="http://www.w3.org/2000/svg"'
        f' style="width:100%;max-width:{W}px;display:block;margin:10px 0">'
        # Background (All)
        f'<rect x="4" y="4" width="{W-8}" height="{H-8}" rx="12"'
        f' fill="#0f172a" stroke="#334155" stroke-width="1.5"/>'
        f'<text x="16" y="20" fill="#475569" font-size="10" font-weight="600"'
        f' font-family="inherit" text-anchor="start">{outer_label}</text>'
        # Left circle (A)
        f'<circle cx="{cx_a}" cy="{cy}" r="{r}"'
        f' fill="{color_a}" fill-opacity="0.18" stroke="{color_a}" stroke-width="1.5"/>'
        # Right circle (B)
        f'<circle cx="{cx_b}" cy="{cy}" r="{r}"'
        f' fill="{color_b}" fill-opacity="0.18" stroke="{color_b}" stroke-width="1.5"/>'
        # Label A (top of left circle)
        f'<text x="{cx_a}" y="{cy - r - 10}" text-anchor="middle" fill="{color_a}"'
        f' font-size="11" font-weight="700" font-family="inherit">{label_a}</text>'
        f'<text x="{cx_a}" y="{cy - r + 2}" text-anchor="middle" fill="{color_a}"'
        f' font-size="9" font-family="inherit" fill-opacity="0.75">{sub_a}</text>'
        # Label B (top of right circle)
        f'<text x="{cx_b}" y="{cy - r - 10}" text-anchor="middle" fill="{color_b}"'
        f' font-size="11" font-weight="700" font-family="inherit">{label_b}</text>'
        f'<text x="{cx_b}" y="{cy - r + 2}" text-anchor="middle" fill="{color_b}"'
        f' font-size="9" font-family="inherit" fill-opacity="0.75">{sub_b}</text>'
        # Count: A only
        f'<text x="{lx_a}" y="{ly_num}" text-anchor="middle" fill="#e2e8f0"'
        f' font-size="18" font-weight="700" font-family="inherit">{n_only_a:,}</text>'
        f'<text x="{lx_a}" y="{ly_lbl}" text-anchor="middle" fill="#94a3b8"'
        f' font-size="9" font-family="inherit">only</text>'
        # Count: overlap
        f'<text x="{lx_ab}" y="{ly_num}" text-anchor="middle" fill="#e2e8f0"'
        f' font-size="18" font-weight="700" font-family="inherit">{n_both:,}</text>'
        f'<text x="{lx_ab}" y="{ly_lbl}" text-anchor="middle" fill="#94a3b8"'
        f' font-size="9" font-family="inherit">shared</text>'
        # Count: B only
        f'<text x="{lx_b}" y="{ly_num}" text-anchor="middle" fill="#e2e8f0"'
        f' font-size="18" font-weight="700" font-family="inherit">{n_only_b:,}</text>'
        f'<text x="{lx_b}" y="{ly_lbl}" text-anchor="middle" fill="#94a3b8"'
        f' font-size="9" font-family="inherit">only</text>'
        f'{outer_text}'
        f'</svg>'
    )


def _network_analysis_tab_html(results_by_pair: dict,
                                condition_data: dict,
                                all_organs: list[str]) -> str:
    """Build the Network Analysis tab content."""

    def _net_stats(rbp: dict) -> tuple[set, set]:
        """Return (edges_set, nodes_set) for pairs with ≥1 paper."""
        edges, nodes = set(), set()
        for v in rbp.values():
            if v.get("n_papers_found", 0) > 0:
                o1, o2 = v.get("organ1", ""), v.get("organ2", "")
                if o1 and o2:
                    edges.add((min(o1, o2), max(o1, o2)))
                    nodes.update([o1, o2])
        return edges, nodes

    edges_all, nodes_all = _net_stats(results_by_pair)
    n_all_organs = len(all_organs)
    max_edges = n_all_organs * (n_all_organs - 1) // 2

    cd_h = condition_data.get("healthy")
    cd_o = condition_data.get("obese")
    edges_h, nodes_h = _net_stats(cd_h[1]) if cd_h else (set(), set())
    edges_o, nodes_o = _net_stats(cd_o[1]) if cd_o else (set(), set())

    def _density(e: int, n: int) -> str:
        return f"{100*e/(n*(n-1)/2):.1f}%" if n > 1 else "—"

    def _isolated(all_organ_list: list, node_set: set) -> list:
        return sorted(o for o in all_organ_list if o not in node_set)

    iso_all  = _isolated(all_organs, nodes_all)
    iso_h    = _isolated(all_organs, nodes_h)
    iso_o    = _isolated(all_organs, nodes_o)

    def _iso_cell(iso: list) -> str:
        return ", ".join(iso) if iso else '<em style="color:#475569">none</em>'

    # Pair-set overlaps
    both_h_o   = edges_h & edges_o
    only_h     = edges_h - edges_o
    only_o     = edges_o - edges_h
    neither    = edges_all - edges_h - edges_o

    # Venn diagram for edge overlap
    edges_venn_svg = _venn_svg(
        n_only_a=len(only_h),  n_only_b=len(only_o),  n_both=len(both_h_o),
        label_a="Healthy",     color_a="#34d399",
        label_b="Obese",       color_b="#fb923c",
        n_outer=len(neither),  outer_label=f"All ({len(edges_all)} edges)",
        subtitle_a=f"{len(edges_h)} edges", subtitle_b=f"{len(edges_o)} edges",
    )

    # Per-organ degree table
    def _degrees(edges: set, organ_list: list) -> dict[str, int]:
        deg = {o: 0 for o in organ_list}
        for o1, o2 in edges:
            if o1 in deg: deg[o1] += 1
            if o2 in deg: deg[o2] += 1
        return deg

    deg_all = _degrees(edges_all, all_organs)
    deg_h   = _degrees(edges_h,   all_organs)
    deg_o   = _degrees(edges_o,   all_organs)

    # Sort organs by degree in "All" view descending
    sorted_organs = sorted(all_organs, key=lambda o: -deg_all[o])

    degree_rows = ""
    for organ in sorted_organs:
        da, dh, do_ = deg_all[organ], deg_h[organ], deg_o[organ]
        bar_a = f'<span style="display:inline-block;width:{int(da/max(deg_all.values(),default=1)*60)}px;height:8px;background:#818cf8;border-radius:2px;vertical-align:middle"></span>'
        bar_h = f'<span style="display:inline-block;width:{int(dh/max(deg_all.values(),default=1)*60)}px;height:8px;background:#34d399;border-radius:2px;vertical-align:middle"></span>'
        bar_o = f'<span style="display:inline-block;width:{int(do_/max(deg_all.values(),default=1)*60)}px;height:8px;background:#fb923c;border-radius:2px;vertical-align:middle"></span>'
        degree_rows += (
            f'<tr>'
            f'<td style="padding:3px 10px 3px 0;color:#cbd5e1;white-space:nowrap">{organ}</td>'
            f'<td style="padding:3px 10px 3px 0;text-align:right;color:#e2e8f0;font-weight:600">{da}</td>'
            f'<td style="padding:3px 10px 3px 0">{bar_a}</td>'
            f'<td style="padding:3px 10px 3px 0;text-align:right;color:#e2e8f0;font-weight:600">{dh}</td>'
            f'<td style="padding:3px 10px 3px 0">{bar_h}</td>'
            f'<td style="padding:3px 10px 3px 0;text-align:right;color:#e2e8f0;font-weight:600">{do_}</td>'
            f'<td style="padding:3px 0">{bar_o}</td>'
            f'</tr>'
        )

    th = 'style="padding:3px 10px 6px 0;color:#64748b;font-weight:600;text-align:right;border-bottom:1px solid #334155"'
    th_l = 'style="padding:3px 10px 6px 0;color:#64748b;font-weight:600;text-align:left;border-bottom:1px solid #334155"'

    return f"""
        <div class="info-h2">Network Summary</div>
        <p class="info-p">Each view is treated as an unweighted undirected network where a
        pair counts as a connection if it has ≥1 literature paper after that view's filter.</p>
        <table style="width:100%;border-collapse:collapse;font-size:12px;margin-bottom:14px">
          <thead>
            <tr>
              <th {th_l}>View</th>
              <th {th}>Organs</th>
              <th {th}>Edges</th>
              <th {th}>Max edges</th>
              <th {th}>Density</th>
              <th style="padding:3px 0 6px;color:#64748b;font-weight:600;text-align:left;border-bottom:1px solid #334155">Isolated organs</th>
            </tr>
          </thead>
          <tbody>
            <tr>
              <td style="padding:4px 10px 4px 0;color:#a5b4fc;font-weight:600">All</td>
              <td style="padding:4px 10px 4px 0;text-align:right;color:#e2e8f0">{len(nodes_all)}</td>
              <td style="padding:4px 10px 4px 0;text-align:right;color:#e2e8f0;font-weight:700">{len(edges_all)}</td>
              <td style="padding:4px 10px 4px 0;text-align:right;color:#64748b">{max_edges}</td>
              <td style="padding:4px 10px 4px 0;text-align:right;color:#e2e8f0">{_density(len(edges_all), len(nodes_all))}</td>
              <td style="padding:4px 0;color:#94a3b8;font-size:11px">{_iso_cell(iso_all)}</td>
            </tr>
            <tr>
              <td style="padding:4px 10px 4px 0;color:#34d399;font-weight:600">Healthy</td>
              <td style="padding:4px 10px 4px 0;text-align:right;color:#e2e8f0">{len(nodes_h)}</td>
              <td style="padding:4px 10px 4px 0;text-align:right;color:#e2e8f0;font-weight:700">{len(edges_h)}</td>
              <td style="padding:4px 10px 4px 0;text-align:right;color:#64748b">{max_edges}</td>
              <td style="padding:4px 10px 4px 0;text-align:right;color:#e2e8f0">{_density(len(edges_h), len(nodes_h))}</td>
              <td style="padding:4px 0;color:#94a3b8;font-size:11px">{_iso_cell(iso_h)}</td>
            </tr>
            <tr>
              <td style="padding:4px 10px 4px 0;color:#fb923c;font-weight:600">Obese</td>
              <td style="padding:4px 10px 4px 0;text-align:right;color:#e2e8f0">{len(nodes_o)}</td>
              <td style="padding:4px 10px 4px 0;text-align:right;color:#e2e8f0;font-weight:700">{len(edges_o)}</td>
              <td style="padding:4px 10px 4px 0;text-align:right;color:#64748b">{max_edges}</td>
              <td style="padding:4px 10px 4px 0;text-align:right;color:#e2e8f0">{_density(len(edges_o), len(nodes_o))}</td>
              <td style="padding:4px 0;color:#94a3b8;font-size:11px">{_iso_cell(iso_o)}</td>
            </tr>
          </tbody>
        </table>
        <div class="info-h2">Edge Overlap Between Conditions</div>
        <p class="info-p">How the {len(edges_h)} healthy edges and {len(edges_o)} obese edges
        relate to each other and to the full network of {len(edges_all)} edges:</p>
        <table style="width:100%;border-collapse:collapse;font-size:12px;margin-bottom:14px">
          <thead>
            <tr>
              <th {th_l}>Edge subset</th>
              <th {th}>Count</th>
              <th style="padding:3px 0 6px;color:#64748b;font-weight:600;text-align:left;border-bottom:1px solid #334155">Description</th>
            </tr>
          </thead>
          <tbody>
            <tr>
              <td style="padding:4px 10px 4px 0;color:#34d399;white-space:nowrap">Healthy only</td>
              <td style="padding:4px 10px 4px 0;text-align:right;color:#e2e8f0;font-weight:700">{len(only_h)}</td>
              <td style="padding:4px 0;color:#94a3b8;font-size:11px">In Healthy but absent from Obese</td>
            </tr>
            <tr>
              <td style="padding:4px 10px 4px 0;color:#fb923c;white-space:nowrap">Obese only</td>
              <td style="padding:4px 10px 4px 0;text-align:right;color:#e2e8f0;font-weight:700">{len(only_o)}</td>
              <td style="padding:4px 0;color:#94a3b8;font-size:11px">In Obese but absent from Healthy</td>
            </tr>
            <tr>
              <td style="padding:4px 10px 4px 0;color:#a5b4fc;white-space:nowrap">Healthy ∩ Obese</td>
              <td style="padding:4px 10px 4px 0;text-align:right;color:#e2e8f0;font-weight:700">{len(both_h_o)}</td>
              <td style="padding:4px 0;color:#94a3b8;font-size:11px">Present in both Healthy and Obese</td>
            </tr>
            <tr>
              <td style="padding:4px 10px 4px 0;color:#64748b;white-space:nowrap">Neither condition</td>
              <td style="padding:4px 10px 4px 0;text-align:right;color:#e2e8f0;font-weight:700">{len(neither)}</td>
              <td style="padding:4px 0;color:#94a3b8;font-size:11px">In full network but pass no condition filter</td>
            </tr>
          </tbody>
        </table>
        <div class="info-h2">Per-Organ Degree</div>
        <p class="info-p">Number of connections each organ has in each view (sorted by All-view degree):</p>
        <table style="width:100%;border-collapse:collapse;font-size:12px;margin-bottom:6px">
          <thead>
            <tr>
              <th {th_l}>Organ</th>
              <th {th} colspan="2" style="color:#a5b4fc">All</th>
              <th {th} colspan="2" style="color:#34d399">Healthy</th>
              <th {th} colspan="2" style="color:#fb923c">Obese</th>
            </tr>
          </thead>
          <tbody>
            {degree_rows}
          </tbody>
        </table>
"""


def _build_info_tabs(organs: list[str], results_by_pair: dict, allowed_pairs: set,
                      condition_data: dict[str, tuple | None] | None = None,
                      oo_data: dict | None = None) -> list[dict]:
    """
    Build info panel content dynamically from the current configuration.

    condition_data (optional): {"healthy": (organs, results_by_pair, pairs) | None,
    "obese": (...) | None} -- when present, the Literature Statistics tab
    gets an All/Healthy/Obese selector using this data for the latter two
    (see build_condition_results() / run_metabolic_lit_search.py
    --condition). A condition mapped to None (pipeline not yet run for it)
    is shown as unavailable rather than omitted, so the selector's presence
    doesn't silently depend on which steps happen to have been run.

    oo_data (optional): output of _compute_oo_data -- adds an "Edge
    Statistics" tab with OO coefficient distribution and elbow threshold.
    """
    condition_data = condition_data or {}

    # ── Compute paper-funnel stats (with duplicate accounting) ──────────────
    _n_raw      = sum(v.get("n_papers_found_raw", v.get("n_papers_found", 0))
                      for v in results_by_pair.values())
    _n_filtered = sum(v.get("n_papers_found", 0) for v in results_by_pair.values())
    _pairs_filtered = sum(1 for v in results_by_pair.values() if v.get("n_papers_found", 0) > 0)

    def _unique_pmids(rbp: dict) -> set:
        pmids = set()
        for v in rbp.values():
            for p in v.get("papers", []):
                pid = p.get("pmid")
                if pid:
                    pmids.add(pid)
        return pmids

    _pmids_filtered     = _unique_pmids(results_by_pair)
    _n_unique_filtered  = len(_pmids_filtered)
    _n_dup_filtered     = _n_filtered - _n_unique_filtered

    def _cond_stats(cond_name: str) -> tuple[int, int, int, int]:
        """Returns (total_papers, unique_papers, pairs_with_papers, dup_count)."""
        cd = condition_data.get(cond_name)
        if cd is None:
            return 0, 0, 0, 0
        _, cond_rbp, _ = cd
        total  = sum(v.get("n_papers_found", 0) for v in cond_rbp.values())
        unique = len(_unique_pmids(cond_rbp))
        pairs  = sum(1 for v in cond_rbp.values() if v.get("n_papers_found", 0) > 0)
        return total, unique, pairs, total - unique

    _n_healthy, _n_uniq_healthy, _pairs_healthy, _n_dup_healthy = _cond_stats("healthy")
    _n_obese,   _n_uniq_obese,   _pairs_obese,   _n_dup_obese   = _cond_stats("obese")

    def _pct(num: int, denom: int) -> str:
        return f"{100*num/denom:.0f}%" if denom else "—"
    def _dup_label(dup: int, total: int) -> str:
        return f"{dup:,} ({100*dup/total:.0f}%)" if total else "—"

    _pct_filt    = _pct(_n_filtered, _n_raw)
    _pct_healthy = _pct(_n_healthy, _n_filtered)
    _pct_obese   = _pct(_n_obese,   _n_filtered)

    # Unique pmids shared between healthy and obese
    _pmids_healthy = (_unique_pmids(condition_data["healthy"][1])
                      if condition_data.get("healthy") else set())
    _pmids_obese   = (_unique_pmids(condition_data["obese"][1])
                      if condition_data.get("obese") else set())
    _n_shared_ho     = len(_pmids_healthy & _pmids_obese)
    _n_only_h_papers = _n_uniq_healthy - _n_shared_ho
    _n_only_o_papers = _n_uniq_obese   - _n_shared_ho
    _n_outer_papers  = _n_unique_filtered - _n_uniq_healthy - _n_uniq_obese + _n_shared_ho

    _paper_venn_svg = _venn_svg(
        n_only_a=_n_only_h_papers, n_only_b=_n_only_o_papers, n_both=_n_shared_ho,
        label_a="Healthy",         color_a="#34d399",
        label_b="Obese",           color_b="#fb923c",
        n_outer=_n_outer_papers,
        outer_label=f"All ({_n_unique_filtered:,} unique papers)",
        subtitle_a=f"{_n_uniq_healthy:,} unique", subtitle_b=f"{_n_uniq_obese:,} unique",
    )

    _funnel_rows = f"""
        <tr>
          <td class="frow-step">PubMed fetch</td>
          <td class="frow-val">{_n_raw:,}</td>
          <td class="frow-uniq" style="color:#64748b;font-style:italic">n/a</td>
          <td class="frow-pairs">{len(results_by_pair)}</td>
          <td class="frow-pct">—</td>
        </tr>
        <tr>
          <td class="frow-step">Same-sentence filter</td>
          <td class="frow-val">{_n_filtered:,}</td>
          <td class="frow-uniq">{_n_unique_filtered:,}</td>
          <td class="frow-pairs">{_pairs_filtered}</td>
          <td class="frow-pct">{_pct_filt}</td>
        </tr>
        <tr>
          <td class="frow-step">Healthy view</td>
          <td class="frow-val">{_n_healthy:,}</td>
          <td class="frow-uniq">{_n_uniq_healthy:,}</td>
          <td class="frow-pairs">{_pairs_healthy}</td>
          <td class="frow-pct">{_pct_healthy}</td>
        </tr>
        <tr>
          <td class="frow-step">Obese view</td>
          <td class="frow-val">{_n_obese:,}</td>
          <td class="frow-uniq">{_n_uniq_obese:,}</td>
          <td class="frow-pairs">{_pairs_obese}</td>
          <td class="frow-pct">{_pct_obese}</td>
        </tr>"""

    _dup_section = f"""
        <div class="info-h2">Duplicate Papers</div>
        <p class="info-p">
          Each paper is fetched independently for every organ pair whose PubMed
          query it matches, so the same PMID can appear in multiple pairs
          (e.g., a paper on liver–gut–brain crosstalk may pass the same-sentence
          filter for both Liver–Colon and Liver–Brain). The counts above reflect
          this redundancy; the table below shows how many entries are duplicates
          (same PMID, different pair):
        </p>
        <table style="width:100%;border-collapse:collapse;font-size:12px;margin-bottom:10px">
          <thead>
            <tr>
              <th style="padding:3px 8px 6px 0;color:#64748b;font-weight:600;text-align:left;border-bottom:1px solid #334155">View</th>
              <th style="padding:3px 12px 6px 0;color:#64748b;font-weight:600;text-align:right;border-bottom:1px solid #334155">Total</th>
              <th style="padding:3px 12px 6px 0;color:#64748b;font-weight:600;text-align:right;border-bottom:1px solid #334155">Unique</th>
              <th style="padding:3px 0 6px;color:#64748b;font-weight:600;text-align:right;border-bottom:1px solid #334155">Duplicates</th>
            </tr>
          </thead>
          <tbody>
            <tr>
              <td class="frow-step">All (same-sentence)</td>
              <td class="frow-val">{_n_filtered:,}</td>
              <td class="frow-uniq">{_n_unique_filtered:,}</td>
              <td class="frow-pct">{_dup_label(_n_dup_filtered, _n_filtered)}</td>
            </tr>
            <tr>
              <td class="frow-step">Healthy view</td>
              <td class="frow-val">{_n_healthy:,}</td>
              <td class="frow-uniq">{_n_uniq_healthy:,}</td>
              <td class="frow-pct">{_dup_label(_n_dup_healthy, _n_healthy)}</td>
            </tr>
            <tr>
              <td class="frow-step">Obese view</td>
              <td class="frow-val">{_n_obese:,}</td>
              <td class="frow-uniq">{_n_uniq_obese:,}</td>
              <td class="frow-pct">{_dup_label(_n_dup_obese, _n_obese)}</td>
            </tr>
          </tbody>
        </table>
        <p class="info-p">
          Additionally, <strong>{_n_shared_ho:,} unique papers</strong> appear in
          both the Healthy and Obese views — papers whose abstracts contain keywords
          from both conditions (e.g., a study comparing lean vs. obese subjects).
        </p>
"""
    # ────────────────────────────────────────────────────────────────────────

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
          This is the <strong>14-organ metabolic reference network</strong> covering
          <strong>{len(allowed_pairs)} organ–organ pairs</strong> across
          <strong>{len(organs)} metabolically active organs</strong>: adrenal
          glands, bone marrow, brain, gut, heart, kidney, liver, lung, muscle,
          pancreas, small intestine, spleen, thyroid, and white adipose tissue (WAT).
        </p>
        <p class="info-p">
          Three views are available via the topbar:
          <strong>No Condition</strong> — all pairs with ≥1 literature paper
          (the full network);
          <strong>Healthy</strong> — pairs whose papers mention healthy/normal
          physiology keywords;
          <strong>Obese</strong> — pairs whose papers mention obesity/adiposity
          keywords.
          An additional comparison overlay (⚖ buttons) lets you highlight how the
          literature-defined condition edges differ from the predefined cohort
          connections from
          <a href="https://doi.org/10.1016/j.medj.2025.100881" target="_blank"
             style="color:#818cf8;text-decoration:underline"
             title="Geist, B. K. et al. The metabolic organ connectome: A novel approach to measure allostatic load during health-to-disease transition. Med 6, 100881 (2025)">Geist et al. (2025)</a>
          (healthy cohort n=241, BMI &lt;24.4 kg/m²).
        </p>
        <p class="info-p">
          For each pair the pipeline issues
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
        <div class="info-h2">Paper Funnel</div>
        <p class="info-p" style="margin-bottom:8px">Papers surviving each step of the filtering pipeline across all {len(results_by_pair)} organ–organ pairs. "Total" counts the same paper once per pair it was assigned to; "Unique" counts distinct PMIDs:</p>
        <table style="width:100%;border-collapse:collapse;font-size:12px;margin-bottom:10px">
          <style>
            .frow-step{{padding:4px 8px 4px 0;color:#cbd5e1;white-space:nowrap}}
            .frow-val{{padding:4px 8px 4px 0;color:#f8fafc;font-weight:700;text-align:right;white-space:nowrap}}
            .frow-uniq{{padding:4px 8px 4px 0;color:#93c5fd;text-align:right;white-space:nowrap}}
            .frow-pairs{{padding:4px 8px 4px 0;color:#94a3b8;text-align:right;white-space:nowrap}}
            .frow-pct{{padding:4px 0;color:#64748b;text-align:right;white-space:nowrap}}
          </style>
          <thead>
            <tr>
              <th style="padding:3px 8px 6px 0;color:#64748b;font-weight:600;text-align:left;border-bottom:1px solid #334155">Step</th>
              <th style="padding:3px 8px 6px 0;color:#64748b;font-weight:600;text-align:right;border-bottom:1px solid #334155">Total</th>
              <th style="padding:3px 8px 6px 0;color:#64748b;font-weight:600;text-align:right;border-bottom:1px solid #334155">Unique</th>
              <th style="padding:3px 8px 6px 0;color:#64748b;font-weight:600;text-align:right;border-bottom:1px solid #334155">Pairs</th>
              <th style="padding:3px 0 6px;color:#64748b;font-weight:600;text-align:right;border-bottom:1px solid #334155">Retained</th>
            </tr>
          </thead>
          <tbody>{_funnel_rows}
          </tbody>
        </table>
        {_dup_section}
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
        {
            "id": "network_analysis",
            "label": "Network Analysis",
            "content": _network_analysis_tab_html(results_by_pair, condition_data, organs),
        },
    ] + ([
        {
            "id": "oo_threshold",
            "label": "Edge Statistics",
            "content": _oo_threshold_tab_html(oo_data),
        },
    ] if oo_data else [])


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


def _compute_oo_data(results_by_pair: dict, allowed_pairs: set,
                      cond_pats: dict | None = None) -> dict:
    """
    Compute Otsuka-Ochiai coefficients and elbow threshold for the
    reference-metabolic network, for the "all" view and (optionally) each
    condition-filtered view.

    Per-organ paper count = union of PMIDs across every pair touching that
    organ (same derivation used for the Sankey diagram).

    Returns a dict keyed by condition name ("all", "healthy", "obese"), each
    containing:
        pairs:     sorted list of (organ1, organ2, oo_value)
        threshold: float  (kneedle elbow on sorted OO values, or 0 if N<3)
        n_edges:   int    (pairs with oo_value >= threshold and > 0)
    """
    from threshold_utils import kneedle_elbow

    def _oo_for_papers(pairs_papers: dict[tuple, list]) -> dict:
        """pairs_papers: {(o1,o2): [paper, ...]}  →  result dict."""
        organ_pmids: dict[str, set] = {}
        for (o1, o2), papers in pairs_papers.items():
            pmids = {p.get("pmid") for p in papers if p.get("pmid")}
            organ_pmids.setdefault(o1, set()).update(pmids)
            organ_pmids.setdefault(o2, set()).update(pmids)
        organ_n = {o: len(s) for o, s in organ_pmids.items()}

        oo_list = []
        for (o1, o2), papers in pairs_papers.items():
            n_co = len([p for p in papers if p.get("pmid")])
            n1 = organ_n.get(o1, 0)
            n2 = organ_n.get(o2, 0)
            denom = (n1 * n2) ** 0.5 if n1 > 0 and n2 > 0 else 0
            oo = (n_co / denom) if denom > 0 else 0.0
            oo_list.append((o1, o2, oo))
        oo_list.sort(key=lambda t: -t[2])

        vals = [t[2] for t in oo_list]
        nonzero = sorted(v for v in vals if v > 0)
        thresh = kneedle_elbow(nonzero) if len(nonzero) >= 3 else 0.0
        n_edges = sum(1 for v in vals if v >= thresh and v > 0)
        return {"pairs": oo_list, "threshold": thresh, "n_edges": n_edges}

    # "all" condition: use already-filtered (same-sentence) papers
    all_pp = {pair: results_by_pair.get(pair, {}).get("papers", [])
              for pair in allowed_pairs}
    result = {"all": _oo_for_papers(all_pp)}

    if cond_pats:
        for cond, pats in cond_pats.items():
            cond_pp = {pair: _filter_by_condition(papers, pats)
                       for pair, papers in all_pp.items()}
            result[cond] = _oo_for_papers(cond_pp)

    return result


def _oo_threshold_tab_html(oo_data: dict) -> str:
    """
    Build the "Edge Statistics" info-panel tab: left column explains the OO
    coefficient and kneedle threshold, right column shows the OO distribution
    and elbow threshold for all conditions at once (compact, hi-DPI canvases).
    """
    import json as _json

    condition_labels = {"all": "No Condition", "healthy": "Healthy", "obese": "Obese"}
    colors = {"all": "#a78bfa", "healthy": "#22d3ee", "obese": "#fbbf24"}

    cards = ""
    chart_data = {}
    for cond in ("all", "healthy", "obese"):
        if cond not in oo_data:
            continue
        d = oo_data[cond]
        oo_list = d["pairs"]
        thresh  = d["threshold"]
        n_edges = d["n_edges"]
        n_total = len(oo_list)
        n_pos   = sum(1 for _, _, v in oo_list if v > 0)
        col     = colors[cond]
        chart_data[cond] = {
            "pairs": [{"label": f"{o1[:3]}–{o2[:3]}", "oo": round(oo, 5)}
                      for o1, o2, oo in oo_list],
            "threshold": round(thresh, 6),
            "color": col,
        }
        cards += f"""
      <div class="oo-card" style="border-left:3px solid {col}">
        <div class="oo-card-head">
          <span class="oo-card-title" style="color:{col}">{condition_labels[cond]}</span>
          <span class="oo-chip"><b>{n_total}</b> pairs</span>
          <span class="oo-chip"><b>{n_pos}</b> with &ge;1 paper</span>
          <span class="oo-chip"><b>{thresh:.4f}</b> elbow (OO)</span>
          <span class="oo-chip"><b>{n_edges}</b> edges &ge; elbow</span>
        </div>
        <div class="oo-chart-wrap"><canvas id="oo-chart-{cond}" class="oo-canvas"></canvas></div>
      </div>
"""

    chart_json = _json.dumps(chart_data)

    return f"""
    <style>
      .oo-split {{
        display: grid;
        grid-template-columns: minmax(260px, 5fr) minmax(320px, 7fr);
        gap: 28px;
        align-items: start;
      }}
      @media (max-width: 900px) {{ .oo-split {{ grid-template-columns: 1fr; }} }}
      .oo-right {{ display: flex; flex-direction: column; gap: 10px; min-width: 0; height: calc(100vh - 175px); min-height: 330px; }}
      .oo-card {{
        flex: 1 1 0; min-height: 0; display: flex; flex-direction: column;
        background: #0f172a; border: 1px solid #334155; border-radius: 8px;
        padding: 8px 10px 6px;
      }}
      .oo-card-head {{
        display: flex; flex-wrap: wrap; align-items: baseline; gap: 4px 12px;
        margin-bottom: 4px;
      }}
      .oo-card-title {{ font-size: 13px; font-weight: 700; margin-right: 4px; }}
      .oo-chip {{ font-size: 11px; color: #94a3b8; white-space: nowrap; }}
      .oo-chip b {{ color: #e2e8f0; font-weight: 600; }}
      .oo-chart-wrap {{ flex: 1 1 0; min-height: 0; position: relative; }}
      .oo-canvas {{ position: absolute; inset: 0; width: 100%; height: 100%; display: block; }}
    </style>
    <div class="oo-split">
      <div class="oo-left">
        <div class="info-h2">Otsuka–Ochiai Coefficient</div>
        <p class="info-p">
          For each organ pair (A, B) the <strong>Otsuka–Ochiai (OO) coefficient</strong>
          normalises the raw co-occurrence paper count by the geometric mean of each
          organ's total paper pool, correcting for organs with a large literature base
          that would otherwise dominate any raw co-mention count:
        </p>
        <div class="info-code">OO(A, B) = n_papers(A↔B) / √( n_papers(A) × n_papers(B) )</div>
        <p class="info-p">
          <em>n_papers(A)</em> is the union of unique PubMed IDs across all pairs
          involving organ A (same same-sentence-filtered paper set used throughout).
          The coefficient lies in [0, 1] and is independent of organ-literature size.
        </p>
        <div class="info-h2">Threshold Selection (Kneedle / Elbow)</div>
        <p class="info-p">
          OO values across all pairs are sorted in ascending order and the
          <strong>kneedle algorithm</strong> (rolling-mean smoothing → log-transform →
          normalise → maximise diagonal distance) identifies the "elbow" — the point
          of maximum curvature where adding more edges yields diminishing evidence
          quality. Pairs with OO ≥ elbow are considered evidence-supported connections.
        </p>
        <p class="info-p" style="color:#94a3b8">
          Charts on the right show the OO distribution (bars, sorted descending) and the
          elbow threshold (red dashed line) for each condition. Coloured bars are
          included as network edges; grey bars fall below the threshold. All three
          charts share one y-axis so they can be compared directly. Hover a bar for
          the pair label.
        </p>
      </div>
      <div class="oo-right">
        {cards}
      </div>
    </div>
    <script>
    (function() {{
      const DATA = {chart_json};
      const sharedMax = Math.max(0.01, ...Object.values(DATA).map(
        d => Math.max(...d.pairs.map(p => p.oo), d.threshold * 1.1)));

      function draw(cond) {{
        const D = DATA[cond];
        const canvas = document.getElementById('oo-chart-' + cond);
        if (!canvas) return;
        const cssW = canvas.clientWidth, cssH = canvas.clientHeight;
        if (!cssW || !cssH) return;
        const dpr = window.devicePixelRatio || 1;
        canvas.width  = Math.round(cssW * dpr);
        canvas.height = Math.round(cssH * dpr);
        const ctx = canvas.getContext('2d');
        ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
        ctx.clearRect(0, 0, cssW, cssH);

        const pad = {{l: 40, r: 8, t: 12, b: 20}};
        const plotW = cssW - pad.l - pad.r, plotH = cssH - pad.t - pad.b;
        const pairs = D.pairs, thresh = D.threshold;
        const slot = plotW / pairs.length;
        const bw = Math.max(1, slot * 0.78);

        ctx.font = '9px system-ui, sans-serif';
        ctx.textBaseline = 'middle';
        ctx.textAlign = 'right';
        ctx.lineWidth = 1;
        [0, 0.5, 1].forEach(f => {{
          const y = Math.round(pad.t + plotH - f * plotH) + 0.5;
          ctx.strokeStyle = '#1e293b';
          ctx.beginPath(); ctx.moveTo(pad.l, y); ctx.lineTo(cssW - pad.r, y); ctx.stroke();
          ctx.fillStyle = '#94a3b8';
          ctx.fillText((sharedMax * f).toFixed(2), pad.l - 4, y);
        }});

        pairs.forEach((p, i) => {{
          const bh = (p.oo / sharedMax) * plotH;
          ctx.fillStyle = (p.oo + 5e-6 >= thresh && p.oo > 0) ? D.color : '#475569';
          ctx.fillRect(pad.l + i * slot, pad.t + plotH - bh, bw, bh);
        }});

        ctx.strokeStyle = '#64748b';
        ctx.beginPath();
        ctx.moveTo(pad.l + 0.5, pad.t); ctx.lineTo(pad.l + 0.5, pad.t + plotH + 0.5);
        ctx.lineTo(cssW - pad.r, pad.t + plotH + 0.5); ctx.stroke();

        if (thresh > 0) {{
          const ty = pad.t + plotH - (thresh / sharedMax) * plotH;
          ctx.setLineDash([4, 3]);
          ctx.strokeStyle = '#f87171'; ctx.lineWidth = 1.2;
          ctx.beginPath(); ctx.moveTo(pad.l, ty); ctx.lineTo(cssW - pad.r, ty); ctx.stroke();
          ctx.setLineDash([]);
          ctx.fillStyle = '#f87171'; ctx.textAlign = 'right'; ctx.textBaseline = 'bottom';
          ctx.fillText('elbow ' + thresh.toFixed(4), cssW - pad.r - 2, ty - 2);
        }}

        ctx.fillStyle = '#94a3b8'; ctx.textBaseline = 'top'; ctx.textAlign = 'left';
        ctx.fillText('organ pairs sorted by OO coefficient (desc); y-axis = OO', pad.l, pad.t + plotH + 6);
      }}

      function drawAll() {{ Object.keys(DATA).forEach(draw); }}

      Object.keys(DATA).forEach(cond => {{
        const canvas = document.getElementById('oo-chart-' + cond);
        if (!canvas) return;
        canvas.addEventListener('mousemove', ev => {{
          const D = DATA[cond];
          const r = canvas.getBoundingClientRect();
          const slot = (r.width - 48) / D.pairs.length;
          const i = Math.floor((ev.clientX - r.left - 40) / slot);
          const p = D.pairs[i];
          canvas.title = p ? p.label + '  OO = ' + p.oo.toFixed(4) : '';
        }});
        if (window.ResizeObserver) new ResizeObserver(() => draw(cond)).observe(canvas);
      }});
      window.addEventListener('resize', drawAll);
      drawAll();
    }})();
    </script>
"""


def _cohort_pairs_to_csv(pairs: set, organs: set) -> str:
    """Render a set of canonical organ pairs as a 0/1 adjacency CSV."""
    organs_sorted = sorted(organs)
    lines = ["," + ",".join(organs_sorted)]
    for row in organs_sorted:
        cells = [row]
        for col in organs_sorted:
            if row == col:
                cells.append("0")
            elif (min(row, col), max(row, col)) in pairs:
                cells.append("1")
            else:
                cells.append("0")
        lines.append(",".join(cells))
    return "\n".join(lines)


def _inject_cohort_comparison(html_path: Path,
                               healthy_cohort: set, obese_cohort: set,
                               all_organs: set) -> None:
    """
    Post-process the generated HTML to add one-click cohort comparison buttons.

    Embeds the healthy- and obese-cohort adjacency matrices as JS constants so
    the user can trigger the existing upload-comparison coloring (green = SHARED,
    blue = ONLY_LIT, orange = ONLY_COHORT) without having to locate and upload a
    file manually.  Uses the DataTransfer API to synthesise a File object from the
    embedded CSV text and dispatches a 'change' event on the hidden file input --
    exactly what the manual upload path does, so no JS internals need exposing.
    """
    h_csv = _cohort_pairs_to_csv(healthy_cohort, all_organs).replace("`", "\\`")
    o_csv = _cohort_pairs_to_csv(obese_cohort,   all_organs).replace("`", "\\`")

    extra_js = f"""<script>
(function(){{
  var _HC = `{h_csv}`;
  var _OC = `{o_csv}`;
  // Use loadExternalComparison (exposed by networkBuilderUtils) so that the
  // file-upload button stays dark — it should only light up when the user
  // actually uploads their own CSV, not when clicking a preset cohort button.
  window.loadHealthyCohort = function(){{ if (window.loadExternalComparison) window.loadExternalComparison(_HC, 'healthy_cohort'); }};
  window.loadObeseCohort   = function(){{ if (window.loadExternalComparison) window.loadExternalComparison(_OC, 'obese_cohort'); }};
}})();
</script>"""

    btn_h = (
        '<button id="healthy-cohort-btn" onclick="loadHealthyCohort()" class="ctrl-btn" '
        'title="Compare current Healthy view vs the predefined healthy-cohort edges">'
        '&#9878; Healthy cohort</button>'
    )
    btn_o = (
        '<button id="obese-cohort-btn" onclick="loadObeseCohort()" class="ctrl-btn" '
        'title="Compare current Obese view vs the predefined obese-cohort edges">'
        '&#9878; Obese cohort</button>'
    )

    # "No Condition" button — resets to the full 'all' view (clears healthy/obese toggle)
    no_cond_js = """<script>
(function(){
  function _noCondBtn() {
    var btn = document.getElementById('topbar-no-condition-btn');
    if (typeof setDisplayMode === 'function') {
      setDisplayMode('all');
    }
    // Deactivate healthy/obese toggle buttons visually
    ['healthy-toggle-btn','obese-toggle-btn'].forEach(function(id){
      var el = document.getElementById(id);
      if (el) el.classList.remove('active');
    });
    if (btn) btn.classList.add('active');
  }
  // Hook into the existing toggle buttons to deactivate ours when they activate
  function _patchToggleButtons() {
    ['healthy-toggle-btn','obese-toggle-btn'].forEach(function(id){
      var el = document.getElementById(id);
      if (!el || el.__ncPatched) return;
      el.__ncPatched = true;
      el.addEventListener('click', function(){
        var ncBtn = document.getElementById('topbar-no-condition-btn');
        if (ncBtn) ncBtn.classList.remove('active');
      });
    });
  }
  document.addEventListener('DOMContentLoaded', function(){
    _patchToggleButtons();
    // Mark No Condition as active by default (page starts in 'all' mode)
    var btn = document.getElementById('topbar-no-condition-btn');
    if (btn) btn.classList.add('active');
  });
  window._noCondBtn = _noCondBtn;
})();
</script>"""

    # CSS patch: style the No Condition active state
    css_patch = """<style>
#topbar-no-condition-btn.active {
  background: #4f46e5 !important;
  border-color: #4f46e5 !important;
  color: #fff !important;
}
</style>"""

    no_cond_btn = (
        '<button id="topbar-no-condition-btn" class="ctrl-btn" onclick="_noCondBtn()" '
        'title="Show all organ connections (no condition filter)">No Condition</button>'
    )

    html = html_path.read_text(encoding="utf-8")

    # 1. Inject "No Condition" button BEFORE the Healthy toggle button in the comparison-box.
    healthy_btn_anchor = 'id="healthy-toggle-btn"'
    if healthy_btn_anchor in html:
        pos = html.index(healthy_btn_anchor)
        btn_start = html.rindex('<button', 0, pos)
        html = html[:btn_start] + no_cond_btn + html[btn_start:]

    # 2. Insert the two cohort comparison buttons after the upload/clear button (inside comparison-box).
    clear_anchor  = 'id="upload-comparison-clear-btn"'
    upload_anchor = 'id="upload-comparison-btn"'
    if clear_anchor in html:
        pos      = html.index(clear_anchor)
        btn_end  = html.index('</button>', pos) + len('</button>')
        html = html[:btn_end] + btn_h + btn_o + html[btn_end:]
    elif upload_anchor in html:
        pos     = html.index(upload_anchor)
        btn_end = html.index('</button>', pos) + len('</button>')
        html = html[:btn_end] + btn_h + btn_o + html[btn_end:]
    else:
        html = html.replace('</div>\n  </div>\n</div>', btn_h + btn_o + '\n  </div>\n  </div>\n</div>', 1)

    # 3. Inject JS and CSS before </body>
    html = html.replace('</body>',
                        css_patch + no_cond_js + extra_js + '\n</body>', 1)
    html_path.write_text(html, encoding="utf-8")


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

        # Cohort CSVs are used ONLY for the comparison toggle buttons --
        # they no longer define which nodes/edges the Healthy/Obese display
        # toggle shows (that is now determined by which pairs have ≥1
        # condition-filtered paper in the literature results).
        healthy_cohort = set(load_edge_filter(HEALTHY_FILTER_CSV)) if HEALTHY_FILTER_CSV.exists() else set()
        obese_cohort   = set(load_edge_filter(OBESE_FILTER_CSV))   if OBESE_FILTER_CSV.exists() else set()

        llm_types   = load_connection_type_classifications(OUTPUT_LLM_TYPES)
        llm_descs   = load_llm_descriptions(LLM_DESCRIPTIONS)
        organ_descs = load_organ_descriptions(ORGAN_DESCRIPTIONS)

        results_by_pair = _index_results_by_pair(results)

        # Pre-compile condition keyword patterns for on-the-fly filtering.
        cond_pats = {
            cond: _compile_patterns(CONDITION_KEYWORDS[cond])
            for cond in ("healthy", "obese")
        }

        # Load condition-specific LLM data (built by --condition runs; typically
        # only covers the cohort pairs, but that's fine -- we fall back to the
        # base "all" LLM data for non-cohort pairs in condition_variants).
        cond_llm: dict[str, tuple] = {}
        for cond in ("healthy", "obese"):
            paths = condition_output_paths(cond)
            c_types = (load_connection_type_classifications(paths["llm_types"])
                       if paths["llm_types"].exists() else {})
            c_descs = (load_llm_descriptions(paths["llm_descs"])
                       if paths["llm_descs"].exists() else {})
            cond_llm[cond] = (c_types, c_descs)

        # Determine literature-defined condition edges: a pair is "healthy" (or
        # "obese") if at least one of its already-fetched papers mentions that
        # condition's keywords anywhere in title/abstract.  This is the same
        # filter applied by _filter_by_condition / build_condition_results but
        # evaluated on-the-fly here so it covers ALL 91 pairs, not just the
        # cohort subset.
        lit_healthy: set[tuple[str, str]] = set()
        lit_obese:   set[tuple[str, str]] = set()
        for pair in allowed_pairs:
            papers = results_by_pair.get(pair, {}).get("papers", [])
            if _filter_by_condition(papers, cond_pats["healthy"]):
                lit_healthy.add(pair)
            if _filter_by_condition(papers, cond_pats["obese"]):
                lit_obese.add(pair)

        lit_healthy_organs = {o for p in lit_healthy for o in p}
        lit_obese_organs   = {o for p in lit_obese   for o in p}

        G = nx.Graph()
        for organ in sorted(allowed_organs):
            organ_entry = organ_descs.get(organ, {})
            G.add_node(organ,
                       llm_description=organ_entry.get("description", ""),
                       llm_papers=organ_entry.get("papers", []),
                       is_healthy=organ in lit_healthy_organs,
                       is_obese=organ in lit_obese_organs)

        for (o1, o2) in allowed_pairs:
            G.add_edge(o1, o2)
            is_h = (o1, o2) in lit_healthy
            is_o = (o1, o2) in lit_obese
            G.edges[o1, o2]['is_healthy'] = is_h
            G.edges[o1, o2]['is_obese']   = is_o

            edge_data = results_by_pair.get((o1, o2), {})
            G.edges[o1, o2]['merged_data'] = _compute_edge_display_data(
                o1, o2, edge_data, llm_types, llm_descs,
                get_connection_types, CONNECTION_TYPES,
            )
            G.edges[o1, o2]['color'] = "#64748b"

            # Condition variants: computed on-the-fly for every pair that has
            # ≥1 condition-filtered paper (not just cohort pairs).  LLM
            # classification and descriptions use the condition-specific files
            # if they exist for this pair, or fall back to the base "all" LLM
            # data so non-cohort pairs still get a meaningful label.
            variants: dict = {}
            for cond in ("healthy", "obese"):
                cond_papers = _filter_by_condition(
                    edge_data.get("papers", []), cond_pats[cond]
                )
                if not cond_papers:
                    continue
                cond_kp = extract_key_players(cond_papers)
                cond_edge_data = {
                    **edge_data,
                    "papers":         cond_papers,
                    "n_papers_found": len(cond_papers),
                    "key_players":    cond_kp,
                }
                c_types, c_descs = cond_llm[cond]
                variants[cond] = _compute_edge_display_data(
                    o1, o2, cond_edge_data, c_types, c_descs,
                    get_connection_types, CONNECTION_TYPES,
                )
            if variants:
                G.edges[o1, o2]['condition_variants'] = variants

        # Build condition data for the Literature Statistics tab (Sankey diagrams
        # under the All / Healthy / Obese selector).  Always available now since
        # we compute it on-the-fly from the base results.
        info_condition_data: dict = {}
        for cond, cond_lit_pairs in (("healthy", lit_healthy), ("obese", lit_obese)):
            cond_rbp: dict = {}
            for pair in cond_lit_pairs:
                base   = results_by_pair.get(pair, {})
                papers = _filter_by_condition(base.get("papers", []), cond_pats[cond])
                cond_rbp[pair] = {**base, "papers": papers, "n_papers_found": len(papers)}
            cond_organs = sorted({o for p in cond_lit_pairs for o in p})
            info_condition_data[cond] = (cond_organs, cond_rbp, cond_lit_pairs)

        # Compute OO coefficient distributions for the Edge Statistics tab
        oo_data = _compute_oo_data(results_by_pair, allowed_pairs, cond_pats)

        # Store per-edge OO coefficients as graph edge attributes so
        # networkBuilderUtils can scale edge widths by evidence strength.
        _oo_lookup: dict[tuple, dict] = {}
        for cond, cdata in oo_data.items():
            for (o1, o2, oo_val) in cdata.get("pairs", []):
                key = (o1, o2) if G.has_edge(o1, o2) else (o2, o1)
                _oo_lookup.setdefault(key, {})[cond] = oo_val
        for pair in allowed_pairs:
            o1, o2 = pair
            if G.has_edge(o1, o2):
                vals = _oo_lookup.get((o1, o2), _oo_lookup.get((o2, o1), {}))
                G.edges[o1, o2]['oo_coeff']         = vals.get("all",     0.0)
                G.edges[o1, o2]['oo_coeff_healthy'] = vals.get("healthy", 0.0)
                G.edges[o1, o2]['oo_coeff_obese']   = vals.get("obese",   0.0)

        oo_thresholds = {
            cond: oo_data[cond]["threshold"]
            for cond in ("all", "healthy", "obese")
            if cond in oo_data
        }

        export_network_to_cytoscape_dashboard(
            graph=G,
            filename=str(output_path),
            include_legend=False,
            title=VIZ_TITLE,
            info_panel_tabs=_build_info_tabs(sorted(allowed_organs), results_by_pair, allowed_pairs,
                                             condition_data=info_condition_data,
                                             oo_data=oo_data),
            healthy_toggle_label="Healthy",
            obese_toggle_label="Obese",
            comparison_upload_label="Upload comparison network",
            oo_thresholds=oo_thresholds,
        )

        # Inject one-click cohort comparison buttons (⚖ Healthy cohort / ⚖ Obese cohort)
        # into the topbar so the user can compare the literature-defined condition edges
        # against the predefined cohort edges without having to upload a file manually.
        _inject_cohort_comparison(output_path, healthy_cohort, obese_cohort, allowed_organs)

        # UI refinements (compact layout, reference-network mode with blue/red signed
        # partial-correlation edges, r badge) and the weighted PC cohort matrices:
        # see Visualisation/ui_patches.py
        from Visualisation.ui_patches import apply_ui_patches, inject_weighted_cohort_matrices
        apply_ui_patches(output_path)
        inject_weighted_cohort_matrices(output_path, OUTPUT_DIR)

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
    print("     LLM types/descriptions are now saved for the cohort pairs.")
    print("     Run --viz-only to rebuild the dashboard with this condition's LLM data.")


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
