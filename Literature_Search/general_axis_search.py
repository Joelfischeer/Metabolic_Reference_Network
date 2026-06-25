"""
Data-driven organ-axis network: search ALL possible organ-organ pairs in PubMed
(last 5 years), then build a network from every pair with ≥ 2 papers found.

No curated edge list is assumed — any pair that has literature support gets an edge.

Results are saved incrementally so the search can be resumed if interrupted.

Usage (standalone):
    uv run python Literature_Search/general_axis_search.py
    uv run python Literature_Search/general_axis_search.py --reset
    uv run python Literature_Search/general_axis_search.py --min-papers 3
"""

import json
import sys
import argparse
import itertools
from pathlib import Path

HERE = Path(__file__).parent.parent
sys.path.insert(0, str(HERE))

DEFAULT_OUTPUT_SEARCH = HERE / "metabolic_data" / "general_axis_results.json"
DEFAULT_OUTPUT_HTML   = HERE / "metabolic_data" / "general_organ_axis_network.html"
DEFAULT_MIN_PAPERS    = 10
DEFAULT_YEARS_BACK    = 5
DEFAULT_MAX_RESULTS   = 200
DEFAULT_DELAY         = 0.4


def all_organ_pairs(organs: list[str]) -> list[tuple[str, str]]:
    """Return all unique sorted pairs of organs."""
    return sorted(
        (a, b) if a < b else (b, a)
        for a, b in itertools.combinations(organs, 2)
    )


def run_general_axis_search(
    organs: list[str],
    output_path: "str | Path" = DEFAULT_OUTPUT_SEARCH,
    min_papers: int = DEFAULT_MIN_PAPERS,
    years_back: int = DEFAULT_YEARS_BACK,
    max_results_per_pair: int = DEFAULT_MAX_RESULTS,
    delay: float = DEFAULT_DELAY,
    resume: bool = True,
    reset: bool = False,
) -> dict:
    """
    Search all organ-organ pairs and save results incrementally.

    Returns the full results dict keyed by "organ1|organ2".
    """
    from Literature_Search.pubmed_search import (
        search_with_cascade, fetch_abstracts,
        extract_key_players, _infer_connection_type,
    )
    from datetime import datetime

    output_path = Path(output_path)

    if reset and output_path.exists():
        output_path.unlink()
        print(f"[i] Cache deleted: {output_path.name}")

    results: dict = {}
    if resume and output_path.exists():
        with open(output_path, encoding="utf-8") as f:
            results = json.load(f)
        n_cached = len(results)
        print(f"[i] Resuming: {n_cached} pairs already searched.")

    pairs = all_organ_pairs(organs)
    total = len(pairs)
    print(f"[i] {total} organ-organ pairs to search ({len(organs)} organs).\n")

    for idx, (o1, o2) in enumerate(pairs):
        key = f"{o1}|{o2}"
        sym = f"{o2}|{o1}"

        if key in results or sym in results:
            existing = results.get(key) or results.get(sym, {})
            n = existing.get("n_papers_found", 0)
            print(f"  [{idx+1}/{total}] Cached ({n} papers): {o1} <-> {o2}")
            continue

        print(f"\n  [{idx+1}/{total}] Searching: {o1} <-> {o2}")
        import time
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
        print(
            f"    => {n} papers | strategy: {strategy}"
            f" | hormones: {len(key_players['hormones'])}"
            f" | metabolites: {len(key_players['metabolites'])}"
            f" | proteins: {len(key_players['proteins'])}"
        )

        results[key] = {
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

    n_with_papers = sum(1 for v in results.values() if v.get("n_papers_found", 0) >= min_papers)
    print(
        f"\n[ok] Done. {len(results)} pairs searched | "
        f"{n_with_papers} pairs with ≥{min_papers} papers → will become edges."
    )
    print(f"     Saved to: {output_path}")
    return results


def build_general_axis_viz(
    results: dict,
    node_metadata: dict,
    output_html: "str | Path" = DEFAULT_OUTPUT_HTML,
    min_papers: int = DEFAULT_MIN_PAPERS,
) -> None:
    """
    Build and export the general organ-axis network HTML.

    Only pairs with ≥ min_papers become edges. Node colour, sidebar, search,
    and layout controls are identical to the reference network.
    """
    import networkx as nx
    from Visualisation.networkBuilderUtils import (
        export_network_to_cytoscape_dashboard, ORGAN_COLORS, DEFAULT_NODE_COLOR,
    )

    G = nx.Graph()

    # Add all organs as nodes
    for organ, desc in node_metadata.items():
        G.add_node(organ, description=desc)

    # Add edges for pairs with enough literature support
    for key, data in results.items():
        o1 = data.get("organ1", "")
        o2 = data.get("organ2", "")
        n  = data.get("n_papers_found", 0)

        if n < min_papers:
            continue
        if not G.has_node(o1) or not G.has_node(o2):
            continue
        if G.has_edge(o1, o2):
            continue

        kp    = data.get("key_players", {})
        papers = data.get("papers", [])

        # Build a merged_data structure matching what networkBuilderUtils expects
        merged_data = {
            "description":     "",
            "connection_type": data.get("connection_type", ""),
            "key_players_raw": [],
            "key_players_merged": {
                "hormones":    kp.get("hormones", []),
                "metabolites": kp.get("metabolites", []),
                "proteins":    kp.get("proteins", []),
            },
            "key_players_counts": {
                "hormones":    kp.get("hormones_counts", {}),
                "metabolites": kp.get("metabolites_counts", {}),
                "proteins":    kp.get("proteins_counts", {}),
            },
            "notes":  "",
            "sources": [],
            "ai_description": "",
            "pubmed": {
                "n_papers": n,
                "papers":   papers[:5],
                "query":    data.get("pubmed_query", ""),
                "strategy": data.get("strategy_used", ""),
            },
        }

        # Edge thickness: scale by number of papers (capped at 20)
        weight = min(n / 5, 4)

        G.add_edge(o1, o2)
        G.edges[o1, o2]["description"]  = ""
        G.edges[o1, o2]["merged_data"]  = merged_data
        G.edges[o1, o2]["color"]        = "#64748b"
        G.edges[o1, o2]["weight"]       = weight

    n_edges = G.number_of_edges()
    print(f"[i] General axis network: {G.number_of_nodes()} organs, {n_edges} evidence-based edges.")

    export_network_to_cytoscape_dashboard(
        graph=G,
        filename=str(output_html),
        include_legend=False,
        title=f"General Organ Axis Network (≥{min_papers} papers, last 5 years)",
    )
    print(f"[ok] General axis network: {output_html}")


def main():
    parser = argparse.ArgumentParser(
        description="Build a data-driven organ-axis network from PubMed."
    )
    parser.add_argument("--reset", action="store_true",
                        help="Delete cache and re-search all pairs from scratch.")
    parser.add_argument("--min-papers", type=int, default=DEFAULT_MIN_PAPERS,
                        help=f"Minimum papers to include an edge (default {DEFAULT_MIN_PAPERS}).")
    parser.add_argument("--years-back", type=int, default=DEFAULT_YEARS_BACK,
                        help=f"How many years of literature to search (default {DEFAULT_YEARS_BACK}).")
    args = parser.parse_args()

    from Data_Loader.load_data import load_node_metadata_from_csv
    node_metadata = load_node_metadata_from_csv(str(HERE / "metabolic_data" / "organ_data.csv"))
    organs = list(node_metadata.keys())

    results = run_general_axis_search(
        organs=organs,
        output_path=DEFAULT_OUTPUT_SEARCH,
        min_papers=args.min_papers,
        years_back=args.years_back,
        resume=not args.reset,
        reset=args.reset,
    )

    build_general_axis_viz(
        results=results,
        node_metadata=node_metadata,
        output_html=DEFAULT_OUTPUT_HTML,
        min_papers=args.min_papers,
    )


if __name__ == "__main__":
    main()
