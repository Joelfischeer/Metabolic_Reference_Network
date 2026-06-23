"""
Run the full pipeline in one command:
  1. PubMed literature search (resumable)
  2. Reference network visualization
  3. Comparison network visualization (if --input is provided)

Usage
-----
  uv run python run_all.py
  uv run python run_all.py --input ../metabolic_network.csv --threshold 0.3
  uv run python run_all.py --force-empty          # re-search edges with 0 papers
  uv run python run_all.py --reset-search         # wipe literature cache and re-search
  uv run python run_all.py --skip-search          # skip PubMed, just rebuild visualizations
"""

import sys
import argparse
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))

ORGAN_DATA          = HERE / "metabolic_data" / "organ_data.csv"
CONNECTION_DATA     = HERE / "metabolic_data" / "connection_data.csv"
LITERATURE_RESULTS  = HERE / "metabolic_data" / "literature_results.json"
LLM_DESCRIPTIONS    = HERE / "metabolic_data" / "llm_descriptions.json"
REFERENCE_HTML      = HERE / "metabolic_data" / "reference_network.html"


def step(label: str) -> None:
    width = 60
    print(f"\n{'='*width}")
    print(f"  {label}")
    print(f"{'='*width}")


def run_search(force_empty: bool, reset: bool) -> None:
    from Data_Loader.load_data import load_edge_metadata_from_csv
    from Literature_Search.pubmed_search import run_literature_search

    if reset and LITERATURE_RESULTS.exists():
        LITERATURE_RESULTS.unlink()
        print(f"[i] Cache deleted: {LITERATURE_RESULTS.name}")

    edge_metadata = load_edge_metadata_from_csv(str(CONNECTION_DATA))
    pairs = sorted({
        (o1, o2) if o1 < o2 else (o2, o1)
        for (o1, o2) in edge_metadata
    })
    print(f"[i] {len(pairs)} organ-organ connections to search.")

    run_literature_search(
        organ_pairs=pairs,
        output_path=LITERATURE_RESULTS,
        max_results_per_pair=25,
        years_back=10,
        min_papers=5,
        delay=0.4,
        resume=not reset,
        force_research_empty=force_empty,
    )


def run_llm(reset: bool) -> None:
    from Data_Loader.load_data import load_edge_metadata_from_csv
    from Literature_Search.pubmed_search import load_literature_results
    from Literature_Search.llm_descriptions import generate_llm_descriptions

    edge_metadata = load_edge_metadata_from_csv(str(CONNECTION_DATA))
    lit_results   = load_literature_results(LITERATURE_RESULTS)

    pairs = sorted({
        (o1, o2) if o1 < o2 else (o2, o1)
        for (o1, o2) in edge_metadata
    })
    print(f"[i] {len(pairs)} organ-organ pairs to describe.")

    generate_llm_descriptions(
        organ_pairs=pairs,
        edge_metadata=edge_metadata,
        literature_results=lit_results,
        output_path=LLM_DESCRIPTIONS,
        resume=not reset,
        reset=reset,
    )


def run_reference_viz() -> None:
    import networkx as nx
    from Data_Loader.load_data import load_node_metadata_from_csv, load_edge_metadata_from_csv
    from Literature_Search.pubmed_search import load_literature_results, merge_with_edge_metadata
    from Literature_Search.llm_descriptions import load_llm_descriptions
    from Visualisation.networkBuilderUtils import export_network_to_cytoscape_dashboard

    node_metadata = load_node_metadata_from_csv(str(ORGAN_DATA))
    edge_metadata = load_edge_metadata_from_csv(str(CONNECTION_DATA))
    lit_results   = load_literature_results(LITERATURE_RESULTS)
    llm_descs     = load_llm_descriptions(LLM_DESCRIPTIONS)
    merged        = merge_with_edge_metadata(edge_metadata, lit_results, llm_descs)

    G = nx.Graph()
    for organ, desc in node_metadata.items():
        G.add_node(organ, description=desc)
    for (o1, o2), text in edge_metadata.items():
        if o1 >= o2:
            continue
        G.add_edge(o1, o2)
        G.edges[o1, o2]['description'] = text
        G.edges[o1, o2]['merged_data'] = merged.get((o1, o2), {})
        G.edges[o1, o2]['color'] = "#64748b"

    export_network_to_cytoscape_dashboard(
        graph=G,
        filename=str(REFERENCE_HTML),
        include_legend=False,
        title="Metabolic Reference Network",
    )
    print(f"[ok] Reference network: {REFERENCE_HTML}")


def run_comparison(input_path: Path, threshold: float) -> None:
    from Matrix_Comparison.Comparison import run_network_comparison
    run_network_comparison(
        given_path=str(input_path),
        organ_data=str(ORGAN_DATA),
        connection_data=str(CONNECTION_DATA),
        threshold=threshold,
        literature_results_path=str(LITERATURE_RESULTS),
        llm_descriptions_path=str(LLM_DESCRIPTIONS),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Full metabolic network pipeline.")
    parser.add_argument("--input",        default=None,
                        help="Path to comparison network CSV (optional).")
    parser.add_argument("--threshold",    type=float, default=0.3,
                        help="Edge threshold for comparison (default 0.3).")
    parser.add_argument("--skip-search",  action="store_true",
                        help="Skip the PubMed literature search step.")
    parser.add_argument("--skip-llm",    action="store_true",
                        help="Skip the LLM description generation step.")
    parser.add_argument("--force-empty",  action="store_true",
                        help="Re-search edges that previously returned 0 papers.")
    parser.add_argument("--reset-search", action="store_true",
                        help="Delete literature cache and search everything from scratch.")
    parser.add_argument("--reset-llm",   action="store_true",
                        help="Delete LLM description cache and regenerate everything.")
    args = parser.parse_args()

    # ── Step 1: Literature search ──────────────────────────────────────────
    if not args.skip_search:
        step("Step 1 / 4  —  PubMed literature search")
        run_search(force_empty=args.force_empty, reset=args.reset_search)
    else:
        print("\n[i] Skipping literature search (--skip-search).")

    # ── Step 2: LLM descriptions ───────────────────────────────────────────
    if not args.skip_llm:
        step("Step 2 / 4  —  LLM connection descriptions")
        run_llm(reset=args.reset_llm)
    else:
        print("\n[i] Skipping LLM descriptions (--skip-llm).")

    # ── Step 3: Reference network ──────────────────────────────────────────
    step("Step 3 / 4  —  Reference network visualization")
    run_reference_viz()

    # ── Step 4: Comparison network ─────────────────────────────────────────
    step("Step 4 / 4  —  Comparison network visualization")
    if args.input:
        input_path = Path(args.input)
        if not input_path.exists():
            print(f"[!] Input file not found: {input_path}")
            sys.exit(1)
        run_comparison(input_path, args.threshold)
    else:
        print("[i] No --input provided — skipping comparison network.")
        print("    Run with:  uv run python run_all.py --input path/to/network.csv")

    step("Done")
    print(f"  Reference network : {REFERENCE_HTML}")
    if args.input:
        input_path = Path(args.input)
        print(f"  Comparison network: {input_path.parent / (input_path.stem + '_comparison.html')}")


if __name__ == "__main__":
    main()
