"""
Run the full pipeline in one command:
  1. PubMed literature search for curated edges (resumable)
  2. LLM connection descriptions (requires ANTHROPIC_API_KEY)
  3. Reference network visualization
  4. General organ-axis network (all-pairs PubMed search, evidence-based edges)
  5. Comparison network visualization (if --input is provided)

Usage
-----
  uv run python run_all.py
  uv run python run_all.py --skip-search --skip-llm   # just rebuild visualizations
  uv run python run_all.py --reset-search              # wipe curated search cache
  uv run python run_all.py --reset-general             # wipe general-axis cache
  uv run python run_all.py --skip-general              # skip general-axis step
  uv run python run_all.py --min-papers 3              # stricter threshold for general network
"""

import sys
import argparse
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))

ORGAN_DATA              = HERE / "metabolic_data" / "organ_data.csv"
CONNECTION_DATA         = HERE / "metabolic_data" / "connection_data.csv"
LITERATURE_RESULTS      = HERE / "metabolic_data" / "literature_results.json"
LLM_DESCRIPTIONS        = HERE / "metabolic_data" / "llm_descriptions.json"
ORGAN_SEARCH_RESULTS    = HERE / "metabolic_data" / "organ_search_results.json"
ORGAN_LLM_DESCRIPTIONS  = HERE / "metabolic_data" / "organ_descriptions.json"
GENERAL_AXIS_RESULTS    = HERE / "metabolic_data" / "general_axis_results.json"
REFERENCE_HTML          = HERE / "metabolic_data" / "reference_network.html"
GENERAL_AXIS_HTML       = HERE / "metabolic_data" / "general_organ_axis_network.html"


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
        max_results_per_pair=200,
        years_back=10,
        min_papers=5,
        delay=0.4,
        resume=not reset,
        force_research_empty=force_empty,
    )


def run_organ_llm(reset: bool, model: str) -> None:
    from Data_Loader.load_data import load_node_metadata_from_csv
    from Literature_Search.organ_descriptions import (
        run_organ_search, generate_organ_llm_descriptions,
    )

    node_metadata = load_node_metadata_from_csv(str(ORGAN_DATA))
    organs = list(node_metadata.keys())

    search_results = run_organ_search(
        organs=organs,
        output_path=ORGAN_SEARCH_RESULTS,
        resume=not reset,
        reset=reset,
    )
    generate_organ_llm_descriptions(
        organs=organs,
        search_results=search_results,
        output_path=ORGAN_LLM_DESCRIPTIONS,
        model=model,
        resume=not reset,
        reset=reset,
    )


def run_llm(reset: bool, model: str) -> None:
    from Data_Loader.load_data import load_edge_metadata_from_csv
    from Literature_Search.pubmed_search import load_literature_results
    from Literature_Search.llm_descriptions import generate_llm_descriptions

    edge_metadata = load_edge_metadata_from_csv(str(CONNECTION_DATA))
    lit_results   = load_literature_results(LITERATURE_RESULTS)

    pairs = sorted({
        (o1, o2) if o1 < o2 else (o2, o1)
        for (o1, o2) in edge_metadata
    })
    print(f"[i] {len(pairs)} organ-organ pairs to summarise.")

    generate_llm_descriptions(
        organ_pairs=pairs,
        literature_results=lit_results,
        output_path=LLM_DESCRIPTIONS,
        model=model,
        resume=not reset,
        reset=reset,
    )


def run_reference_viz() -> None:
    import networkx as nx
    from Data_Loader.load_data import load_node_metadata_from_csv, load_edge_metadata_from_csv
    from Literature_Search.pubmed_search import (
        load_literature_results, merge_with_edge_metadata, export_vocabulary_to_excel,
    )
    from Literature_Search.llm_descriptions import load_llm_descriptions
    from Literature_Search.organ_descriptions import load_organ_descriptions
    from Visualisation.networkBuilderUtils import export_network_to_cytoscape_dashboard

    export_vocabulary_to_excel(HERE / "metabolic_data" / "key_player_vocabulary.xlsx")

    node_metadata  = load_node_metadata_from_csv(str(ORGAN_DATA))
    edge_metadata  = load_edge_metadata_from_csv(str(CONNECTION_DATA))
    lit_results    = load_literature_results(LITERATURE_RESULTS)
    llm_descs      = load_llm_descriptions(LLM_DESCRIPTIONS)
    organ_descs    = load_organ_descriptions(ORGAN_LLM_DESCRIPTIONS)
    merged         = merge_with_edge_metadata(edge_metadata, lit_results, llm_descs)

    G = nx.Graph()
    for organ, desc in node_metadata.items():
        organ_entry = organ_descs.get(organ, {})
        G.add_node(organ,
                   description=desc,
                   llm_description=organ_entry.get("description", ""),
                   llm_papers=organ_entry.get("papers", []))
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


def run_general_axis(reset: bool, min_papers: int) -> None:
    from Data_Loader.load_data import load_node_metadata_from_csv
    from Literature_Search.general_axis_search import (
        run_general_axis_search, build_general_axis_viz,
    )

    node_metadata = load_node_metadata_from_csv(str(ORGAN_DATA))
    organs = list(node_metadata.keys())

    results = run_general_axis_search(
        organs=organs,
        output_path=GENERAL_AXIS_RESULTS,
        min_papers=min_papers,
        resume=not reset,
        reset=reset,
    )

    build_general_axis_viz(
        results=results,
        node_metadata=node_metadata,
        output_html=GENERAL_AXIS_HTML,
        min_papers=min_papers,
    )


def run_comparison(input_path: Path, threshold: float) -> None:
    from Matrix_Comparison.Comparison import run_network_comparison
    run_network_comparison(
        given_path=str(input_path),
        organ_data=str(ORGAN_DATA),
        connection_data=str(CONNECTION_DATA),
        threshold=threshold,
        literature_results_path=str(LITERATURE_RESULTS),
        llm_descriptions_path=str(LLM_DESCRIPTIONS),
        organ_descriptions_path=str(ORGAN_LLM_DESCRIPTIONS),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Full metabolic network pipeline.")
    parser.add_argument("--input",        default=None,
                        help="Path to comparison network CSV (optional).")
    parser.add_argument("--threshold",    type=float, default=0.3,
                        help="Edge threshold for comparison (default 0.3).")
    parser.add_argument("--skip-search",  action="store_true",
                        help="Skip the PubMed literature search step.")
    parser.add_argument("--skip-llm",        action="store_true",
                        help="Skip the LLM edge description step.")
    parser.add_argument("--skip-organ-llm",  action="store_true",
                        help="Skip the LLM organ description step.")
    parser.add_argument("--skip-general",    action="store_true",
                        help="Skip the general organ-axis network step.")
    parser.add_argument("--force-empty",     action="store_true",
                        help="Re-search edges that previously returned 0 papers.")
    parser.add_argument("--reset-search",    action="store_true",
                        help="Delete literature cache and search everything from scratch.")
    parser.add_argument("--reset-llm",       action="store_true",
                        help="Delete edge LLM cache and regenerate.")
    parser.add_argument("--reset-organ-llm", action="store_true",
                        help="Delete organ LLM cache and regenerate.")
    parser.add_argument("--llm-model",       default="north-mini-code-1.0",
                        help="Ollama model for LLM summaries (default: north-mini-code-1.0).")
    parser.add_argument("--reset-general",   action="store_true",
                        help="Delete general-axis cache and re-search all pairs.")
    parser.add_argument("--min-papers",      type=int, default=10,
                        help="Min papers required for an edge in the general network (default 10).")
    args = parser.parse_args()

    # ── Step 1: Curated literature search ─────────────────────────────────
    if not args.skip_search:
        step("Step 1 / 6  —  PubMed literature search (curated edges)")
        run_search(force_empty=args.force_empty, reset=args.reset_search)
    else:
        print("\n[i] Skipping literature search (--skip-search).")

    # ── Step 2: LLM edge descriptions ─────────────────────────────────────
    if not args.skip_llm:
        step("Step 2 / 6  —  LLM edge summaries")
        run_llm(reset=args.reset_llm, model=args.llm_model)
    else:
        print("\n[i] Skipping edge LLM descriptions (--skip-llm).")

    # ── Step 3: LLM organ descriptions ────────────────────────────────────
    if not args.skip_organ_llm:
        step("Step 3 / 6  —  LLM organ descriptions")
        run_organ_llm(reset=args.reset_organ_llm, model=args.llm_model)
    else:
        print("\n[i] Skipping organ LLM descriptions (--skip-organ-llm).")

    # ── Step 4: Reference network ──────────────────────────────────────────
    step("Step 4 / 6  —  Reference network visualization")
    run_reference_viz()

    # ── Step 5: General organ-axis network ────────────────────────────────
    if not args.skip_general:
        step("Step 5 / 6  —  General organ-axis network (all-pairs search)")
        run_general_axis(reset=args.reset_general, min_papers=args.min_papers)
    else:
        print("\n[i] Skipping general organ-axis network (--skip-general).")

    # ── Step 6: Comparison network ─────────────────────────────────────────
    step("Step 6 / 6  —  Comparison network visualization")
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
    print(f"  Reference network    : {REFERENCE_HTML}")
    print(f"  General axis network : {GENERAL_AXIS_HTML}")
    if args.input:
        input_path = Path(args.input)
        print(f"  Comparison network   : {input_path.parent / (input_path.stem + '_comparison.html')}")


if __name__ == "__main__":
    main()
