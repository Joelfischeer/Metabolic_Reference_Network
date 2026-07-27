"""
Build the literature-derived multi-layer reference network.

  Step 1 — PubMed search: for every organ pair × 5 communication types
  Step 2 — Visualization: build the interactive multi-layer HTML

Usage
-----
  # Full run (search + viz):
  uv run python run_lit_ref.py

  # Skip search (rebuild viz from existing cache):
  uv run python run_lit_ref.py --skip-search

  # Reset search cache and start fresh:
  uv run python run_lit_ref.py --reset-search

  # Stricter edge threshold (default 5):
  uv run python run_lit_ref.py --min-papers 10
"""

import sys
import argparse
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))

ORGAN_DATA      = HERE / "metabolic_data" / "organ_data.csv"
LIT_REF_RESULTS = HERE / "metabolic_data" / "lit_ref_results.json"
LIT_REF_HTML    = HERE / "metabolic_data" / "literature_reference_network.html"


def step(label: str) -> None:
    width = 60
    print(f"\n{'='*width}\n  {label}\n{'='*width}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Literature reference network pipeline (search + viz)."
    )
    parser.add_argument("--skip-search",  action="store_true",
                        help="Skip PubMed search, use existing cache.")
    parser.add_argument("--reset-search", action="store_true",
                        help="Delete search cache and re-run all queries.")
    parser.add_argument("--min-papers",   type=int, default=5,
                        help="Min papers per communication type to draw an edge (default 5).")
    parser.add_argument("--max-results",  type=int, default=200,
                        help="Max papers fetched per query (default 200).")
    parser.add_argument("--years-back",   type=int, default=10,
                        help="PubMed look-back window in years (default 10).")
    args = parser.parse_args()

    from Data_Loader.load_data import load_node_metadata_from_csv
    node_metadata = load_node_metadata_from_csv(str(ORGAN_DATA))
    organs = list(node_metadata.keys())
    print(f"[i] {len(organs)} organs loaded.")

    # ── Step 1: Literature search ─────────────────────────────────────────────
    if not args.skip_search:
        step("Step 1 / 2  —  PubMed search (organ pairs × 5 communication types)")
        from Literature_Reference_Network.lit_ref_search import run_lit_ref_search
        run_lit_ref_search(
            organs      = organs,
            output_path = LIT_REF_RESULTS,
            min_papers  = args.min_papers,
            max_results = args.max_results,
            years_back  = args.years_back,
            resume      = not args.reset_search,
            reset       = args.reset_search,
        )
    else:
        print(f"\n[i] Skipping search (--skip-search).  Using: {LIT_REF_RESULTS.name}")

    # ── Step 2: Visualization ─────────────────────────────────────────────────
    step("Step 2 / 2  —  Build interactive HTML visualization")

    if not LIT_REF_RESULTS.exists():
        print(f"[!] Search results not found: {LIT_REF_RESULTS}")
        print("    Run without --skip-search to generate them.")
        sys.exit(1)

    import json
    with open(LIT_REF_RESULTS, encoding="utf-8") as f:
        results = json.load(f)

    from Literature_Reference_Network.lit_ref_viz import build_lit_ref_viz
    build_lit_ref_viz(
        search_results = results,
        node_metadata  = node_metadata,
        output_html    = LIT_REF_HTML,
        min_papers     = args.min_papers,
    )

    step("Done")
    print(f"  Literature reference network : {LIT_REF_HTML}")


if __name__ == "__main__":
    main()
