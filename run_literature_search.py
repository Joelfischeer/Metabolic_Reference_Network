"""
Standalone script to run PubMed literature search for all organ-organ connections.

Saves results to metabolic_data/literature_results.json.
Auto-resumes from where it left off if interrupted.

Run from any directory:
    python run_literature_search.py               # first run or resume
    python run_literature_search.py --force-empty # re-search edges with 0 papers
    python run_literature_search.py --reset       # wipe cache and start over
"""

import sys
import argparse
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))

from Data_Loader.load_data import load_edge_metadata_from_csv
from Literature_Search.pubmed_search import run_literature_search

CONNECTION_DATA = HERE / "metabolic_data" / "connection_data.csv"
OUTPUT_PATH     = HERE / "metabolic_data" / "literature_results.json"

MAX_PAPERS  = 200
YEARS_BACK  = 10
MIN_PAPERS  = 5
DELAY       = 0.4


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--force-empty", action="store_true",
                        help="Re-search edges that previously returned 0 papers.")
    parser.add_argument("--reset", action="store_true",
                        help="Delete cache and start from scratch.")
    args = parser.parse_args()

    if args.reset and OUTPUT_PATH.exists():
        OUTPUT_PATH.unlink()
        print(f"[i] Cache deleted: {OUTPUT_PATH}")

    print("[i] Loading edge metadata...")
    edge_metadata = load_edge_metadata_from_csv(str(CONNECTION_DATA))

    pairs = sorted({
        (o1, o2) if o1 < o2 else (o2, o1)
        for (o1, o2) in edge_metadata
    })
    print(f"[i] {len(pairs)} organ-organ connections to search.\n")

    results = run_literature_search(
        organ_pairs=pairs,
        output_path=OUTPUT_PATH,
        max_results_per_pair=MAX_PAPERS,
        years_back=YEARS_BACK,
        min_papers=MIN_PAPERS,
        delay=DELAY,
        resume=not args.reset,
        force_research_empty=args.force_empty,
    )

    empty = [k for k, v in results.items() if v.get("n_papers_found", 0) == 0]
    if empty:
        print(f"\n[!] {len(empty)} edges still have 0 papers — try --force-empty.")
        for k in empty:
            print(f"    {k}")

    print(f"\n[ok] Run 'python build_reference_viz.py' to regenerate the visualization.")


if __name__ == "__main__":
    main()
