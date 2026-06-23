"""
Compare a given metabolic network CSV against the reference network.

Usage:
    python Analysis.py
    python Analysis.py --input path/to/metabolic_network.csv --threshold 0.3
"""

import sys
import argparse
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))

ORGAN_DATA         = HERE / "metabolic_data" / "organ_data.csv"
CONNECTION_DATA    = HERE / "metabolic_data" / "connection_data.csv"
LITERATURE_RESULTS = HERE / "metabolic_data" / "literature_results.json"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default=str(HERE.parent / "metabolic_network.csv"),
                        help="Path to the metabolic network CSV to compare against the reference.")
    parser.add_argument("--threshold", type=float, default=0.3,
                        help="Edge presence threshold (default 0.3).")
    args = parser.parse_args()

    from Matrix_Comparison.Comparison import run_network_comparison
    run_network_comparison(
        given_path=args.input,
        organ_data=str(ORGAN_DATA),
        connection_data=str(CONNECTION_DATA),
        threshold=args.threshold,
        literature_results_path=str(LITERATURE_RESULTS),
    )


if __name__ == "__main__":
    main()
