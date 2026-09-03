"""
run_llm_descriptions.py
========================
Generates LLM-based literature summaries for this network's organs (node
sidebar "Metabolic Summary") and organ-organ connections (edge sidebar
"Literature Summary"), using a locally running Ollama model.

Also classifies each robust pair into 1-3 connection types from
config's CONNECTION_TYPES (same categories as
reference_network_only_metabolic), via 3 independent LLM passes per pair
with majority voting — see Literature_Search/llm_connection_type.py.

No new PubMed queries are issued — everything is grounded entirely in
data already fetched for edge definition:
  - Organ summaries      <- search_results_{condition}.json, i.e. the same
    per-organ pool fetched by (MeSH_ORGAN OR aliases) AND METABOLIC_FILTER
    AND CONDITION_FILTER.
  - Connection summaries
    and connection types  <- bootstrap_results_{condition}.json's "papers"
    field per pair, which is already restricted to same-sentence
    cross-mention evidence (the union of papers that triggered a
    same-sentence hit across the bootstrap iterations — see run_network.py's
    run_bootstrap()). If a condition's bootstrap cache predates the
    same-sentence rule, rebuild it first (delete bootstrap_results_*.json
    and re-run run_network.py) before generating connection summaries here.

Requires Ollama running locally with the model pulled:
    ollama pull llama3.2

Run from the Metabolic_Reference_Network/ directory:
    uv run -m Edge_cosine_general_reference_network.run_llm_descriptions --condition healthy
    uv run -m Edge_cosine_general_reference_network.run_llm_descriptions --condition obese

    --reset   regenerate all summaries instead of resuming from cache
    --model   Ollama model to use (default: llama3.2)
"""

import sys
import json
import argparse
import importlib.util
from pathlib import Path
from itertools import combinations

HERE = Path(__file__).parent          # Edge_cosine_general_reference_network/
ROOT = HERE.parent                     # Metabolic_Reference_Network/
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from threshold_utils import Elbow
from Literature_Search.llm_descriptions import (
    generate_llm_descriptions,
    generate_llm_organ_descriptions,
    DEFAULT_MODEL,
)
from Literature_Search.llm_connection_type import generate_connection_type_classifications

CONDITION_CONFIGS = {
    "healthy": HERE / "healthy" / "config_healthy.py",
    "obese":   HERE / "obese"   / "config_obese.py",
}


def _load_config(condition: str):
    path = CONDITION_CONFIGS[condition]
    spec = importlib.util.spec_from_file_location(f"_cfg_{condition}", path)
    mod  = importlib.util.module_from_spec(spec)
    mod.Elbow = Elbow  # lets the config file write `MIN_BOOTSTRAP_MEAN = Elbow`
    spec.loader.exec_module(mod)
    return mod


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--condition", required=True,
                        choices=list(CONDITION_CONFIGS),
                        help="Which condition to generate summaries for.")
    parser.add_argument("--reset", action="store_true",
                        help="Delete cache and regenerate all summaries.")
    parser.add_argument("--model", default=DEFAULT_MODEL,
                        help=f"Ollama model to use (default: {DEFAULT_MODEL}).")
    args = parser.parse_args()

    cfg = _load_config(args.condition)
    out_dir     = HERE / cfg.CONDITION_NAME
    search_json = out_dir / f"search_results_{cfg.CONDITION_NAME}.json"
    boot_json   = out_dir / f"bootstrap_results_{cfg.CONDITION_NAME}.json"
    organ_out   = out_dir / f"llm_organ_descriptions_{cfg.CONDITION_NAME}.json"
    pair_out    = out_dir / f"llm_pair_descriptions_{cfg.CONDITION_NAME}.json"
    type_out    = out_dir / f"connection_types_{cfg.CONDITION_NAME}.json"

    if not search_json.exists() or not boot_json.exists():
        print(f"[!] Search/bootstrap results not found in {out_dir}/. "
              f"Run run_network.py for this condition first:\n"
              f"    uv run -m Edge_cosine_general_reference_network.run_network "
              f"--condition {cfg.CONDITION_NAME}")
        sys.exit(1)

    with open(search_json, encoding="utf-8") as f:
        search_results = json.load(f)
    with open(boot_json, encoding="utf-8") as f:
        boot_results = json.load(f)

    # Derive the organ set from the bootstrap results themselves (not
    # search_results.keys(), which may still hold organs outside the current
    # cohort from an earlier/broader run) — this matches exactly which
    # organs/pairs were actually bootstrapped.
    organs = sorted({v["organ1"] for v in boot_results.values()} |
                    {v["organ2"] for v in boot_results.values()})
    pairs  = list(combinations(organs, 2))

    print(f"\n[i] Condition: {cfg.VIZ_LABEL}")
    print(f"[i] {len(organs)} organs, {len(pairs)} pairs")
    print(f"[i] Model: {args.model}\n")

    print("=" * 60)
    print("Organ-level summaries (search_results — per-organ query)")
    print("=" * 60)
    generate_llm_organ_descriptions(
        organs=organs,
        search_results=search_results,
        output_path=organ_out,
        model=args.model,
        resume=not args.reset,
        reset=args.reset,
    )

    print("\n" + "=" * 60)
    print("Connection-level summaries (bootstrap results — same-sentence evidence)")
    print("=" * 60)
    literature_results = {
        key: {"papers": v.get("papers", []), "n_papers_found": v.get("n_cooccur_total", 0)}
        for key, v in boot_results.items()
    }
    generate_llm_descriptions(
        organ_pairs=pairs,
        literature_results=literature_results,
        output_path=pair_out,
        model=args.model,
        resume=not args.reset,
        reset=args.reset,
    )

    print("\n" + "=" * 60)
    print("Connection-type classification (bootstrap results — same-sentence evidence)")
    print("=" * 60)
    generate_connection_type_classifications(
        organ_pairs=pairs,
        literature_results=literature_results,
        connection_types=cfg.CONNECTION_TYPES,
        output_path=type_out,
        model=args.model,
        max_papers=cfg.LLM_MAX_PAPERS,
        resume=not args.reset,
        reset=args.reset,
    )

    print(f"\n[ok] Rebuild the dashboard to pick these up:\n"
          f"    uv run -m Edge_cosine_general_reference_network.run_network "
          f"--condition {cfg.CONDITION_NAME} --viz-only")


if __name__ == "__main__":
    main()
