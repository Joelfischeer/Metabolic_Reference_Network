"""
reprocess_key_players.py
=========================
Re-derive key_players/connection_type for every cached pair in
reference_network_only_metabolic/metabolic_literature_results.json from the
already-cached papers — no PubMed calls. Use this after a fix to
extract_key_players()'s synonym-merging logic (or vocab) to bring an
already-fully-cached results file up to date without a full --reset
re-search (which would needlessly re-hit PubMed for every pair).

Run from the Metabolic_Reference_Network/ directory:
    uv run python reprocess_key_players.py

Then rebuild the dashboard from the updated JSON:
    uv run python run_metabolic_lit_search.py --viz-only
"""

import json
import shutil
from pathlib import Path

from Literature_Search.pubmed_search import extract_key_players

HERE = Path(__file__).parent
RESULTS_JSON = HERE / "reference_network_only_metabolic" / "metabolic_literature_results.json"


def main():
    if not RESULTS_JSON.exists():
        print(f"[!] Not found: {RESULTS_JSON}")
        return

    backup = RESULTS_JSON.with_suffix(".json.bak")
    shutil.copy2(RESULTS_JSON, backup)
    print(f"[i] Backup saved: {backup}")

    with open(RESULTS_JSON, encoding="utf-8") as f:
        results = json.load(f)

    total = len(results)
    changed = 0
    for idx, (key, entry) in enumerate(results.items()):
        papers = entry.get("papers", [])
        if not papers:
            continue
        old_counts = {
            cat: entry.get("key_players", {}).get(f"{cat}_counts", {})
            for cat in ("hormones", "metabolites", "proteins")
        }
        key_players = extract_key_players(papers)
        entry["key_players"] = key_players

        new_counts = {cat: key_players.get(f"{cat}_counts", {}) for cat in old_counts}
        if new_counts != old_counts:
            changed += 1

        if (idx + 1) % 20 == 0 or idx + 1 == total:
            print(f"  [{idx+1}/{total}] reprocessed")

    with open(RESULTS_JSON, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)

    print(f"\n[ok] Reprocessed {total} pairs ({changed} had different key-player "
          f"counts after the synonym-merge fix). Saved: {RESULTS_JSON}")


if __name__ == "__main__":
    main()
