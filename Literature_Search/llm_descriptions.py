"""
Generate LLM-based descriptions for each organ-organ connection pair.

Uses the Anthropic API to write 2-3 sentences explaining why a given
metabolic or hormonal connection makes biological sense.

Results are saved incrementally to metabolic_data/llm_descriptions.json
so the process can be resumed if interrupted.

Usage (standalone):
    uv run python Literature_Search/llm_descriptions.py
    uv run python Literature_Search/llm_descriptions.py --reset
"""

import json
import os
import sys
import argparse
from pathlib import Path
from datetime import datetime

HERE = Path(__file__).parent.parent
sys.path.insert(0, str(HERE))

DEFAULT_OUTPUT = HERE / "metabolic_data" / "llm_descriptions.json"
MODEL = "claude-haiku-4-5"


def _build_prompt(organ1: str, organ2: str, edge_text: str, key_players: list[str]) -> str:
    kp_str = ", ".join(key_players[:10]) if key_players else "unknown"
    context = edge_text.strip() if edge_text else "No curated notes available."
    return (
        f"You are a metabolic biology expert. "
        f"Write exactly 2-3 sentences explaining why the connection between "
        f"'{organ1}' and '{organ2}' makes metabolic or hormonal sense. "
        f"Be specific and mention key molecules or hormones where relevant. "
        f"Do not use bullet points — write flowing, informative prose.\n\n"
        f"Curated notes: {context}\n"
        f"Key players from literature: {kp_str}"
    )


def generate_llm_descriptions(
    organ_pairs: list[tuple[str, str]],
    edge_metadata: dict,
    literature_results: dict,
    output_path: "str | Path" = DEFAULT_OUTPUT,
    resume: bool = True,
    reset: bool = False,
) -> dict[str, str]:
    """
    Generate LLM descriptions for all organ pairs and save incrementally.

    Parameters
    ----------
    organ_pairs       : sorted list of (organ1, organ2) tuples
    edge_metadata     : dict keyed by (organ1, organ2) -> description text
    literature_results: raw literature JSON (keyed "organ1|organ2")
    output_path       : where to save/load llm_descriptions.json
    resume            : if True, skip pairs already in the cache file
    reset             : if True, delete the cache and regenerate everything
    """
    try:
        import anthropic
    except ImportError:
        print("[!] anthropic package not installed — skipping LLM descriptions.")
        print("    Install with: uv add anthropic")
        return descriptions

    api_key = os.environ.get("ANTHROPIC_API_KEY")
    if not api_key:
        print("[!] ANTHROPIC_API_KEY not set — skipping LLM descriptions.")
        print("    Set it with: export ANTHROPIC_API_KEY='sk-ant-...'")
        print("    Get a key at: https://console.anthropic.com")
        return descriptions

    output_path = Path(output_path)

    if reset and output_path.exists():
        output_path.unlink()
        print(f"[i] Cache deleted: {output_path.name}")

    descriptions: dict[str, str] = {}
    if resume and output_path.exists():
        with open(output_path, encoding="utf-8") as f:
            descriptions = json.load(f)
        print(f"[i] Resuming: {len(descriptions)} pairs already done.")

    client = anthropic.Anthropic(api_key=api_key)

    total = len(organ_pairs)
    for idx, (o1, o2) in enumerate(organ_pairs):
        key = f"{o1}|{o2}"
        sym = f"{o2}|{o1}"
        if key in descriptions or sym in descriptions:
            print(f"  [{idx+1}/{total}] Cached: {o1} <-> {o2}")
            continue

        edge_text = (
            edge_metadata.get((o1, o2), "")
            or edge_metadata.get((o2, o1), "")
        )

        # Pull key players from literature results
        lit = (
            literature_results.get(key)
            or literature_results.get(sym)
            or {}
        )
        kp = lit.get("key_players", {})
        all_kp = (
            kp.get("hormones", [])
            + kp.get("metabolites", [])
            + kp.get("proteins", [])
        )

        prompt = _build_prompt(o1, o2, edge_text, all_kp)

        print(f"  [{idx+1}/{total}] Generating: {o1} <-> {o2} … ", end="", flush=True)
        try:
            response = client.messages.create(
                model=MODEL,
                max_tokens=300,
                messages=[{"role": "user", "content": prompt}],
            )
            text = response.content[0].text.strip()
            descriptions[key] = {
                "description": text,
                "generated_at": datetime.now().isoformat(),
                "model": MODEL,
            }
            print("ok")
        except Exception as exc:
            print(f"ERROR: {exc}")
            descriptions[key] = {
                "description": "",
                "generated_at": datetime.now().isoformat(),
                "model": MODEL,
                "error": str(exc),
            }

        with open(output_path, "w", encoding="utf-8") as f:
            json.dump(descriptions, f, indent=2, ensure_ascii=False)

    print(f"\n[ok] {len(descriptions)} pairs. Saved to: {output_path}")
    return descriptions


def load_llm_descriptions(path: "str | Path") -> dict:
    path = Path(path)
    if not path.exists():
        return {}
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def get_description(descriptions: dict, organ1: str, organ2: str) -> str:
    """Return the LLM description text for an organ pair, or empty string."""
    entry = (
        descriptions.get(f"{organ1}|{organ2}")
        or descriptions.get(f"{organ2}|{organ1}")
        or {}
    )
    return entry.get("description", "")


def main():
    parser = argparse.ArgumentParser(
        description="Generate LLM descriptions for organ-organ connections."
    )
    parser.add_argument("--reset", action="store_true",
                        help="Delete cache and regenerate all descriptions.")
    args = parser.parse_args()

    from Data_Loader.load_data import load_edge_metadata_from_csv
    from Literature_Search.pubmed_search import load_literature_results

    connection_data = HERE / "metabolic_data" / "connection_data.csv"
    literature_path = HERE / "metabolic_data" / "literature_results.json"

    edge_metadata = load_edge_metadata_from_csv(str(connection_data))
    lit_results = load_literature_results(literature_path)

    pairs = sorted({
        (o1, o2) if o1 < o2 else (o2, o1)
        for (o1, o2) in edge_metadata
    })
    print(f"[i] {len(pairs)} organ-organ pairs to describe.")

    generate_llm_descriptions(
        organ_pairs=pairs,
        edge_metadata=edge_metadata,
        literature_results=lit_results,
        output_path=DEFAULT_OUTPUT,
        resume=not args.reset,
        reset=args.reset,
    )


if __name__ == "__main__":
    main()
