"""
Generate LLM-based literature summaries for each organ-organ connection pair.

Uses a locally running Ollama model to write a cited, evidence-based summary
(up to 10 sentences) for each edge, grounded in the PubMed papers already
fetched for that pair.

Requirements
------------
  1. Install Ollama: https://ollama.com
  2. Pull the model:  ollama pull north-mini-code-1.0
  3. Ollama must be running (it starts automatically on most systems)

Results are saved incrementally to metabolic_data/llm_descriptions.json.

Usage (standalone):
    uv run python Literature_Search/llm_descriptions.py
    uv run python Literature_Search/llm_descriptions.py --reset
    uv run python Literature_Search/llm_descriptions.py --model llama3.2
"""

import json
import sys
import argparse
import textwrap
from pathlib import Path
from datetime import datetime

import requests as _requests   # already a project dependency

HERE = Path(__file__).parent.parent
sys.path.insert(0, str(HERE))

DEFAULT_OUTPUT     = HERE / "metabolic_data" / "llm_descriptions.json"
DEFAULT_MODEL      = "north-mini-code-1.0"
OLLAMA_URL         = "http://localhost:11434/api/chat"
MAX_PAPERS_PROMPT  = 8     # how many papers to include in the prompt
MAX_ABSTRACT_CHARS = 600   # truncate long abstracts to keep prompt lean
OLLAMA_TIMEOUT     = 120   # seconds — local inference can be slow


# ---------------------------------------------------------------------------
# Prompt builder
# ---------------------------------------------------------------------------

def _build_prompt(organ1: str, organ2: str, papers: list[dict]) -> str:
    """
    Build a prompt that includes the actual paper titles and abstracts so the
    model can write a grounded, cited summary.
    """
    if not papers:
        paper_block = "No papers available for this pair."
    else:
        lines = []
        for i, p in enumerate(papers[:MAX_PAPERS_PROMPT], 1):
            pmid    = p.get("pmid", "?")
            title   = p.get("title", "No title").strip()
            year    = p.get("year", "")
            abstract = p.get("abstract", "").strip()
            if len(abstract) > MAX_ABSTRACT_CHARS:
                abstract = abstract[:MAX_ABSTRACT_CHARS].rsplit(" ", 1)[0] + "…"
            lines.append(
                f"[{i}] PMID {pmid} ({year})\n"
                f"    Title: {title}\n"
                f"    Abstract: {abstract}"
            )
        paper_block = "\n\n".join(lines)

    return textwrap.dedent(f"""
        You are an expert in metabolic physiology and inter-organ communication.

        Task: Write a scientific summary (maximum 10 sentences) explaining the
        metabolic and hormonal basis of the {organ1}–{organ2} axis based solely
        on the research papers listed below. Cite papers using [PMID XXXXX]
        inline. Focus on: key signalling molecules, metabolic substrates,
        hormones, physiological mechanisms, and any disease relevance.
        Write flowing prose — no bullet points, no headings.

        Papers:
        {paper_block}

        Summary:
    """).strip()


# ---------------------------------------------------------------------------
# Ollama interface
# ---------------------------------------------------------------------------

def _check_ollama(model: str) -> bool:
    """Return True if Ollama is reachable and the model is available."""
    try:
        resp = _requests.get("http://localhost:11434/api/tags", timeout=5)
        if resp.status_code != 200:
            return False
        available = [m["name"] for m in resp.json().get("models", [])]
        # Allow partial match (e.g. "llama3.2" matches "llama3.2:latest")
        return any(model in m for m in available)
    except Exception:
        return False


def _call_ollama(prompt: str, model: str) -> str:
    """Send a chat request to Ollama and return the response text."""
    payload = {
        "model":    model,
        "messages": [{"role": "user", "content": prompt}],
        "stream":   False,
        "options":  {"temperature": 0.3, "num_predict": 600},
    }
    resp = _requests.post(OLLAMA_URL, json=payload, timeout=OLLAMA_TIMEOUT)
    resp.raise_for_status()
    return resp.json()["message"]["content"].strip()


# ---------------------------------------------------------------------------
# Main generation function
# ---------------------------------------------------------------------------

def generate_llm_descriptions(
    organ_pairs: list[tuple[str, str]],
    literature_results: dict,
    output_path: "str | Path" = DEFAULT_OUTPUT,
    model: str = DEFAULT_MODEL,
    resume: bool = True,
    reset: bool = False,
) -> dict:
    """
    Generate literature-grounded LLM summaries for all organ pairs.

    Parameters
    ----------
    organ_pairs        : sorted list of (organ1, organ2) tuples
    literature_results : raw literature JSON (keyed "organ1|organ2")
    output_path        : where to save/load llm_descriptions.json
    model              : Ollama model name
    resume             : skip pairs already in the cache file
    reset              : delete cache and regenerate everything
    """
    output_path = Path(output_path)

    # Load existing cache
    descriptions: dict = {}
    if reset and output_path.exists():
        output_path.unlink()
        print(f"[i] Cache deleted: {output_path.name}")
    elif resume and output_path.exists():
        with open(output_path, encoding="utf-8") as f:
            descriptions = json.load(f)
        # Count only entries that were generated by this model and have content
        usable = sum(
            1 for v in descriptions.values()
            if v.get("description") and v.get("model") == model
        )
        stale = len(descriptions) - usable
        print(f"[i] Cache: {usable} usable entries for model '{model}'"
              + (f", {stale} stale/wrong-model (will regenerate)" if stale else "") + ".")

    # Verify Ollama is running and model is available
    if not _check_ollama(model):
        print(f"[!] Ollama not reachable or model '{model}' not found.")
        print(f"    Start Ollama and run:  ollama pull {model}")
        print(f"    Then re-run with:      uv run python run_all.py --skip-search")
        return descriptions

    print(f"[i] Using model: {model}")

    total = len(organ_pairs)
    for idx, (o1, o2) in enumerate(organ_pairs):
        key = f"{o1}|{o2}"
        sym = f"{o2}|{o1}"

        existing = descriptions.get(key) or descriptions.get(sym)
        if existing and existing.get("description") and existing.get("model") == model:
            print(f"  [{idx+1}/{total}] Cached: {o1} <-> {o2}")
            continue

        # Retrieve papers for this pair
        lit = (
            literature_results.get(key)
            or literature_results.get(sym)
            or {}
        )
        papers = lit.get("papers", [])
        n_papers = lit.get("n_papers_found", len(papers))

        print(f"  [{idx+1}/{total}] {o1} <-> {o2}  ({n_papers} papers) … ",
              end="", flush=True)

        if n_papers == 0:
            print("skipped (no papers)")
            descriptions[key] = {
                "description":  "",
                "generated_at": datetime.now().isoformat(),
                "model":        model,
                "n_papers_used": 0,
            }
        else:
            prompt = _build_prompt(o1, o2, papers)
            try:
                text = _call_ollama(prompt, model)
                descriptions[key] = {
                    "description":   text,
                    "generated_at":  datetime.now().isoformat(),
                    "model":         model,
                    "n_papers_used": min(len(papers), MAX_PAPERS_PROMPT),
                }
                print(f"ok ({len(text.split())} words)")
            except Exception as exc:
                print(f"ERROR: {exc}")
                descriptions[key] = {
                    "description":   "",
                    "generated_at":  datetime.now().isoformat(),
                    "model":         model,
                    "n_papers_used": 0,
                    "error":         str(exc),
                }

        # Save after every pair so progress is never lost
        with open(output_path, "w", encoding="utf-8") as f:
            json.dump(descriptions, f, indent=2, ensure_ascii=False)

    done = sum(1 for v in descriptions.values() if v.get("description"))
    print(f"\n[ok] {done}/{len(descriptions)} pairs have summaries. Saved to: {output_path}")
    return descriptions


# ---------------------------------------------------------------------------
# Helpers used by other modules
# ---------------------------------------------------------------------------

def load_llm_descriptions(path: "str | Path") -> dict:
    path = Path(path)
    if not path.exists():
        return {}
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def get_description(descriptions: dict, organ1: str, organ2: str) -> str:
    """Return the LLM summary text for an organ pair, or empty string."""
    entry = (
        descriptions.get(f"{organ1}|{organ2}")
        or descriptions.get(f"{organ2}|{organ1}")
        or {}
    )
    return entry.get("description", "")


# ---------------------------------------------------------------------------
# Standalone entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Generate literature-grounded LLM summaries for organ-organ connections."
    )
    parser.add_argument("--reset", action="store_true",
                        help="Delete cache and regenerate all descriptions.")
    parser.add_argument("--model", default=DEFAULT_MODEL,
                        help=f"Ollama model to use (default: {DEFAULT_MODEL}).")
    args = parser.parse_args()

    from Literature_Search.pubmed_search import load_literature_results
    from Data_Loader.load_data import load_edge_metadata_from_csv

    connection_data = HERE / "metabolic_data" / "connection_data.csv"
    literature_path = HERE / "metabolic_data" / "literature_results.json"

    edge_metadata = load_edge_metadata_from_csv(str(connection_data))
    lit_results   = load_literature_results(literature_path)

    pairs = sorted({
        (o1, o2) if o1 < o2 else (o2, o1)
        for (o1, o2) in edge_metadata
    })
    print(f"[i] {len(pairs)} organ-organ pairs to summarise.")

    generate_llm_descriptions(
        organ_pairs=pairs,
        literature_results=lit_results,
        output_path=DEFAULT_OUTPUT,
        model=args.model,
        resume=not args.reset,
        reset=args.reset,
    )


if __name__ == "__main__":
    main()
