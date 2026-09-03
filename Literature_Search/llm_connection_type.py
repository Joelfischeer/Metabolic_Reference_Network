"""
LLM-based connection-type classification for organ-organ pairs.

For each pair the LLM reads up to MAX_PAPERS paper abstracts and picks
between 1 and 3 connection types (ordered by relevance) from a user-defined
list. To make the classification robust against single-sample noise, each
pair is classified 3 times independently and only the types that appear in
at least 2 of the 3 runs (majority vote) are kept; if no type reaches a
majority, the single most-voted type is kept as a fallback so every pair
always ends up with at least one type.

Results are saved incrementally to a JSON file so the run can be resumed.
"""

import json
import textwrap
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import requests as _requests

OLLAMA_URL       = "http://localhost:11434/api/chat"
OLLAMA_TIMEOUT   = 180
OLLAMA_RETRIES   = 2
MAX_ABSTRACT_CHARS = 400
VOTE_RUNS        = 3     # independent classification passes per pair
VOTE_TEMPERATURE = 0.5   # higher than a single-shot call, so the 3 runs can
                          # actually disagree — a majority vote over 3 near-
                          # identical low-temperature replies is meaningless
MAX_TYPES        = 3


# ── Ollama helpers ─────────────────────────────────────────────────────────────

def _check_ollama(model: str) -> bool:
    try:
        resp = _requests.get("http://localhost:11434/api/tags", timeout=5)
        if resp.status_code != 200:
            return False
        available = [m["name"] for m in resp.json().get("models", [])]
        return any(model in m for m in available)
    except Exception:
        return False


def _call_ollama(prompt: str, model: str, temperature: float = 0.1) -> str:
    payload = {
        "model":    model,
        "messages": [{"role": "user", "content": prompt}],
        "stream":   False,
        "options":  {"temperature": temperature, "num_predict": 512},
    }
    for attempt in range(1, OLLAMA_RETRIES + 1):
        try:
            resp = _requests.post(OLLAMA_URL, json=payload, timeout=OLLAMA_TIMEOUT)
            resp.raise_for_status()
            msg = resp.json()["message"]
            # Thinking models (e.g. north-mini-code) put output in 'thinking';
            # standard models use 'content'.
            content = msg.get("content", "").strip()
            if not content:
                content = msg.get("thinking", "").strip()
            return content
        except _requests.exceptions.Timeout:
            if attempt < OLLAMA_RETRIES:
                print(f"\n    [!] Timeout (attempt {attempt}/{OLLAMA_RETRIES}), retrying…",
                      end=" ", flush=True)
            else:
                raise


# ── Prompt builder ─────────────────────────────────────────────────────────────

def _build_prompt(
    organ1: str,
    organ2: str,
    papers: list[dict],
    connection_types: dict[str, dict],
    max_papers: int,
) -> str:
    # Build numbered paper block
    if not papers:
        paper_block = "No papers available."
    else:
        lines = []
        for i, p in enumerate(papers[:max_papers], 1):
            title    = p.get("title", "").strip()
            abstract = p.get("abstract", "").strip()
            if len(abstract) > MAX_ABSTRACT_CHARS:
                abstract = abstract[:MAX_ABSTRACT_CHARS].rsplit(" ", 1)[0] + "…"
            lines.append(f"[{i}] {title}\n    {abstract}")
        paper_block = "\n\n".join(lines)

    # Build numbered type list
    type_lines = "\n".join(
        f"  {i+1}. {key} — {info['label']}: {info['description']}"
        for i, (key, info) in enumerate(connection_types.items())
    )
    keys = list(connection_types.keys())
    valid_keys = ", ".join(keys)

    return textwrap.dedent(f"""
        You are an expert in metabolic physiology and inter-organ signalling.

        Task: Based on the papers below about the {organ1}–{organ2} connection,
        choose which connection type(s) from the list below best fit the
        evidence. Pick the SINGLE most-feasible type as TYPE_1. Only add
        TYPE_2 and/or TYPE_3 if they are ALSO clearly and independently
        supported by the papers — do not pad the list to reach 3. It is
        fine and common to list only 1 type.

        Connection types:
        {type_lines}

        Papers:
        {paper_block}

        Instructions:
        - Reply with one to three lines and nothing else.
        - Line 1 (required): TYPE_1: <key>
        - Line 2 (optional):  TYPE_2: <key>
        - Line 3 (optional):  TYPE_3: <key>
        - <key> must be one of: {valid_keys}
        - Order matters: TYPE_1 must be the best-matching type.
        - Omit TYPE_2/TYPE_3 entirely (do not write the line) if no further
          type clearly applies.

        Reply:
    """).strip()


# ── Parser ─────────────────────────────────────────────────────────────────────

def _parse_response(text: str, valid_keys: list[str]) -> list[str]:
    """Extract an ordered, deduplicated list of up to MAX_TYPES keys."""
    types: list[str] = []
    for line in text.splitlines():
        line = line.strip()
        if not line.upper().startswith("TYPE_"):
            continue
        if ":" not in line:
            continue
        candidate = line.split(":", 1)[1].strip().lower().replace(" ", "_")
        if candidate in valid_keys and candidate not in types:
            types.append(candidate)

    # Fallback: scan for any valid key mentioned anywhere in the reply,
    # in the order they appear, in case the model didn't follow the format.
    if not types:
        lowered = text.lower()
        found = [(lowered.index(key), key) for key in valid_keys if key in lowered]
        types = [key for _, key in sorted(found)]

    return types[:MAX_TYPES]


# ── Majority vote across runs ────────────────────────────────────────────────

def _majority_vote(runs: list[list[str]], max_types: int = MAX_TYPES) -> list[str]:
    """
    Combine N independent ordered type-lists into one final list.

    A type is kept if it appears in a majority of the runs (more than half).
    If no type reaches a majority (e.g. all runs disagree completely), the
    single most-voted type is kept instead, so every pair with at least one
    valid run always ends up with >=1 type. Ties are broken by earlier
    average rank (types listed first / more consistently are more relevant).
    Result is capped at `max_types`.
    """
    votes: dict[str, int] = defaultdict(int)
    rank_sum: dict[str, int] = defaultdict(int)
    order: list[str] = []
    for run in runs:
        for rank, key in enumerate(run):
            votes[key] += 1
            rank_sum[key] += rank
            if key not in order:
                order.append(key)

    if not votes:
        return []

    def sort_key(k: str):
        return (-votes[k], rank_sum[k] / votes[k])

    majority_needed = len(runs) // 2 + 1
    majority = sorted((k for k in order if votes[k] >= majority_needed), key=sort_key)

    if not majority:
        best = min(order, key=sort_key)
        return [best]

    return majority[:max_types]


def _classify_once(
    organ1: str, organ2: str, papers: list[dict],
    connection_types: dict[str, dict], max_papers: int,
    model: str, valid_keys: list[str],
) -> tuple[list[str], str]:
    """Run one independent classification pass. Returns (types, raw_reply)."""
    prompt = _build_prompt(organ1, organ2, papers, connection_types, max_papers)
    raw = _call_ollama(prompt, model, temperature=VOTE_TEMPERATURE)
    return _parse_response(raw, valid_keys), raw


# ── Main function ──────────────────────────────────────────────────────────────

def generate_connection_type_classifications(
    organ_pairs: list[tuple[str, str]],
    literature_results: dict,
    connection_types: dict[str, dict],
    output_path: "str | Path",
    model: str = "north-mini-code-1.0",
    max_papers: int = 5,
    resume: bool = True,
    reset: bool = False,
) -> dict:
    """
    Classify each organ pair into 1-3 connection types using an Ollama LLM.

    Each pair is classified VOTE_RUNS (3) times independently; only types
    appearing in a majority of the runs are kept (falls back to the single
    most-voted type if no type reaches a majority), guaranteeing every
    classified pair gets at least one type and never more than MAX_TYPES.

    Parameters
    ----------
    organ_pairs         : sorted list of (organ1, organ2) tuples
    literature_results  : raw search JSON (keyed "organ1|organ2")
    connection_types    : dict of {key: {label, description}} from config
    output_path         : where to save/load the classification JSON
    model               : Ollama model name
    max_papers          : how many abstracts to include in the prompt
    resume              : skip pairs already in the cache
    reset               : delete cache and reclassify everything
    """
    output_path = Path(output_path)
    cache: dict = {}

    if reset and output_path.exists():
        output_path.unlink()
        print(f"[i] LLM type cache deleted: {output_path.name}")
    elif resume and output_path.exists():
        with open(output_path, encoding="utf-8") as f:
            cache = json.load(f)
        print(f"[i] Resuming LLM type classification: {len(cache)} pairs cached.")

    if not _check_ollama(model):
        print(f"[!] Ollama not reachable or model '{model}' not found.")
        print(f"    Start Ollama and run:  ollama pull {model}")
        return cache

    print(f"[i] LLM type classification using model: {model}")
    valid_keys = list(connection_types.keys())
    total = len(organ_pairs)

    for idx, (o1, o2) in enumerate(organ_pairs):
        key = f"{o1}|{o2}"
        sym = f"{o2}|{o1}"

        if (key in cache or sym in cache) and resume and not reset:
            entry = cache.get(key) or cache.get(sym, {})
            shown = "+".join(entry.get("types", [])) or entry.get("primary", "?")
            print(f"  [{idx+1}/{total}] Cached ({shown}): {o1} <-> {o2}")
            continue

        lit    = literature_results.get(key) or literature_results.get(sym) or {}
        papers = lit.get("papers", [])
        n      = lit.get("n_papers_found", len(papers))

        print(f"  [{idx+1}/{total}] {o1} <-> {o2}  ({n} papers) … ",
              end="", flush=True)

        if not papers:
            print("skipped (no papers)")
            cache[key] = {
                "types": [],
                "generated_at": datetime.now().isoformat(),
                "model": model,
            }
        else:
            runs: list[list[str]] = []
            raw_responses: list[str] = []
            errors: list[str] = []
            for _ in range(VOTE_RUNS):
                try:
                    types, raw = _classify_once(
                        o1, o2, papers, connection_types, max_papers, model, valid_keys
                    )
                    runs.append(types)
                    raw_responses.append(raw)
                except Exception as exc:
                    runs.append([])
                    raw_responses.append("")
                    errors.append(str(exc))

            final_types = _majority_vote(runs)
            cache[key] = {
                "types":         final_types,
                "runs":          runs,
                "raw_responses": raw_responses,
                "generated_at":  datetime.now().isoformat(),
                "model":         model,
            }
            if errors:
                cache[key]["errors"] = errors
            if final_types:
                print(f"types={'+'.join(final_types)}  (votes: {[r for r in runs]})")
            else:
                print(f"ERROR: all {VOTE_RUNS} runs failed" if len(errors) == VOTE_RUNS
                      else "no type matched")

        with open(output_path, "w", encoding="utf-8") as f:
            json.dump(cache, f, indent=2, ensure_ascii=False)

    classified = sum(1 for v in cache.values() if v.get("types"))
    print(f"\n[ok] {classified}/{len(cache)} pairs classified. Saved: {output_path}")
    return cache


# ── Loader helper ──────────────────────────────────────────────────────────────

def load_connection_type_classifications(path: "str | Path") -> dict:
    path = Path(path)
    if not path.exists():
        return {}
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def get_connection_types(
    classifications: dict,
    organ1: str,
    organ2: str,
    connection_types: dict[str, dict],
) -> list[str]:
    """
    Return an ordered list of 0-3 human-readable type labels for an organ
    pair, most-feasible first. Empty list if not classified.

    Also reads old-format caches (single "primary"/"secondary" keys, from
    before the multi-type majority-vote scheme) so a stale cache still
    renders sensibly until it's reclassified with --reset-llm-type.
    """
    entry = (
        classifications.get(f"{organ1}|{organ2}")
        or classifications.get(f"{organ2}|{organ1}")
        or {}
    )
    keys = entry.get("types")
    if keys is None:
        keys = [k for k in (entry.get("primary"), entry.get("secondary")) if k]

    return [connection_types.get(k, {}).get("label", k) for k in keys if k]
