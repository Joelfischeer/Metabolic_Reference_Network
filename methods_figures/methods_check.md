# Methods text vs. implementation (checked 2026-10-09)

Status: **OK** = matches the code; **CHANGED** = text is outdated or wrong; **UNVERIFIED** = cannot be checked from the repository.

## Introduction paragraph

| Sentence | Status | Finding |
|---|---|---|
| "The method entails a structured literature search for each organ-organ pair ..." | OK | 14 organs, 91 pairs (`all_organ_connections.csv`). |
| "In the visualization the organ-organ connections can be selected based on healthy or obese subjects." | OK | Healthy/Obese switch exists. |
| "These connections are selected based on the PET based partial correlation network ..." | **CHANGED** | The switch no longer uses the PET edges. A connection is shown for a condition if it has at least one paper matching the condition keywords (70 of 91 pairs healthy, 72 obese; `build_viz`, `lit_healthy`/`lit_obese`). The PET cohort edges (16 healthy, 9 obese) are only used for the comparison buttons and for the condition-specific LLM files. |
| "healthy subjects (n = )" | to fill | PET file: 241 records, 197 with BMI; 112 with BMI <= 25 and 85 with BMI > 25 (`compute_cohort_PC_networks.py`). |
| "BMI under 25 ... above 25" | **CHANGED** | Code uses BMI <= 25 for healthy and BMI > 25 for "obese". Only 53 of the 85 have BMI >= 30, so "obese" includes overweight. |
| PET networks from Geist et al. [9] | **UNVERIFIED** | The repository computes its own partial-correlation networks (`compute_cohort_PC_networks.py`, p < 0.05 uncorrected, 13 organs, Colon absent). Edge counts are 16 (healthy) and 9 (obese), not 35 and 11. Please confirm whether these are the Geist et al. networks. |

## Step 1: Paper retrieval

| Sentence | Status | Finding |
|---|---|---|
| MeSH term plus text aliases in Title/Abstract | OK | `_alias_clause`: `(MeSH OR "alias"[Title/Abstract] ...)`. |
| Metabolic clause 22 terms; crosstalk filter 21 terms; query structure | OK | `config.py` lists exactly 22 and 21 terms; `build_query`. |
| Local filters rule-based, not LLM | OK | `_filter_same_sentence_crosstalk`. |
| Same-sentence co-location, split on . ! ? ; | OK | `re.split(r"(?<=[.!?;])\s+", ...)`. |
| Whole-word and negation-aware | OK | Negative look-behind for "non-", "non ", "non", "not ". |
| Multi-organ disambiguation (list groups, other organ before the later of the two) | OK | `_list_groups`, `_passes_multi_organ_check` implement exactly the described rule. |
| "Maximum 1,000 papers per pair" | **CHANGED** | `MAX_PAPERS = 20000`. |
| "10-year look-back window" | **CHANGED** | `YEARS_BACK = 20`. |
| "0.4 s delay" | OK | `DELAY = 0.4`. |
| "PMID, title, abstract, year and author keywords" | **CHANGED** | The stored `keywords` are MeSH descriptors (`<DescriptorName>`), not author keywords. |
| (not in text) | add | Result sizes: 187,757 hits retrieved (summed over pairs); 27,416 pair-paper assignments kept; 22,451 unique papers; 87 of 91 pairs with at least one paper; 36 pairs with more than 100 papers; maximum 4,986 (Brain-Small Intestine). Search date 2026-09-17. |

## Step 2: Subselection based on condition

| Sentence | Status | Finding |
|---|---|---|
| Edge-level restriction to 35 healthy / 11 obese pairs | **CHANGED** | Not applied to the dashboard any more. Condition views cover all 91 pairs with at least one condition paper (70 healthy, 72 obese). The cohort edge sets (16 / 9) only define the pairs that get a separate condition-specific LLM classification and description. |
| Paper-level restriction: at least one keyword anywhere in title/abstract | OK | `_filter_by_condition`. |
| Healthy 7 terms, obese 11 terms | OK | `config.py` (obese list includes four "overweight" terms). |
| Not mutually exclusive | OK | 200 papers match both; unique papers: healthy 1,585, obese 3,594. |
| (not in text) | add | Sums over pairs (pair-paper assignments): healthy 1,973, obese 4,373. The flowchart currently uses these assignment counts for A/B/C (27,416 / 1,973 / 4,373); the unique counts are 22,451 / 1,585 / 3,594. |

## Step 3: Key player extraction

| Sentence | Status | Finding |
|---|---|---|
| Titles, abstracts, MeSH descriptors scanned; number of papers per term counted | OK | `extract_key_players`; each term once per paper. |
| Hormone vocabulary 170 terms | **CHANGED** | Now 171 (alias "glucagon-like peptide-1" added). |
| Metabolites 267 -> 254 terms | OK / check | Final 254 matches. The "267" and "8 removed" do not add up (267 - 8 = 259). |
| Proteins 45,236 -> 45,138 | OK / check | Final 45,138 matches the current vocabulary. 45,236 - 84 = 45,152, so the intermediate numbers do not reconcile. |
| Hormones removed from the other lists | OK | Vocabulary overlaps are 0. |
| Top-10,000 English-word filter (first20hours list) | **UNVERIFIED** | No script or reference in the repository. |
| Exceptions kept: ACC, ACE, CLOCK, INS, glucose, oxygen, cholesterol, AMP, GDP, dUMP | **CHANGED** | CLOCK was removed (alias matched the ordinary word "clock"). All others are in the vocabulary. |
| HPA (protein-level evidence) and Stanford pathway map as sources | **UNVERIFIED** | Not documented in the code. |
| Hormone and metabolite synonyms created with ChatGPT | **UNVERIFIED** | Hormone groups are hand-written; metabolite dictionary has only 6 groups. |
| HPA synonym rule (name, symbol, hyphenated variants) | OK | e.g. `il2`, `il-2`, `interleukin 2`, `interleukin-2`; no historical aliases. |
| Grouping applied before counting | OK | Merging happens in `_find_canonical_terms`. |
| Curation paragraphs, hyphen rule | OK | Implemented in the curation block of `pubmed_search.py`. |
| (not in text) | add | Network-level counts sum per-pair counts, so a paper matching several pairs counts once for each. |

## Step 4: LLM connection-type classification

| Sentence | Status | Finding |
|---|---|---|
| "100 papers are sampled if more than 100" | **CHANGED** | The first 100 papers in stored order are used (`papers[:max_papers]`), not a random sample. |
| Abstracts truncated to 400 characters | OK | `MAX_ABSTRACT_CHARS = 400` (cut at a word boundary). |
| llama3.2 via Ollama | OK | `LLM_MODEL`. |
| "rank the connection types ... at least one, 2 or 3 only if well supported" | OK (wording) | The model returns TYPE_1 (best) and optionally TYPE_2/3. |
| "run 3 times; types in at least two of the three runs" | **CHANGED** | `VOTE_RUNS = 5`, majority = more than half (at least 3 of 5), temperature 0.5. The comments in `config.py` and the module docstring still say 3 runs. |
| Fallback: most voted type if all disagree | OK | `_majority_vote`. |
| 11 connection types | OK | "Hormonal influence" is labelled "Other Hormonal influence" in the code. |

## Step 5: LLM edge descriptions

| Sentence | Status | Finding |
|---|---|---|
| 25 random papers (all if fewer) | OK | `SAMPLE_POOL = 25`, `random.sample`. |
| LLM picks the 5 most relevant | OK | `TOP_N_PAPERS = 5`; skipped if 5 or fewer; padded if the LLM returns fewer than 5 valid indices. |
| 3 to 4 sentence cited summary | OK | The prompt requires every selected paper to be cited. The module docstring still says "up to 10 sentences". |

## Step 6: LLM organ descriptions

| Sentence | Status | Finding |
|---|---|---|
| Up to 50 papers, last 5 years, same metabolic keywords | OK | `MAX_RESULTS = 50`, `YEARS_BACK = 5`; 14 x 50 = 700 papers. |
| "six are selected" | clarify | The first six retrieved papers are used (`papers[:6]`), not an LLM selection. |
| 400 characters, exactly 5 sentences | OK | All 14 descriptions have 5 sentences. |
| Refusal detection with fallback prompt | OK | `_is_refusal`. |

## Step 7: Visualization

| Sentence | Status | Finding |
|---|---|---|
| "four JSON caches from steps 1 to 4" | **CHANGED** | Caches used: results (steps 1 and 3), connection types (4), connection descriptions (5), organ descriptions (6), plus the condition-specific files. |
| NetworkX undirected graph | OK | `nx.Graph()`. |
| "edges are the 35 organ-organ pairs" | **CHANGED** | 14 nodes and 91 edges. |
| "self-contained single HTML file ... Cytoscape.js" | partly | One HTML file with embedded data, but Cytoscape.js 3.29.2 is loaded from unpkg, so an internet connection is needed to open it. |
| (not in text) | add | OO-based edge weights with a per-condition elbow threshold; Literature and Edge Statistics tabs; comparison with the PET cohort networks. |

## Flowchart (Figure 4)

* "Pubmed database last 20 years": matches the code (the text says 10 years).
* n = 27,416 / 1,973 / 4,373 are sums over pairs. Unique papers are 22,451 / 1,585 / 3,594.
* "only healthy" / "only obese": the subsets are not exclusive (200 papers in both). Use "healthy keywords" / "obese keywords".
* "Extracted papers based on organs n = 700": matches (14 organs x 50).
