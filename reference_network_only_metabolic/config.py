# =============================================================================
#  REFERENCE NETWORK — ONLY METABOLIC
#  Configuration file
# =============================================================================
#  Edit the values below to change how the network is built.
# =============================================================================

# -----------------------------------------------------------------------------
#  PubMed search parameters
# -----------------------------------------------------------------------------

MAX_PAPERS = 1000   # Maximum papers fetched per organ pair.
                    # Higher = better recall, slower run.

YEARS_BACK = 10     # How many years back to search in PubMed.

DELAY = 0.4         # Seconds between NCBI API requests.
                    # Keep at ≥ 0.34 to avoid rate-limiting.

# -----------------------------------------------------------------------------
#  Metabolic keyword filter
# -----------------------------------------------------------------------------
#  Every retrieved paper MUST contain at least one of these words or phrases
#  in its title or abstract to be included.
#
#  - Single words:         no quotes needed  →  "glucose"
#  - Multi-word phrases:   keep as-is        →  "glucose uptake"
#  - Hyphenated terms:     keep as-is        →  "beta-oxidation"
#
#  Add a term to broaden the search.
#  Remove a term to make the search more specific.
# -----------------------------------------------------------------------------

METABOLIC_KEYWORDS: list[str] = [
    # General metabolic
    "metabolism",
    "metabolic",
    "catabolism",
    "anabolism",
    "biosynthesis",
    "metabolite",
    "flux",
    "substrate",
    "nutrient",
    "fuel",
    # Carbohydrate / glucose
    "glucose",
    "glucose uptake",
    "glycolysis",
    "gluconeogenesis",
    # Energy / mitochondria
    "energy",
    "bioenergetics",
    "mitochondria",
    "mitochondrial",
    "thermogenesis",
    "Krebs cycle",
    "electron transport",
    "oxidation",
]

# -----------------------------------------------------------------------------
#  Crosstalk keyword filter
# -----------------------------------------------------------------------------
#  At least ONE of these must appear in title/abstract in addition to the
#  metabolic filter above.
#
#  Goal: restrict papers to those actually discussing inter-organ
#  communication/relationships, not just papers that happen to mention both
#  organ names and a metabolic term without relating them to each other.
# -----------------------------------------------------------------------------

CROSSTALK_KEYWORDS: list[str] = [
    "network",
    "interplay",
    "cross-talk",
    "crosstalk",
    "connection",
    "axis",
    "inter-organ",
    "interorgan",
    "cross-organ",
    "communication",
    "bidirectional",
    "feedback loop",
    "signaling",
    "regulation",
    "mediates",
    "modulates",
    "interaction",
    "link",
    "coupling",
    "cross-regulation",
    "pathway",
]

# -----------------------------------------------------------------------------
#  LLM connection-type classification
# -----------------------------------------------------------------------------
#  For each organ pair the LLM reads up to LLM_MAX_PAPERS paper abstracts and
#  assigns between 1 and 3 connection types from the list below (only the
#  most-feasible one is required; a second and third are added only if
#  clearly supported). Classification runs 3 independent times per pair and
#  keeps only the types that appear in at least 2 of the 3 runs (majority
#  vote) — see Literature_Search/llm_connection_type.py.
#
#  Add or remove types freely — the LLM prompt is built automatically from
#  this dict.  Keys are short identifiers used internally; values are the
#  human-readable label + description shown in the visualization.
#
#  LLM_MODEL must match an Ollama model you have pulled locally.
#  Run   ollama pull <model-name>   if it is missing.
# -----------------------------------------------------------------------------

CONNECTION_TYPES: dict[str, dict] = {
    "direct_glucose_transfer": {
        "label":       "Direct glucose transfer",
        "description": "One organ releases glucose into circulation that the "
                       "other directly takes up as primary fuel.",
    },
    "hormonal_glucose_control": {
        "label":       "Hormonal glucose control",
        "description": "One organ secretes a hormone (insulin, glucagon, GLP-1, GIP) "
                       "that directly regulates glucose handling in the other.",
    },
    "Other_hormonal_control": {
            "label":       "Other Hormonal influence",
            "description": "One organ secretes a hormone that influences "
                           "the metabolism of function of another organ while "
                           "not directly affecting glucose handling.",
    },
    "nerve_connection": {
        "label":       "Nerve connection",
        "description": "The interaction is via an actual physical nerve connection"
                       " (e.g. vagus nerve in the gut-brain-axis).",
    },
    "other_physical_connection": {
        "label":       "Other physical connection",
        "description": "Non nerve but still physical connection "
                       "Could e.g. be mechanical influence of one organ on the other",
    },
    "ketone_body_exchange": {
        "label":       "Ketone body exchange",
        "description": "Ketone bodies (beta-hydroxybutyrate, acetoacetate) produced "
                       "by hepatic ketogenesis are exported and oxidized as fuel "
                       "by the other organ, especially during fasting or low-carb states.",
    },
    "mineral_electrolyte_regulation": {
        "label":       "Mineral/electrolyte-coupled metabolic regulation",
        "description": "One organ regulates a mineral or electrolyte (calcium, "
                       "phosphate, iron) whose circulating level feeds back on "
                       "metabolic enzyme activity or energy handling in the other.",
    },
    "adipokine_myokine_signaling": {
        "label":       "Adipokine/myokine signaling",
        "description": "One organ secretes an adipokine or myokine (leptin, "
                       "adiponectin, irisin, FGF21) that acts on metabolic pathways "
                       "in the other, beyond core glucose-regulating hormones.",
    },
    "inflammatory_cytokine_crosstalk": {
        "label":       "Inflammatory/cytokine crosstalk",
        "description": "Cytokines released by one organ drive "
                       "inflammation that alters metabolism or "
                       "function in the other.",
    },
    "Metabolite_crosstalk": {
        "label":       "Metabolite crosstalk",
        "description": "Metabolites produced by one organ"
                       "act directly on the other organ to influence its metabolism.",
    },
    "circadian_clock_synchronization": {
        "label":       "Circadian/clock synchronization",
        "description": "The molecular clock in one organ entrains or is entrained "
                       "by the other's clock, and misalignment between them is "
                       "linked to disrupted metabolic rhythms.",
    },
}

LLM_MODEL      = "llama3.2"   # Ollama model name
LLM_MAX_PAPERS = 100            # Max abstracts sent per pair

# -----------------------------------------------------------------------------
#  Visualization
# -----------------------------------------------------------------------------

VIZ_TITLE = "Metabolic Reference Network — Metabolic Query"


# =============================================================================
#  HOW TO RUN THIS NETWORK
#  All commands are run from the Metabolic_Reference_Network/ folder.
#  Requires Ollama running locally with LLM_MODEL pulled (see above).
#
#  This pipeline has four independent steps:
#    1. PubMed search   — fetches papers for each organ-organ pair
#    2. LLM edge types  — classifies each pair into a connection type
#    3. LLM edge summaries — writes a prose summary per pair from abstracts
#    4. LLM organ descriptions — fetches papers per organ, writes a summary
#  Steps 2–4 are all LLM steps and can be run in any order after step 1.
#  The visualization reads all cached outputs and can be rebuilt at any time.
# =============================================================================
#
#  ── FULL RUN — everything from scratch ──────────────────────────────────────
#
#      python run_metabolic_lit_search.py --reset
#      python Literature_Search/llm_descriptions.py --reset
#      python Literature_Search/organ_descriptions.py --reset
#
#      Wipes all caches and reruns search + all LLM steps + visualization.
#      Expect several hours for a full run (search + 18 GB model).
#
#  ── LITERATURE SEARCH ONLY (no LLM, no viz) ─────────────────────────────────
#
#      python run_metabolic_lit_search.py --skip-llm-type --viz-only
#
#      Fetches PubMed papers for all pairs in healthy_cohort_connections.csv.
#      Already-cached pairs are skipped automatically.
#      Add --reset to wipe the search cache and re-fetch everything.
#      Add --force-empty to re-search pairs that previously returned 0 papers.
#
#  ── LLM STEPS ONLY (no search, no viz) ──────────────────────────────────────
#
#      Requires the search cache to exist (run the search step first).
#
#      # Classify each organ-pair into a connection type:
#      python run_metabolic_lit_search.py --skip-llm-type
#      # (omit --skip-llm-type to run classification; cached pairs are skipped)
#      python run_metabolic_lit_search.py
#
#      # Write a prose summary for each organ-pair from its abstracts:
#      python Literature_Search/llm_descriptions.py \
#          --edge-filter reference_network_only_metabolic/healthy_cohort_connections.csv \
#          --literature  reference_network_only_metabolic/metabolic_literature_results.json \
#          --output      reference_network_only_metabolic/metabolic_llm_descriptions.json
#
#      # Fetch papers per organ and write an organ-level metabolic description:
#      python Literature_Search/organ_descriptions.py
#
#      To redo a specific LLM step from scratch, add --reset to that command.
#      To reclassify all pairs after editing CONNECTION_TYPES above:
#      python run_metabolic_lit_search.py --reset-llm-type
#
#  ── VISUALIZATION ONLY (no search, no LLM) ───────────────────────────────────
#
#      python run_metabolic_lit_search.py --viz-only
#
#      Reads all existing JSON caches and regenerates the HTML only.
#      Use after editing VIZ_TITLE, CONNECTION_TYPES labels, or layout.
#      Output: reference_network_only_metabolic/metabolic_literature_network.html
#
# =============================================================================
