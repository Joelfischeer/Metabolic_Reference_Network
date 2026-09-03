# =============================================================================
#  EDGE COSINE GENERAL REFERENCE NETWORK  (OBESE / HIGH-BMI CONDITION)
#  Configuration file
#
#  Per-organ query covers BOTH the metabolic and hormonal layer (like
#  Edge_general_reference_network):
#      (MeSH_ORGAN OR aliases_ORGAN)
#      AND [METABOLIC_FILTER | HORMONAL_FILTER]
#      AND CONDITION_FILTER
#  Edges are weighted by the Otsuka–Ochiai (cosine) coefficient instead of a
#  raw cross-mention count, so organs with a much larger literature pool
#  don't dominate edges just by volume. See run_network.py's module
#  docstring for details.
# =============================================================================
#  Edit the values below, then re-run:
#      uv run -m Edge_cosine_general_reference_network.run_network --condition obese
#
#  To rebuild visualizations only (no re-search, no re-bootstrap):
#      uv run -m Edge_cosine_general_reference_network.run_network --condition obese --viz-only
#
#  To wipe cache and restart from scratch:
#      uv run -m Edge_cosine_general_reference_network.run_network --condition obese --reset
# =============================================================================

CONDITION_NAME  = "obese"
VIZ_LABEL       = "Obese / High BMI"

# -----------------------------------------------------------------------------
#  PubMed search parameters
# -----------------------------------------------------------------------------

MAX_PAPERS  = 50000  # Maximum papers fetched per organ.
YEARS_BACK  = 10     # How many years back to search in PubMed.
DELAY       = 0.4    # Seconds between NCBI API requests (keep ≥ 0.34).

# -----------------------------------------------------------------------------
#  Edge threshold (applied to full search results before bootstrapping)
# -----------------------------------------------------------------------------

MIN_COOCCUR = 3      # Minimum co-occurring papers to include a pair in the
                     # bootstrap analysis and overview figure.

# -----------------------------------------------------------------------------
#  Metabolic keyword filter  (Layer: metabolic)
# -----------------------------------------------------------------------------
#  At least ONE of these must appear in title/abstract.
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
#  Hormonal keyword filter  (Layer: hormonal)
# -----------------------------------------------------------------------------
#  At least ONE of these must appear in title/abstract.
# -----------------------------------------------------------------------------

HORMONAL_KEYWORDS: list[str] = [
    # General endocrine
    "hormone",
    "hormonal",
    "endocrine",
    "endocrinology",
    "signaling",
    "signal transduction",
    # Pancreatic / glucose-regulating
    "insulin",
    "glucagon",
    "C-peptide",
    "amylin",
    "somatostatin",
    # Gut hormones
    "GLP-1",
    "glucagon-like peptide",
    "GIP",
    "glucose-dependent insulinotropic",
    "peptide YY",
    "PYY",
    "ghrelin",
    "cholecystokinin",
    "CCK",
    "secretin",
    "motilin",
    # Adipokines
    "leptin",
    "adiponectin",
    "resistin",
    "visfatin",
    "adipsin",
    "adipokine",
    # Thyroid
    "thyroid hormone",
    "thyroxine",
    "triiodothyronine",
    "T3",
    "T4",
    "TSH",
    "thyroid-stimulating hormone",
    # Adrenal / stress
    "cortisol",
    "glucocorticoid",
    "mineralocorticoid",
    "aldosterone",
    "adrenaline",
    "epinephrine",
    "norepinephrine",
    "noradrenaline",
    "catecholamine",
    "ACTH",
    "corticotropin",
    # Growth / anabolic
    "growth hormone",
    "IGF-1",
    "insulin-like growth factor",
    # Reproductive
    "testosterone",
    "estrogen",
    "estradiol",
    "progesterone",
    "DHEA",
    "androgen",
    # Calcium / bone
    "parathyroid hormone",
    "PTH",
    "calcitonin",
    "vitamin D",
    "calcitriol",
    # Myokines / hepatokines / organokines
    "irisin",
    "myokine",
    "hepatokine",
    "FGF21",
    "FGF19",
    "fetuin",
    # Hypothalamic / pituitary
    "prolactin",
    "oxytocin",
    "vasopressin",
    "ADH",
    "antidiuretic hormone",
    "GnRH",
    "CRH",
    "TRH",
]

# -----------------------------------------------------------------------------
#  Condition-specific keyword filter  (OBESE / HIGH BMI)
# -----------------------------------------------------------------------------
#  At least ONE of these must appear in title/abstract in addition to the
#  metabolic OR hormonal filter above.
#
#  Goal: restrict papers to those studying obesity, high BMI, or related
#  metabolic disease states.
# -----------------------------------------------------------------------------

CONDITION_KEYWORDS: list[str] = [
    "obesity",
    "obese subjects",
    "obese patients",
    "obese volunteers",
    "obese adults",
    "adiposity",
    "patients with adiposity",
]

# -----------------------------------------------------------------------------
#  Crosstalk keyword filter
# -----------------------------------------------------------------------------
#  At least ONE of these must appear in title/abstract in addition to the
#  metabolic/hormonal and condition filters above.
#
#  Goal: restrict papers to those actually discussing inter-organ
#  communication/relationships, not just papers that happen to mention an
#  organ name and a metabolic/hormonal term without relating them to each other.
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
#  Bootstrap parameters
# -----------------------------------------------------------------------------

N_BOOTSTRAP      = 50     # Number of bootstrap iterations.
SAMPLE_FRACTION  = 0.10   # Fraction of each organ's papers sampled per iteration (10%).

# -----------------------------------------------------------------------------
#  Robust network threshold
# -----------------------------------------------------------------------------
#  Organ pairs with mean Otsuka–Ochiai coefficient ≥ this value appear in the
#  robust network visualization.
# -----------------------------------------------------------------------------

MIN_BOOTSTRAP_MEAN = Elbow   # mean Otsuka–Ochiai coefficient threshold.
                              # `Elbow` auto-selects the kneedle-elbow value from
                              # this condition's own bootstrap distribution each
                              # run (see bootstrap_overview_obese.html).
                              # Set a number instead (e.g. 0.02) to pin it.

# -----------------------------------------------------------------------------
#  Visualization titles
# -----------------------------------------------------------------------------

VIZ_TITLE_OVERVIEW = "Obese / High-BMI Cosine General Network — Bootstrap Strength Overview"
VIZ_TITLE_ROBUST   = "Obese / High-BMI Cosine General Network — Robust Connections"

# -----------------------------------------------------------------------------
#  LLM connection-type classification
# -----------------------------------------------------------------------------
#  For each robust organ pair the LLM reads up to LLM_MAX_PAPERS of this
#  pair's same-sentence cross-mention papers (bootstrap "papers" field) and
#  assigns between 1 and 3 connection types from the list below (only the
#  most-feasible one is required; a second and third are added only if
#  clearly supported). Classification runs 3 independent times per pair and
#  keeps only the types that appear in at least 2 of the 3 runs (majority
#  vote) — see Literature_Search/llm_connection_type.py and
#  Edge_cosine_general_reference_network/run_llm_descriptions.py.
#
#  Kept identical to reference_network_only_metabolic/config.py's
#  CONNECTION_TYPES so the same label means the same thing across all five
#  dashboards. Edit both together if you want to add/remove a type.
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
        "label":       "Hormonal influence",
        "description": "One organ sectretes a hormone that influences "
                       "the metabolism of funciton of another organ while "
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

LLM_MAX_PAPERS = 100   # Max same-sentence cross-mention papers sent per pair
