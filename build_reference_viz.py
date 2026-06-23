"""
Build the reference network HTML visualization (standalone — no input matrix needed).
Merges existing curated edge metadata with any available PubMed literature results.

Run from any directory:
    python build_reference_viz.py

Output: <this script's folder>/metabolic_data/reference_network.html
"""

import sys
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))

import networkx as nx
from Data_Loader.load_data import load_node_metadata_from_csv, load_edge_metadata_from_csv
from Literature_Search.pubmed_search import load_literature_results, merge_with_edge_metadata
from Visualisation.networkBuilderUtils import export_network_to_cytoscape_dashboard

ORGAN_DATA        = HERE / "metabolic_data" / "organ_data.csv"
CONNECTION_DATA   = HERE / "metabolic_data" / "connection_data.csv"
LITERATURE_RESULTS = HERE / "metabolic_data" / "literature_results.json"
OUTPUT_HTML       = HERE / "metabolic_data" / "reference_network.html"


def main():
    node_metadata = load_node_metadata_from_csv(str(ORGAN_DATA))
    edge_metadata = load_edge_metadata_from_csv(str(CONNECTION_DATA))

    lit_results = load_literature_results(LITERATURE_RESULTS)
    if lit_results:
        n_with_papers = sum(1 for v in lit_results.values() if v.get("n_papers_found", 0) > 0)
        print(f"[i] Literature results: {len(lit_results)} edges, {n_with_papers} with papers.")
    else:
        print("[i] No literature results found — run run_literature_search.py first.")

    merged_edge_data = merge_with_edge_metadata(edge_metadata, lit_results)

    G = nx.Graph()

    for organ, desc in node_metadata.items():
        G.add_node(organ, description=desc)

    for (o1, o2), text in edge_metadata.items():
        if o1 >= o2:
            continue  # add each edge once
        G.add_edge(o1, o2)
        G.edges[o1, o2]['description'] = text
        G.edges[o1, o2]['merged_data'] = merged_edge_data.get((o1, o2), {})
        G.edges[o1, o2]['color'] = "#64748b"

    export_network_to_cytoscape_dashboard(
        graph=G,
        filename=str(OUTPUT_HTML),
        include_legend=False,
        title="Metabolic Reference Network",
    )
    print(f"[ok] Saved to: {OUTPUT_HTML}")


if __name__ == "__main__":
    main()
