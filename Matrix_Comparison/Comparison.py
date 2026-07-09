from pathlib import Path


def run_network_comparison(given_path: str,
                           organ_data: str,
                           connection_data: str,
                           threshold: float = 0.3,
                           literature_results_path: str | None = None,
                           llm_descriptions_path: str | None = None,
                           organ_descriptions_path: str | None = None):

    given_path = Path(given_path)
    print(f"[i] Input:     {given_path}")
    print(f"[i] Threshold: {threshold}")

    # ── Load metadata ──────────────────────────────────────────────────────
    import sys
    sys.path.insert(0, str(Path(__file__).parent.parent))

    from Data_Loader import load_data
    node_metadata = load_data.load_node_metadata_from_csv(organ_data)
    edge_metadata = load_data.load_edge_metadata_from_csv(connection_data)

    # ── Literature results ─────────────────────────────────────────────────
    lit_results = {}
    if literature_results_path:
        from Literature_Search.pubmed_search import load_literature_results
        lit_results = load_literature_results(literature_results_path)
        n_papers = sum(1 for v in lit_results.values() if v.get("n_papers_found", 0) > 0)
        print(f"[i] Literature results: {len(lit_results)} edges, {n_papers} with papers.")

    from Literature_Search.pubmed_search import merge_with_edge_metadata
    llm_descriptions = {}
    if llm_descriptions_path:
        from Literature_Search.llm_descriptions import load_llm_descriptions
        llm_descriptions = load_llm_descriptions(llm_descriptions_path)
    merged_edge_data = merge_with_edge_metadata(edge_metadata, lit_results, llm_descriptions)

    # ── Build reference binary matrix ─────────────────────────────────────
    ref = load_data.metadata_dict_to_binary_table(edge_metadata)

    # ── Load & align input matrix ──────────────────────────────────────────
    given = load_data.load_csv(str(given_path))
    from Matrix_Comparison.Alignment import align_matrices
    ref_aligned, given_aligned = align_matrices(ref, given)

    # ── Compare ────────────────────────────────────────────────────────────
    from Matrix_Comparison.Alignment import compare_networks
    comparison_matrix, graph = compare_networks(ref_aligned, given_aligned, threshold)

    organ_descs = {}
    if organ_descriptions_path:
        from Literature_Search.organ_descriptions import load_organ_descriptions
        organ_descs = load_organ_descriptions(organ_descriptions_path)

    for node in graph.nodes:
        organ_entry = organ_descs.get(node, {})
        graph.nodes[node]['description']     = node_metadata.get(node, "")
        graph.nodes[node]['llm_description'] = organ_entry.get("description", "")
        graph.nodes[node]['llm_papers']      = organ_entry.get("papers", [])
    for u, v in graph.edges:
        graph.edges[u, v]['description'] = edge_metadata.get((u, v), "")
        graph.edges[u, v]['merged_data'] = merged_edge_data.get((u, v), {})

    # ── Output paths — alongside the input CSV ─────────────────────────────
    out_folder   = given_path.parent
    out_matrix   = out_folder / f"{given_path.stem}_comparison_matrix.csv"
    out_html     = out_folder / f"{given_path.stem}_comparison.html"

    # Reference network always goes into metabolic_data/ next to source data
    ref_html = Path(connection_data).parent / "reference_network.html"

    comparison_matrix.to_csv(out_matrix)
    print(f"[ok] Comparison matrix: {out_matrix}")

    from Visualisation import networkBuilderUtils

    networkBuilderUtils.export_network_to_cytoscape_dashboard(
        graph=graph,
        filename=str(out_html),
        title="Metabolic Network Comparison",
    )
    print(f"[ok] Comparison network: {out_html}")

    # ── Reference network (all connections, not filtered to input organs) ──
    import networkx as nx
    ref_graph = nx.Graph()

    for organ, desc in node_metadata.items():
        ref_graph.add_node(organ, description=desc)

    for (o1, o2), text in edge_metadata.items():
        if o1 >= o2:
            continue
        ref_graph.add_edge(o1, o2)
        ref_graph.edges[o1, o2]['description'] = text
        ref_graph.edges[o1, o2]['merged_data']  = merged_edge_data.get((o1, o2), {})
        ref_graph.edges[o1, o2]['color'] = "#64748b"

    networkBuilderUtils.export_network_to_cytoscape_dashboard(
        graph=ref_graph,
        filename=str(ref_html),
        include_legend=False,
        title="Metabolic Reference Network",
    )
    print(f"[ok] Reference network: {ref_html}")

    return comparison_matrix
