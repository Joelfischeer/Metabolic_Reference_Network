
def run_network_comparison(given_path: str, 
                           organ_data: str,
                           connection_data: str,
                           threshold: float = 0.3):
    
    from pathlib import Path

    # Convert to Path objects
    given_path = Path(given_path)

    print(f"[ℹ] Using threshold = {threshold}")


    # --- Load metabolic metadata using separate functions ---
    from Data_Loader import load_data
    node_metadata = load_data.load_node_metadata_from_csv(organ_data)
    edge_metadata = load_data.load_edge_metadata_from_csv(connection_data)

    #Convert the edge metadata to a binary matrix:
    ref = load_data.metadata_dict_to_binary_table(edge_metadata)
    # Load given matrix:
    given = load_data.load_csv(str(given_path))

    # Align organs automatically
    from Matrix_Comparison.Alignment import align_matrices
    ref_aligned, given_aligned = align_matrices(ref, given)

    # Compare with threshold
    from Matrix_Comparison.Alignment import compare_networks 
    comparison_matrix, graph = compare_networks(
        ref_aligned,
        given_aligned,
        threshold
    )

    # --- Attach metadata to the graph ---
    for node in graph.nodes:
        graph.nodes[node]['description'] = node_metadata.get(node, "")
    for u, v in graph.edges:
        graph.edges[u, v]['description'] = edge_metadata.get((u, v), "")


    # --- Define output paths in SAME folder as given network ---
    output_folder = given_path.parent
    output_html = output_folder / f"{given_path.stem}_comparison.html"
    output_matrix = output_folder / f"{given_path.stem}_comparison_matrix.csv"

    # Save comparison matrix
    comparison_matrix.to_csv(output_matrix)
    print(f"[✔] Comparison matrix saved to {output_matrix}")

    # Export network HTML
    from Visualisation import networkBuilderUtils
    networkBuilderUtils.export_network_to_cytoscape_dashboard(graph=graph,
                                                filename=str(output_html)
                                            )

    print(f"[✔] Comparison Network saved to {output_html}")

     # --- Export reference network HTML ---
    reference_graph = graph.__class__()  # create same type of graph (Graph or DiGraph)

    # Add nodes
    for node in ref_aligned.index:
        reference_graph.add_node(node)
        reference_graph.nodes[node]['description'] = node_metadata.get(node, "")

    # Add edges exactly as in reference
    for u, v in ref_aligned.stack().items():
        if v == 1:
            reference_graph.add_edge(u[0], u[1])
            reference_graph.edges[u[0], u[1]]['description'] = edge_metadata.get((u[0], u[1]), "")
            reference_graph.edges[u[0], u[1]]['color'] = "green"  # or keep from metadata if exists

    # Export reference network HTML
    reference_html = output_folder / f"reference_network.html"
    networkBuilderUtils.export_network_to_cytoscape_dashboard(
        graph=reference_graph,
        filename=str(reference_html),
        include_legend=False
    )
    print(f"[✔] Reference network HTML visualization saved to {reference_html}")

    return comparison_matrix
