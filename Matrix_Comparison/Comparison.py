
def run_network_comparison(reference_path: str, 
                           given_path: str, 
                           metabolic_data_folder: str,
                           threshold: float = 0.3):
    
    from pathlib import Path

    # Convert to Path objects
    reference_path = Path(reference_path)
    given_path = Path(given_path)

    print(f"[ℹ] Using threshold = {threshold}")

    # Load matrices
    from Data_Loader import load_data
    ref = load_data.load_csv(str(reference_path))
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

    # --- Load metabolic metadata using separate functions ---
    node_metadata = load_data.load_node_metadata(f"{metabolic_data_folder}/organ_data")
    edge_metadata = load_data.load_edge_metadata(f"{metabolic_data_folder}/connection_data")

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
    networkBuilderUtils.export_network_to_html(graph=graph,
                                                filename=str(output_html)
                                            )


    print(f"[✔] HTML visualization saved to {output_html}")

    return comparison_matrix
