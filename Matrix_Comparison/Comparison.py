
def run_network_comparison(reference_path: str, given_path: str, threshold: float = 0.3):
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

    # --- Define output paths in SAME folder as given network ---
    output_folder = given_path.parent
    output_html = output_folder / f"{given_path.stem}_comparison.html"
    output_matrix = output_folder / f"{given_path.stem}_comparison_matrix.csv"

    # Save comparison matrix
    comparison_matrix.to_csv(output_matrix)
    print(f"[✔] Comparison matrix saved to {output_matrix}")

    # Export network HTML
    from Visualisation import networkBuilderUtils
    networkBuilderUtils.export_network_to_html(graph, filename=str(output_html))

    print(f"[✔] HTML visualization saved to {output_html}")

    return comparison_matrix
