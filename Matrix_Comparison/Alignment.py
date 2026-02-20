import numpy as np
import pandas as pd
import networkx as nx


def align_matrices(reference: pd.DataFrame, given: pd.DataFrame):
    """
    Align two adjacency matrices by common organs (row/column names).
    Automatically removes organs not shared.
    """

    # Ensure same index/column format
    reference = reference.copy()
    given = given.copy()

    reference.index = reference.columns
    given.index = given.columns

    common_organs = reference.columns.intersection(given.columns)

    ref_aligned = reference.loc[common_organs, common_organs]
    given_aligned = given.loc[common_organs, common_organs]

    print(f"[✔] Comparing {len(common_organs)} common organs")

    return ref_aligned, given_aligned


def compare_networks(ref: pd.DataFrame,
                     given: pd.DataFrame,
                     threshold: float,
                     node_metadata=None,
                     edge_metadata=None):
    """
    Compare two aligned adjacency matrices using:

    Reference:
        edge exists if value == 1

    Given:
        edge exists if value >= threshold

    Returns:
        comparison_matrix
        colored NetworkX graph
    """

    if node_metadata is None:
        node_metadata = {}

    if edge_metadata is None:
        edge_metadata = {}


    comparison = pd.DataFrame(
        0,
        index=ref.index,
        columns=ref.columns
    )

    G = nx.Graph()

    for node in ref.index:
        description = node_metadata.get(node, "")
        G.add_node(node, description=description)


    for i in ref.index:
        for j in ref.columns:
            if i >= j:
                continue  # undirected

            # --- Binary interpretation ---
            ref_present = pd.notna(ref.loc[i, j]) and ref.loc[i, j] == 1
            given_present = pd.notna(given.loc[i, j]) and given.loc[i, j] >= threshold

            if ref_present and given_present:
                comparison.loc[i, j] = 1
                comparison.loc[j, i] = 1
                edge_text = edge_metadata.get((i, j), "")
                G.add_edge(i, j, weight=1, color="green", description=edge_text)    

            elif not ref_present and given_present:
                comparison.loc[i, j] = -1
                comparison.loc[j, i] = -1
                edge_text = edge_metadata.get((i, j), "")
                G.add_edge(i, j, weight=-1, color="red", description=edge_text)   

            elif ref_present and not given_present:
                comparison.loc[i, j] = 2
                comparison.loc[j, i] = 2
                G.add_edge(i, j, weight=2, color="orange")

    return comparison, G


