import pandas as pd
import numpy as np
from pathlib import Path



def load_excel(file_path: str, decimal=","):
    return pd.read_excel(file_path, decimal=decimal)

""" 
def load_csv(file_path: str, decimal=","):
    df = pd.read_csv(file_path, sep=',', decimal=decimal)
    print("Data loaded")
    return df
"""

def load_odf(file_path: str):
    return pd.read_excel(file_path, engine='odf', header=0)



def load_csv(file_path: str) -> pd.DataFrame:
    """
    Load a CSV safely, ensuring all numeric columns are correctly parsed.
    Automatically handles common numeric issues (commas, strings, 'NA', etc.)
    """

    # Step 1: Read raw CSV
    df = pd.read_csv(file_path, dtype=str)  # load everything as strings

    # Step 2: Clean numeric columns
    for col in df.columns:
        # Remove commas, spaces, and non-numeric chars (except '.' and '-')
        df[col] = df[col].str.replace(',', '.', regex=False)   # handle European decimals
        df[col] = df[col].str.replace(r'[^0-9\.\-eE]', '', regex=True)  # keep numeric patterns

        # Convert to numeric where possible
        df[col] = pd.to_numeric(df[col], errors='coerce')

    # Step 3: Drop completely empty columns if needed
    empty_cols = df.columns[df.isna().all()]
    if len(empty_cols) > 0:
        print(f"[⚠] Dropping {len(empty_cols)} completely empty columns: {list(empty_cols)}")
        df = df.drop(columns=empty_cols)

    print(f"[✔] Loaded {df.shape[0]} rows and {df.shape[1]} columns from {file_path}")
    print(f"[ℹ] Numeric columns detected: {df.select_dtypes(include=np.number).shape[1]}")
    return df


def load_node_metadata_from_csv(csv_path: str):
    """
    Returns dict: {node_name: text}
    CSV format:
        column 1 = node_name
        column 2 = text
    Leading/trailing quotes in text are removed.
    """
    metadata = {}
    csv_path = Path(csv_path)

    if not csv_path.exists():
        print("[ℹ] CSV file not found.")
        return metadata

    df = pd.read_csv(csv_path, encoding="utf-8").fillna("")

    # remove leading/trailing quotes from text column (2nd column)
    df.iloc[:, 1] = df.iloc[:, 1].apply(lambda x: x.strip('"') if isinstance(x, str) else x)

    # assume first column = name, second column = text
    for _, row in df.iterrows():
        node_name = str(row.iloc[0])
        text = str(row.iloc[1])
        metadata[node_name] = text

    print(f"[✔] Loaded metadata for {len(metadata)} nodes")
    return metadata


def load_edge_metadata_from_csv(csv_path: str):
    """
    Returns dict: {(node1, node2): text}
    Assumes symmetric matrix where ONLY upper triangle is filled.
    Leading/trailing quotes in text are removed.
    """
    metadata = {}
    csv_path = Path(csv_path)

    if not csv_path.exists():
        print("[ℹ] No edge metadata CSV found.")
        return metadata

    df = pd.read_csv(csv_path, index_col=0, encoding="utf-8").fillna("")

    # remove leading/trailing quotes from all cells
    df = df.fillna("")  # replace NaN with empty string
    df = df.astype(str).apply(lambda col: col.str.strip('"'))

    nodes = list(df.index)

    for i, n1 in enumerate(nodes):
        for j in range(i + 1, len(nodes)):  # only upper triangle
            n2 = nodes[j]
            text = df.loc[n1, n2]

            if str(text).strip() == "":
                continue

            metadata[(n1, n2)] = text
            metadata[(n2, n1)] = text  # make symmetric

    print(f"[✔] Loaded metadata for {len(metadata)//2} edges")
    return metadata

def metadata_dict_to_binary_table(metadata: dict):
    """
    Convert metadata dict {(node1,node2): text} into a binary table:
      - 1 if a text exists in the upper triangle
      - 0 otherwise
    Only upper triangle (row < col) will have 1s; diagonal and lower triangle remain 0.
    
    Returns:
        binary_df: pd.DataFrame
    """
    if not metadata:
        return pd.DataFrame()
    
    # Determine all nodes from dict keys
    nodes = sorted(set([n for edge in metadata.keys() for n in edge]))
    n = len(nodes)
    
    # Create empty DataFrame
    binary_df = pd.DataFrame(0, index=nodes, columns=nodes, dtype=int)
    
    # Fill upper triangle with 1 where text exists
    for (n1, n2), text in metadata.items():
        if text.strip() == "":
            continue
        i = nodes.index(n1)
        j = nodes.index(n2)
        if i < j:  # only upper triangle
            binary_df.iloc[i, j] = 1
    
    print(f"[✔] Binary table created with {binary_df.values.sum()} positive connections (upper triangle only)")
    return binary_df

