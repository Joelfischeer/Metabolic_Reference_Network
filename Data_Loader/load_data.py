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

def load_node_metadata(folder: str):
    """
    Returns dict: {node_name: text}
    """
    metadata = {}
    folder = Path(folder)

    if not folder.exists():
        print("[ℹ] No node metadata folder found.")
        return metadata

    for file in folder.glob("*.txt"):
        node_name = file.stem  # filename without .txt
        metadata[node_name] = file.read_text(encoding="utf-8")

    print(f"[✔] Loaded metadata for {len(metadata)} nodes")
    return metadata


def load_edge_metadata(folder: str):
    """
    Returns dict: {(node1, node2): text}
    Order independent.
    """
    metadata = {}
    folder = Path(folder)

    if not folder.exists():
        print("[ℹ] No edge metadata folder found.")
        return metadata

    for file in folder.glob("*.txt"):
        name = file.stem
        parts = name.split("_")

        if len(parts) != 2:
            continue

        n1, n2 = parts
        text = file.read_text(encoding="utf-8")

        metadata[(n1, n2)] = text
        metadata[(n2, n1)] = text  # make symmetric

    print(f"[✔] Loaded metadata for {len(metadata)//2} edges")
    return metadata
