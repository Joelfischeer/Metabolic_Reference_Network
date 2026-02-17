import pandas as pd
import numpy as np


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


