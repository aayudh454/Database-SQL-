"""Build a SQLite antibody repertoire database from a CSV or Excel file."""

from __future__ import annotations

import argparse
import pathlib
import sqlite3
from typing import Iterable

import pandas as pd

REQUIRED_COLS = [
    "sequence_id",
    "species",
    "chain_type",
    "v_gene",
    "d_gene",
    "j_gene",
    "framework1",
    "cdr1",
    "framework2",
    "cdr2",
    "framework3",
    "cdr3",
    "framework4",
]

REGION_COLS = [
    "framework1",
    "cdr1",
    "framework2",
    "cdr2",
    "framework3",
    "cdr3",
    "framework4",
]


def load_input_table(input_path: str) -> pd.DataFrame:
    """Load CSV or Excel into a pandas DataFrame."""
    path = pathlib.Path(input_path)
    suffix = path.suffix.lower()

    if suffix == ".csv":
        df = pd.read_csv(path)
    elif suffix in {".xlsx", ".xls"}:
        df = pd.read_excel(path)
    else:
        raise ValueError("Input file must be .csv, .xlsx, or .xls")

    return df



def clean_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Standardize headers from notebook-style input files."""
    df = df.copy()
    df.columns = [str(c).replace("_aa", "") for c in df.columns]
    df.columns = [str(c).strip().lower().replace(" ", "_") for c in df.columns]
    return df



def validate_columns(df: pd.DataFrame, required_cols: Iterable[str]) -> None:
    """Raise an error if required columns are missing."""
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")



def build_full_variable_region(df: pd.DataFrame) -> pd.DataFrame:
    """Concatenate framework and CDR regions into a full variable region."""
    df = df.copy()
    df["full_variable_region"] = ""
    for col in REGION_COLS:
        df["full_variable_region"] += df[col].fillna("").astype(str)
    return df



def write_sqlite(df: pd.DataFrame, output_db: str, table_name: str = "antibodies") -> None:
    """Write the DataFrame into SQLite."""
    with sqlite3.connect(output_db) as conn:
        df.to_sql(table_name, conn, if_exists="replace", index=False)



def main() -> None:
    parser = argparse.ArgumentParser(description="Build antibody SQLite database")
    parser.add_argument("--input", required=True, help="Path to CSV or Excel file")
    parser.add_argument("--output", default="antibody_repertoire.db", help="Output SQLite DB path")
    args = parser.parse_args()

    df = load_input_table(args.input)
    df = clean_columns(df)
    validate_columns(df, REQUIRED_COLS)
    df = build_full_variable_region(df)
    write_sqlite(df, args.output)

    print(f"Database created successfully: {args.output}")
    print(f"Rows written: {len(df)}")


if __name__ == "__main__":
    main()
