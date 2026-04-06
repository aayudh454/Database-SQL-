"""Utilities for searching and aligning an antibody repertoire SQLite database."""

from __future__ import annotations

import sqlite3
from typing import Optional

import pandas as pd
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

ALLOWED_REGIONS = [
    "framework1",
    "cdr1",
    "framework2",
    "cdr2",
    "framework3",
    "cdr3",
    "framework4",
    "full_variable_region",
]



def _validate_region(region: str) -> None:
    if region not in ALLOWED_REGIONS:
        raise ValueError(f"Invalid region: {region}. Allowed: {ALLOWED_REGIONS}")



def search_sequences(
    db_path: str,
    region: str = "cdr3",
    motif: str = "",
    species: Optional[str] = None,
    chain_type: Optional[str] = None,
) -> pd.DataFrame:
    """Search one antibody region for a motif with optional filters."""
    _validate_region(region)

    query = f"SELECT * FROM antibodies WHERE {region} LIKE ?"
    params = [f"%{motif}%"]

    if species:
        query += " AND lower(species) = ?"
        params.append(species.lower())

    if chain_type:
        query += " AND lower(chain_type) = ?"
        params.append(chain_type.lower())

    with sqlite3.connect(db_path) as conn:
        return pd.read_sql_query(query, conn, params=params)



def search_all_regions(
    db_path: str,
    motif: str,
    species: Optional[str] = None,
    chain_type: Optional[str] = None,
) -> pd.DataFrame:
    """Search all framework and CDR regions for a motif."""
    region_list = [
        "framework1",
        "cdr1",
        "framework2",
        "cdr2",
        "framework3",
        "cdr3",
        "framework4",
    ]

    query = "SELECT * FROM antibodies WHERE ("
    query += " OR ".join([f"{r} LIKE ?" for r in region_list])
    query += ")"
    params = [f"%{motif}%" for _ in region_list]

    if species:
        query += " AND lower(species) = ?"
        params.append(species.lower())

    if chain_type:
        query += " AND lower(chain_type) = ?"
        params.append(chain_type.lower())

    with sqlite3.connect(db_path) as conn:
        return pd.read_sql_query(query, conn, params=params)



def get_sequence_by_id(db_path: str, sequence_id: str) -> pd.DataFrame:
    """Retrieve one antibody row by sequence ID."""
    query = "SELECT * FROM antibodies WHERE sequence_id = ?"
    with sqlite3.connect(db_path) as conn:
        return pd.read_sql_query(query, conn, params=[sequence_id])



def align_regions(db_path: str, seq_id1: str, seq_id2: str, region: str = "cdr3") -> str:
    """Align one region between two antibody records."""
    _validate_region(region)

    if region == "full_variable_region":
        # Full variable region alignment is allowed, but still works as a string alignment.
        pass

    df1 = get_sequence_by_id(db_path, seq_id1)
    df2 = get_sequence_by_id(db_path, seq_id2)

    if df1.empty:
        return f"{seq_id1} not found"
    if df2.empty:
        return f"{seq_id2} not found"

    seq1 = str(df1.iloc[0][region]).strip()
    seq2 = str(df2.iloc[0][region]).strip()

    if not seq1 or seq1.lower() == "nan":
        return f"{seq_id1} has no sequence in {region}"
    if not seq2 or seq2.lower() == "nan":
        return f"{seq_id2} has no sequence in {region}"

    alignments = pairwise2.align.globalxx(seq1, seq2)
    if not alignments:
        return "No alignment found"

    return format_alignment(*alignments[0])
