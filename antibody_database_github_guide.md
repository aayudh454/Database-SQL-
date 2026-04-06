# Antibody Repertoire Database: Step-by-Step Guide and Code

This document explains how to build a simple antibody repertoire database from a spreadsheet or CSV file and turn it into a searchable local database with sequence alignment support.

The workflow below follows a practical prototype design for a canine, feline, and equine antibody repertoire database. The goal is to let a user:

- load antibody records from a flat file
- clean and standardize the columns
- validate that all important framework and CDR fields are present
- create a `full_variable_region` sequence
- store the data in an SQLite database
- search frameworks or CDRs by motif
- filter by species and chain type
- retrieve sequences by ID
- align any framework or CDR region between two antibodies

---

## 1. Install the required packages

First install the Python libraries needed for:

- data loading and cleaning: `pandas`, `openpyxl`
- sequence alignment: `biopython`
- optional web interface: `gradio`

```python
!pip install biopython gradio openpyxl pandas
```

---

## 2. Load the antibody repertoire file

This step reads the source dataset into a pandas DataFrame. In this example the file is a CSV called `mock_antibody_repertoire.csv`.

```python
import pandas as pd

file_path = "/content/mock_antibody_repertoire.csv"

df = pd.read_csv(file_path)

print("Shape:", df.shape)
print("Columns:", df.columns.tolist())
df.head()
```

---

## 3. Clean the column names

In many repertoire files, amino acid columns may end with `_aa`. This step removes that suffix. Then we normalize all column names so they are:

- lowercase
- stripped of extra spaces
- converted to underscore format

This makes the database code much more stable.

```python
df.columns = df.columns.str.replace('_aa', '', regex=False)
df.head()
```

```python
df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]
print(df.columns.tolist())
```

---

## 4. Check that all required columns exist

For this database, the minimum fields needed are:

- sequence metadata: `sequence_id`, `species`, `chain_type`
- VDJ annotation: `v_gene`, `d_gene`, `j_gene`
- framework regions: `framework1`, `framework2`, `framework3`, `framework4`
- complementarity determining regions: `cdr1`, `cdr2`, `cdr3`

If any of these are missing, the search and alignment functions will not work correctly.

```python
required_cols = [
    "sequence_id", "species", "chain_type",
    "v_gene", "d_gene", "j_gene",
    "framework1", "cdr1", "framework2", "cdr2",
    "framework3", "cdr3", "framework4"
]

missing = [c for c in required_cols if c not in df.columns]

if missing:
    print("Missing columns:", missing)
else:
    print("All required columns present ✅")
```

---

## 5. Build the full variable region sequence

An antibody variable region can be reconstructed by joining the framework and CDR segments in order:

`FR1 + CDR1 + FR2 + CDR2 + FR3 + CDR3 + FR4`

This gives one complete amino acid string for the variable region, which can be useful for downstream searching, export, and future alignment or clustering tasks.

```python
df["full_variable_region"] = (
    df["framework1"].fillna("") +
    df["cdr1"].fillna("") +
    df["framework2"].fillna("") +
    df["cdr2"].fillna("") +
    df["framework3"].fillna("") +
    df["cdr3"].fillna("") +
    df["framework4"].fillna("")
)

df.head()
```

---

## 6. Create the SQLite database

SQLite is a good starting choice because it is:

- lightweight
- easy to share
- easy to query with SQL
- good for a prototype or local searchable database

This step writes the DataFrame into an SQLite database called `antibody_repertoire.db`.

```python
import sqlite3

conn = sqlite3.connect("antibody_repertoire.db")

df.to_sql("antibodies", conn, if_exists="replace", index=False)

print("Database created successfully ✅")
```

---

## 7. Create a basic motif search function

This function searches one selected region, such as `cdr3` or `framework3`, for a motif. It can also optionally filter by:

- species
- chain type

For example, this allows queries like:

- find all canine antibodies with `YY` in CDR3
- find all feline heavy chains with a motif in framework2

```python
def search_sequences(region="cdr3", motif="", species=None, chain_type=None):
    query = f"SELECT * FROM antibodies WHERE {region} LIKE ?"
    params = [f"%{motif}%"]

    if species:
        query += " AND lower(species) = ?"
        params.append(species.lower())

    if chain_type:
        query += " AND lower(chain_type) = ?"
        params.append(chain_type.lower())

    return pd.read_sql_query(query, conn, params=params)
```

Example:

```python
search_sequences(region="cdr3", motif="YY")
```

---

## 8. Create an advanced search function

A more flexible search function can combine multiple filters such as:

- species
- chain type
- V gene
- J gene
- motif within a chosen region

This is closer to what a user-friendly database backend would need.

```python
def advanced_search(
    species=None,
    chain_type=None,
    v_gene=None,
    j_gene=None,
    region=None,
    motif=None
):
    query = "SELECT * FROM antibodies WHERE 1=1"
    params = []

    if species:
        query += " AND lower(species) = ?"
        params.append(species.lower())

    if chain_type:
        query += " AND lower(chain_type) = ?"
        params.append(chain_type.lower())

    if v_gene:
        query += " AND v_gene LIKE ?"
        params.append(f"%{v_gene}%")

    if j_gene:
        query += " AND j_gene LIKE ?"
        params.append(f"%{j_gene}%")

    if region and motif:
        query += f" AND {region} LIKE ?"
        params.append(f"%{motif}%")

    return pd.read_sql_query(query, conn, params=params)
```

Example:

```python
advanced_search(species="canine", region="cdr3", motif="GG")
```

---

## 9. Retrieve a sequence by ID

Every database should allow direct lookup by sequence identifier. This is useful when a user already knows the antibody record they want to inspect or compare.

```python
def get_sequence_by_id(sequence_id):
    query = "SELECT * FROM antibodies WHERE sequence_id = ?"
    return pd.read_sql_query(query, conn, params=[sequence_id])
```

Example:

```python
get_sequence_by_id("ABR0013")
```

---

## 10. Align a selected region between two antibodies

Alignment is one of the most useful features in an antibody repertoire browser. Here we use Biopython's `pairwise2` to perform a simple global alignment between two selected regions.

This can be applied to:

- `cdr1`
- `cdr2`
- `cdr3`
- `framework1`
- `framework2`
- `framework3`
- `framework4`

```python
!pip install biopython gradio pandas
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

def align_regions(seq_id1, seq_id2, region="cdr3"):
    df1 = get_sequence_by_id(seq_id1)
    df2 = get_sequence_by_id(seq_id2)

    if df1.empty or df2.empty:
        return "Sequence not found"

    seq1 = str(df1.iloc[0][region])
    seq2 = str(df2.iloc[0][region])

    alignments = pairwise2.align.globalxx(seq1, seq2)

    if not alignments:
        return "No alignment"

    return format_alignment(*alignments[0])
```

Example:

```python
print(align_regions("ABR0005", "ABR0010", region="cdr3"))
```

---

## 11. Compute an alignment score

Sometimes you do not need the full formatted alignment. You may only want a similarity score for ranking or comparison. This function returns the score from the best alignment.

```python
def alignment_score(seq_id1, seq_id2, region="cdr3"):
    df1 = get_sequence_by_id(seq_id1)
    df2 = get_sequence_by_id(seq_id2)

    if df1.empty or df2.empty:
        return None

    seq1 = str(df1.iloc[0][region])
    seq2 = str(df2.iloc[0][region])

    alignments = pairwise2.align.globalxx(seq1, seq2)

    if not alignments:
        return None

    return alignments[0].score
```

---

## 12. Search across all frameworks and CDRs at once

Instead of limiting the search to one region, sometimes a user wants to know whether a motif appears anywhere across the full variable region segmentation.

This function checks:

- `framework1`
- `cdr1`
- `framework2`
- `cdr2`
- `framework3`
- `cdr3`
- `framework4`

```python
def search_all_regions(motif, species=None, chain_type=None):
    regions = [
        "framework1", "cdr1", "framework2",
        "cdr2", "framework3", "cdr3", "framework4"
    ]

    query = "SELECT * FROM antibodies WHERE ("
    query += " OR ".join([f"{r} LIKE ?" for r in regions])
    query += ")"

    params = [f"%{motif}%" for _ in regions]

    if species:
        query += " AND lower(species) = ?"
        params.append(species.lower())

    if chain_type:
        query += " AND lower(chain_type) = ?"
        params.append(chain_type.lower())

    return pd.read_sql_query(query, conn, params=params)
```

Example:

```python
search_all_regions("AAT")
```

---

## 13. Download the SQLite database in Google Colab

If you are building this in Google Colab, you can directly download the finished SQLite database file.

```python
from google.colab import files
files.download("antibody_repertoire.db")
```

---

## 14. Store the database path for app development

This is a small but useful step if you later build a frontend app, such as with Gradio, Streamlit, or Flask.

```python
DB_PATH = "antibody_repertoire.db"
```

---

## 15. Optional next step: build a simple Gradio frontend

Once the backend works, a lightweight frontend can be added with tabs for:

- motif search
- advanced filtering
- sequence ID lookup
- framework/CDR alignment

A minimal app structure would look like this:

```python
import gradio as gr


def gradio_search(region, motif, species, chain_type):
    species = species.strip() if species else None
    chain_type = chain_type.strip() if chain_type else None
    motif = motif.strip()

    if not motif:
        return pd.DataFrame()

    result = search_sequences(
        region=region,
        motif=motif,
        species=species if species else None,
        chain_type=chain_type if chain_type else None
    )
    return result


def gradio_align(seq_id1, seq_id2, region):
    return align_regions(seq_id1, seq_id2, region)


search_regions = [
    "framework1", "cdr1", "framework2", "cdr2",
    "framework3", "cdr3", "framework4"
]

with gr.Blocks() as demo:
    gr.Markdown("# Antibody Repertoire Database Prototype")
    gr.Markdown("Search canine, feline, and equine antibody records and align CDR/framework regions.")

    with gr.Tab("Search"):
        region = gr.Dropdown(search_regions, value="cdr3", label="Region")
        motif = gr.Textbox(label="Motif")
        species = gr.Textbox(label="Species")
        chain_type = gr.Textbox(label="Chain Type")
        search_btn = gr.Button("Search")
        search_output = gr.Dataframe()

        search_btn.click(
            fn=gradio_search,
            inputs=[region, motif, species, chain_type],
            outputs=search_output
        )

    with gr.Tab("Align"):
        seq_id1 = gr.Textbox(label="Sequence ID 1")
        seq_id2 = gr.Textbox(label="Sequence ID 2")
        align_region = gr.Dropdown(search_regions, value="cdr3", label="Region")
        align_btn = gr.Button("Align")
        align_output = gr.Textbox(lines=12, label="Alignment")

        align_btn.click(
            fn=gradio_align,
            inputs=[seq_id1, seq_id2, align_region],
            outputs=align_output
        )

demo.launch(debug=True)
```

---

## 16. Summary of the complete workflow

This database prototype follows a simple progression:

1. install libraries  
2. load antibody repertoire data  
3. clean and standardize columns  
4. validate required framework/CDR fields  
5. create a full variable region sequence  
6. save the records to SQLite  
7. enable motif search by region  
8. enable advanced multi-filter search  
9. retrieve sequences by ID  
10. align selected framework or CDR regions  
11. compute alignment scores  
12. search all regions at once  
13. download the database  
14. optionally add a Gradio frontend  

---

## 17. Full code in one place

```python
!pip install biopython gradio openpyxl pandas

import pandas as pd
import sqlite3
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

file_path = "/content/mock_antibody_repertoire.csv"
df = pd.read_csv(file_path)

# Clean columns
(df.columns)
df.columns = df.columns.str.replace('_aa', '', regex=False)
df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]

# Validate columns
required_cols = [
    "sequence_id", "species", "chain_type",
    "v_gene", "d_gene", "j_gene",
    "framework1", "cdr1", "framework2", "cdr2",
    "framework3", "cdr3", "framework4"
]
missing = [c for c in required_cols if c not in df.columns]
if missing:
    raise ValueError(f"Missing columns: {missing}")

# Build full variable region
df["full_variable_region"] = (
    df["framework1"].fillna("") +
    df["cdr1"].fillna("") +
    df["framework2"].fillna("") +
    df["cdr2"].fillna("") +
    df["framework3"].fillna("") +
    df["cdr3"].fillna("") +
    df["framework4"].fillna("")
)

# Save to SQLite
conn = sqlite3.connect("antibody_repertoire.db")
df.to_sql("antibodies", conn, if_exists="replace", index=False)


def search_sequences(region="cdr3", motif="", species=None, chain_type=None):
    query = f"SELECT * FROM antibodies WHERE {region} LIKE ?"
    params = [f"%{motif}%"]

    if species:
        query += " AND lower(species) = ?"
        params.append(species.lower())

    if chain_type:
        query += " AND lower(chain_type) = ?"
        params.append(chain_type.lower())

    return pd.read_sql_query(query, conn, params=params)


def advanced_search(
    species=None,
    chain_type=None,
    v_gene=None,
    j_gene=None,
    region=None,
    motif=None
):
    query = "SELECT * FROM antibodies WHERE 1=1"
    params = []

    if species:
        query += " AND lower(species) = ?"
        params.append(species.lower())

    if chain_type:
        query += " AND lower(chain_type) = ?"
        params.append(chain_type.lower())

    if v_gene:
        query += " AND v_gene LIKE ?"
        params.append(f"%{v_gene}%")

    if j_gene:
        query += " AND j_gene LIKE ?"
        params.append(f"%{j_gene}%")

    if region and motif:
        query += f" AND {region} LIKE ?"
        params.append(f"%{motif}%")

    return pd.read_sql_query(query, conn, params=params)


def get_sequence_by_id(sequence_id):
    query = "SELECT * FROM antibodies WHERE sequence_id = ?"
    return pd.read_sql_query(query, conn, params=[sequence_id])


def align_regions(seq_id1, seq_id2, region="cdr3"):
    df1 = get_sequence_by_id(seq_id1)
    df2 = get_sequence_by_id(seq_id2)

    if df1.empty or df2.empty:
        return "Sequence not found"

    seq1 = str(df1.iloc[0][region])
    seq2 = str(df2.iloc[0][region])

    alignments = pairwise2.align.globalxx(seq1, seq2)

    if not alignments:
        return "No alignment"

    return format_alignment(*alignments[0])


def alignment_score(seq_id1, seq_id2, region="cdr3"):
    df1 = get_sequence_by_id(seq_id1)
    df2 = get_sequence_by_id(seq_id2)

    if df1.empty or df2.empty:
        return None

    seq1 = str(df1.iloc[0][region])
    seq2 = str(df2.iloc[0][region])

    alignments = pairwise2.align.globalxx(seq1, seq2)

    if not alignments:
        return None

    return alignments[0].score


def search_all_regions(motif, species=None, chain_type=None):
    regions = [
        "framework1", "cdr1", "framework2",
        "cdr2", "framework3", "cdr3", "framework4"
    ]

    query = "SELECT * FROM antibodies WHERE ("
    query += " OR ".join([f"{r} LIKE ?" for r in regions])
    query += ")"

    params = [f"%{motif}%" for _ in regions]

    if species:
        query += " AND lower(species) = ?"
        params.append(species.lower())

    if chain_type:
        query += " AND lower(chain_type) = ?"
        params.append(chain_type.lower())

    return pd.read_sql_query(query, conn, params=params)


# Example usage
print(search_sequences(region="cdr3", motif="YY"))
print(advanced_search(species="canine", region="cdr3", motif="GG"))
print(get_sequence_by_id("ABR0013"))
print(align_regions("ABR0005", "ABR0010", region="cdr3"))
print(alignment_score("ABR0005", "ABR0010", region="cdr3"))
print(search_all_regions("AAT"))
```
