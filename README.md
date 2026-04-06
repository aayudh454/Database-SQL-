# Antibody Repertoire Database

A simple GitHub-ready project for building a searchable antibody repertoire database from a CSV or Excel file. This repo is based on the workflow in the uploaded notebook and turns it into a cleaner project structure that is easier to run, reuse, and share.

## What this project does

This database lets you:

- load a canine, feline, and equine antibody repertoire table
- standardize column names
- validate required antibody fields
- build a SQLite database
- search CDRs and frameworks by motif
- filter by species and chain type
- align framework or CDR sequences between two records
- launch a simple Gradio web interface

---

## Repository structure

```text
antibody_database_repo/
├── README.md
├── requirements.txt
├── build_database.py
├── database_utils.py
├── app.py
└── example_usage.py
```

---

## Step-by-step explanation

### Step 1. Prepare the input data
Your source file should contain one row per antibody sequence.

Expected columns:

- `sequence_id`
- `species`
- `chain_type`
- `v_gene`
- `d_gene`
- `j_gene`
- `framework1`
- `cdr1`
- `framework2`
- `cdr2`
- `framework3`
- `cdr3`
- `framework4`

If your original file has headers like `framework1_aa`, this repo automatically removes the `_aa` suffix and standardizes column names.

---

### Step 2. Clean and validate the dataset
The build script:

- loads CSV or Excel
- normalizes column names
- checks whether all required columns are present
- creates a `full_variable_region` field by concatenating framework and CDR segments

This is important because once the database is built, you can search not only individual regions but also the full variable region.

---

### Step 3. Build the SQLite database
The cleaned dataset is written into a SQLite database named `antibody_repertoire.db`.

Why SQLite?

- lightweight
- no server setup needed
- easy to query with Python
- great for a prototype or small production app

---

### Step 4. Search the database
The search utilities allow you to:

- search a single region such as `cdr3`
- search all CDR/framework regions at once
- filter by species
- filter by chain type
- retrieve a sequence by `sequence_id`

Example use cases:

- find all canine antibodies with `YY` in CDR3
- search feline heavy chains for a motif in framework3
- retrieve one specific antibody for comparison

---

### Step 5. Align antibody regions
Using Biopython, the project can align any framework or CDR region between two antibodies.

Example:

- compare CDR3 between `ABR0005` and `ABR0010`
- compare framework2 between two canine sequences

This helps with motif comparison and sequence similarity review.

---

### Step 6. Launch the web app
A simple Gradio app is included.

It has two tabs:

1. **Search**: search motifs in a region with optional filters
2. **Align**: align a selected region between two sequence IDs

This makes the database easier for non-programmers to use.

---

## Installation

```bash
git clone <your-repo-url>
cd antibody_database_repo
python -m venv .venv
source .venv/bin/activate   # Mac/Linux
# .venv\Scripts\activate    # Windows
pip install -r requirements.txt
```

---

## Build the database

### From CSV
```bash
python build_database.py --input mock_antibody_repertoire.csv --output antibody_repertoire.db
```

### From Excel
```bash
python build_database.py --input mock_antibody_repertoire.xlsx --output antibody_repertoire.db
```

---

## Run the Gradio app

```bash
python app.py
```

Then open the local URL shown in the terminal.

---

## Example Python usage

```python
from database_utils import search_sequences, align_regions

result = search_sequences(
    db_path="antibody_repertoire.db",
    region="cdr3",
    motif="YY",
    species="canine",
    chain_type="heavy"
)

print(result.head())

print(align_regions(
    db_path="antibody_repertoire.db",
    seq_id1="ABR0005",
    seq_id2="ABR0010",
    region="cdr3"
))
```

---

## Notes for GitHub

For a real public repo, you should also add:

- `.gitignore`
- a small sample dataset
- screenshots of the app
- a license file
- unit tests

---

## Future improvements

Possible next upgrades:

- FASTA export
- BLAST-style similarity search
- multiple alignment support
- user authentication
- deployment with Docker
- PostgreSQL backend for larger datasets
- filtering by V/D/J gene in the UI
- downloadable search results

