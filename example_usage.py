"""Small examples showing how to use the antibody repertoire database."""

from database_utils import align_regions, get_sequence_by_id, search_all_regions, search_sequences

DB_PATH = "antibody_repertoire.db"

print("\n--- Search one region ---")
result = search_sequences(
    db_path=DB_PATH,
    region="cdr3",
    motif="YY",
    species="canine",
    chain_type="heavy",
)
print(result.head())

print("\n--- Search across all regions ---")
result2 = search_all_regions(db_path=DB_PATH, motif="AAT")
print(result2.head())

print("\n--- Fetch one sequence ---")
print(get_sequence_by_id(db_path=DB_PATH, sequence_id="ABR0013"))

print("\n--- Align two sequences ---")
print(align_regions(db_path=DB_PATH, seq_id1="ABR0005", seq_id2="ABR0010", region="cdr3"))
