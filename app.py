"""Gradio app for antibody repertoire search and alignment."""

from __future__ import annotations

import pandas as pd
import gradio as gr

from database_utils import ALLOWED_REGIONS, align_regions, search_sequences

DB_PATH = "antibody_repertoire.db"



def gradio_search(region: str, motif: str, species: str, chain_type: str) -> pd.DataFrame:
    try:
        species = species.strip() if species else None
        chain_type = chain_type.strip() if chain_type else None
        motif = motif.strip() if motif else ""

        if motif == "":
            return pd.DataFrame()

        return search_sequences(
            db_path=DB_PATH,
            region=region,
            motif=motif,
            species=species,
            chain_type=chain_type,
        )
    except Exception as exc:
        return pd.DataFrame({"error": [str(exc)]})



def gradio_align(seq_id1: str, seq_id2: str, region: str) -> str:
    try:
        return align_regions(
            db_path=DB_PATH,
            seq_id1=seq_id1.strip(),
            seq_id2=seq_id2.strip(),
            region=region,
        )
    except Exception as exc:
        return f"Error: {exc}"


with gr.Blocks() as demo:
    gr.Markdown("# Antibody Repertoire Database Prototype")
    gr.Markdown("Search CDR/framework motifs and align antibody regions.")

    with gr.Tab("Search"):
        region = gr.Dropdown(ALLOWED_REGIONS, value="cdr3", label="Region")
        motif = gr.Textbox(label="Motif / partial sequence")
        species = gr.Textbox(label="Species (optional)")
        chain_type = gr.Textbox(label="Chain type (optional)")
        search_btn = gr.Button("Search")
        search_output = gr.Dataframe()

        search_btn.click(
            fn=gradio_search,
            inputs=[region, motif, species, chain_type],
            outputs=search_output,
        )

    with gr.Tab("Align"):
        seq_id1 = gr.Textbox(label="Sequence ID 1")
        seq_id2 = gr.Textbox(label="Sequence ID 2")
        align_region = gr.Dropdown(ALLOWED_REGIONS, value="cdr3", label="Region")
        align_btn = gr.Button("Align")
        align_output = gr.Textbox(label="Alignment", lines=15)

        align_btn.click(
            fn=gradio_align,
            inputs=[seq_id1, seq_id2, align_region],
            outputs=align_output,
        )


if __name__ == "__main__":
    demo.launch(debug=True)
