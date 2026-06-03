import dash_bootstrap_components as dbc
from dash import dcc, html, Input, Output, State, callback_context, callback, register_page, dash_table, no_update
import pandas as pd
import pathlib
import os
import io
import hashlib
import subprocess

register_page(__name__, "/target_prediction")

################
# PATH
################
# Project path
BASE_PATH = pathlib.Path(__file__).parent.parent.parent.resolve()
# Dashboard path 
DASHBOARD_PATH = BASE_PATH.joinpath("dashboard")
# Input path 
INPUT_PATH = BASE_PATH.joinpath("input")
# Output path
OUTPUT_PATH = BASE_PATH.joinpath("output")

#################
# IMPORT DATASETS
################# 
# Read the meta file to get details of experimental design
metadata_df = pd.read_csv(DASHBOARD_PATH.joinpath("metadata.csv"))

# species-groups dict from metadata, for dropdown 
groups_by_species = metadata_df.groupby('species').apply(lambda x: set(x['group'])).to_dict()

# Species-alias dict from metadata, for dropdown 
species_list = dict(zip(metadata_df['species'].unique(), metadata_df['alias'].unique()))

# 12 isomiR types 
isomir_types = [
    'mirna_exact', 
    'iso_5p-iso_snp-iso_3p', 
    'iso_5p-iso_snp',
    'iso_5p-iso_multi_snp-iso_3p',
    'iso_5p-iso_multi_snp',
    'iso_5p-iso_3p',
    'iso_5p_only',
    'iso_snp-iso_3p',
    'iso_snp_only',
    'iso_multi_snp-iso_3p',
    'iso_multi_snp_only',
    'iso_3p_only'
]

#################
# UI core components
#################
# Side bar - Descritipon
def description_card():
    """

    :return: A Div containing dashboard title & descriptions.
    """
    return html.Div(
        id="description-card",
        children=[
            html.H4("IsomiR Target Prediction")
        ],
    )

# Side bar - filter 
def generate_control_card():
    """

    :return: A Div containing controls for graphs.
    """
    return html.Div(
        id="control-card",
        children=[            
            html.P("Select species", className="select-title"),
            dcc.Dropdown(
                id="species-select-target",
                options=[{"label": v, "value": k} for k,v in species_list.items()],
                multi=False,
            ),
            html.Br(),
            html.P("Select groups", className="select-title"),
            dcc.Dropdown(
                id="group-select-target",
                multi=False,
            ),  
            html.Br(),
            html.P("Select canonical miRNA", className="select-title"),
            dcc.Dropdown(
                id="canonical-select",
                multi=True,
            ), 
            html.Br(),
            html.P("Select isomiR type", className="select-title"),
            dcc.Dropdown(
                id="isomir-type-select",
                options=[{"label": v, "value": v} for v in isomir_types],
                multi=True,
            ), 
            html.Br(),

            html.Div(
                id="btn-container",
                children=[
                    html.Button('Predict targets', id='export-btn', className='target-btn'),
                    dcc.Loading(
                        type="circle",
                        color="#169873",
                        children=[
                            dbc.Modal(
                                [
                                    dbc.ModalHeader(dbc.ModalTitle("Export")),
                                    dbc.ModalBody(id="modal-body-msg", children="Target prediction done. Please check the /output/<species code>/9_target_prediction folders !"),
                                    dbc.ModalFooter(
                                        dbc.Button("Close", id="close", className="ms-auto", n_clicks=0, color="success")
                                    ),
                                ],
                                id="modal-target",
                                is_open=False,
                            )
                        ]
                    ),
                    html.Button('Visualise', id='visualise-btn', className='target-btn', n_clicks=0)  
                ]
            ),
        ],
    )

# MAIN
layout = html.Div(
        id="app-container",
        children=[
            # Storage 
            dcc.Store(id="data"),

            # Side bar
            html.Div(
                id="left-column",
                className="three columns",
                children=[description_card(), generate_control_card()]
            ),
            # Right column
            html.Div(
                id="right-column-target",
                className="nine columns",
                style={
                    "overflowY": "scroll"
                }
            )
        ],
    )

def create_isomirs_fasta(data, selected_canonical, selected_isomir_type, output_path): 
    data = data[data['mirna_name'].isin(selected_canonical)]
    data = data[data['type'].isin(selected_isomir_type)]

    # Create folder if not exist 
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # Create fasta file 
    with open(f"{output_path}/isomiRs.fa", "w+") as fa_file: 
        for _, r in data.iterrows(): 
            fa_file.write(f">{r['annotation'].replace('U', 'T')} {r['type']}\n{r['tag_sequence'].replace('U', 'T')}\n")

def get_miranda_output_path(selected_species, selected_group, selected_canonical, selected_isomir_type):
    sorted_canonical = sorted(selected_canonical)
    sorted_types = sorted(selected_isomir_type)

    output_path = f"{OUTPUT_PATH}/{selected_species}/9_target_prediction/{selected_group}/{'+'.join(selected_canonical)}/{'+'.join(selected_isomir_type)}"
    longest_path = f"{output_path}/skipped_invalid_UTR.txt"
    path_parts = pathlib.Path(longest_path).parts
    long_parts = [part for part in path_parts if len(part.encode()) > 240]

    if len(longest_path.encode()) > 240 or long_parts:
        raise RuntimeError(
            "The selected canonical miRNA and isomiR types create an output path that is too long. "
            "Please select fewer canonical miRNAs or fewer isomiR types, then run target prediction again."
        )

    return output_path

def create_clean_utr_fasta(input_utr_path, output_path):
    clean_utr_path = f"{output_path}/UTR.cleaned.fa"
    skipped_log_path = f"{output_path}/skipped_invalid_UTR.txt"
    valid_bases = set("ACGTUNacgtun")
    skipped_records = []
    written_count = 0

    def write_record(clean_file, header, sequence_lines):
        nonlocal written_count

        if not header:
            return

        sequence = "".join(line.strip() for line in sequence_lines)
        invalid_chars = sorted(set(sequence) - valid_bases)

        if not sequence or invalid_chars:
            skipped_records.append((header, "".join(invalid_chars) or "empty sequence"))
            return

        clean_file.write(f"{header}\n")
        for line in sequence_lines:
            clean_file.write(f"{line.strip().replace('U', 'T').replace('u', 't')}\n")
        written_count += 1

    with open(input_utr_path) as input_file, open(clean_utr_path, "w") as clean_file:
        header = None
        sequence_lines = []

        for line in input_file:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                write_record(clean_file, header, sequence_lines)
                header = line
                sequence_lines = []
            else:
                sequence_lines.append(line)

        write_record(clean_file, header, sequence_lines)

    if skipped_records:
        with open(skipped_log_path, "w") as skipped_log:
            skipped_log.write("header\tinvalid_characters\n")
            for header, invalid_chars in skipped_records:
                skipped_log.write(f"{header[1:]}\t{invalid_chars}\n")

        print(
            f"Skipped {len(skipped_records)} invalid UTR record(s). "
            f"Details written to {skipped_log_path}",
            flush=True
        )
    elif os.path.exists(skipped_log_path):
        os.remove(skipped_log_path)

    if written_count == 0:
        raise RuntimeError("No valid UTR records found after cleaning.")

    return clean_utr_path

def parse_miranda_output(output_path):
    """Parse miranda original_output.txt into perTranscript.txt and perHit.txt using Python."""
    transcript_header = "Seq1\tSeq2\tTot Score\tTot Energy\tMax Score\tMax Energy\tStrand\tLen1\tLen2\tPositions\n"
    hit_header = "Seq1\tSeq2\tScore\tEnergy\tSeq1 Start\tSeq1 End\tSeq2 Start\tSeq2 End\tLen\tSeq1 Identity %\tSeq2 Identity %\n"

    transcript_lines = []
    hit_lines = []

    with open(f"{output_path}/original_output.txt") as f:
        for line in f:
            if line.startswith('>>'):
                transcript_lines.append(line[2:])
            elif line.startswith('>'):
                hit_lines.append(line[1:])

    with open(f"{output_path}/perTranscript.txt", 'w') as f:
        f.write(transcript_header)
        f.writelines(transcript_lines)

    with open(f"{output_path}/perHit.txt", 'w') as f:
        f.write(hit_header)
        f.writelines(hit_lines)


def predict_target(data, selected_species, selected_group, selected_canonical, selected_isomir_type):
    # Output path
    output_path = get_miranda_output_path(selected_species, selected_group, selected_canonical, selected_isomir_type)

    # Create isomiRs fasta file filtered by canonical and isomiR type
    create_isomirs_fasta(data, selected_canonical, selected_isomir_type, output_path)

    clean_utr_path = create_clean_utr_fasta(
        f"{INPUT_PATH}/{selected_species}/UTR.fa",
        output_path
    )

    miranda_cmd = (
        f"miranda {os.path.abspath(f'{output_path}/isomiRs.fa')}"
        f" {os.path.abspath(clean_utr_path)}"
        f" -out {os.path.abspath(f'{output_path}/original_output.txt')}"
    )
    result = subprocess.run(miranda_cmd, shell=True, capture_output=True, text=True)
    if result.stderr:
        print("miRanda stderr:", result.stderr[:500])

    original_out = f"{output_path}/original_output.txt"
    if result.returncode != 0:
        error_message = result.stderr.strip() or "miRanda failed without an error message."
        print(f"miRanda failed with exit code {result.returncode}: {error_message}", flush=True)

        if os.path.exists(original_out):
            os.remove(original_out)
            print(f"Deleted partial miRanda output: {original_out}", flush=True)

        raise RuntimeError(f"miRanda failed with exit code {result.returncode}: {error_message}")

    if not os.path.exists(original_out) or os.path.getsize(original_out) == 0:
        raise RuntimeError(f"miRanda produced no output (exit {result.returncode}): {result.stderr.strip()}")

    parse_miranda_output(output_path)

@callback(
    Output("group-select-target", "options"),
    Input("species-select-target", "value"),
    prevent_initial_call=True,
)
def update_group_options(selected_species):
    if not selected_species:
        return []  # Empty options if no species selected

    groups = groups_by_species.get(selected_species, set())

    return [{"label": g, "value": g} for g in sorted(groups)]

@callback(
    [
        Output('data', 'data'),
        Output('canonical-select', 'options'),
    ],
    [
        Input('species-select-target', 'value'),
        Input("group-select-target", 'value')
    ],
    prevent_initial_call=True
) 
def load_data_and_canonical_list(selected_species, selected_group):
    # list of canonical
    mirnas = []
    # group df
    group_df = pd.DataFrame()

    # load data 
    group_df_list = []
    
    if selected_species and selected_group:
        for rep_file in os.listdir(f"{OUTPUT_PATH}/{selected_species}/1_summarised_isomiRs/{selected_group}"):
            rep_df = pd.read_csv(f'{OUTPUT_PATH}/{selected_species}/1_summarised_isomiRs/{selected_group}/{rep_file}')
            rep_df = rep_df[['mirna_name', 'tag_sequence', 'type', 'annotation']]
            group_df_list.append(rep_df)
        group_df = pd.concat(group_df_list, ignore_index=True) if group_df_list else pd.DataFrame()
        group_df = group_df.drop_duplicates()
        
        if not group_df.empty:
            mirnas = group_df['mirna_name'].unique() 

    return group_df.to_json(), mirnas

@callback(
    Output('export-btn', 'disabled'), 
    [
        Input('species-select-target', 'value'),
        Input('group-select-target', 'value'),
        Input('canonical-select', 'value'),
        Input('isomir-type-select', 'value')
    ],
    prevent_initial_call=True
) 
def disable_btn(selected_species, selected_group, selected_canonical, selected_isomir_type):
    if selected_species and selected_group and selected_canonical and selected_isomir_type:
        return False

@callback(
    [Output("modal-target", "is_open"), Output("modal-body-msg", "children")],
    [
        Input('export-btn', 'n_clicks'),
        Input("close", "n_clicks")
    ],
    [
        State('data', 'data'),
        State('species-select-target', 'value'),
        State('group-select-target', 'value'),
        State('canonical-select', 'value'),
        State('isomir-type-select', 'value'),
        State("modal-target", "is_open")
    ],
    prevent_initial_call=True
)
def export(n_clicks_open, n_clicks_close, data, selected_species, selected_group, selected_canonical, selected_isomir_type, is_open):
    ctx = callback_context

    if not ctx.triggered:
        return is_open, no_update

    triggered_id = ctx.triggered[0]['prop_id'].split('.')[0]

    if triggered_id == 'export-btn' and selected_species and selected_group and selected_canonical and selected_isomir_type:
        try:
            predict_target(pd.read_json(io.StringIO(data)), selected_species, selected_group, selected_canonical, selected_isomir_type)
            return True, "Target prediction done. Please check the /output/<species code>/9_target_prediction folders !"
        except Exception as e:
            return True, f"Error running target prediction: {e}"

    elif triggered_id == 'close':
        return False, no_update

    return is_open, no_update

@callback(
    Output('right-column-target', 'children'),
    Input('visualise-btn', 'n_clicks'),
    [
        State('species-select-target', 'value'),
        State('group-select-target', 'value'),
        State('canonical-select', 'value'),
        State('isomir-type-select', 'value')
    ],
    prevent_initial_call=True
)
def visualise_miranda_output(n_clicks, selected_species, selected_groups, selected_canonical, selected_isomir_type):
    if not selected_species or not selected_groups or not selected_canonical or not selected_isomir_type:
        return []
    else: 
        output_path = get_miranda_output_path(selected_species, selected_groups, selected_canonical, selected_isomir_type) + '/perTranscript.txt'
        if not os.path.exists(output_path):
            return html.P(f'Miranda target prediction not found for {selected_species}, {selected_groups}, {selected_canonical}, {selected_isomir_type}. Predict Target first then try again.')
        else:
            targets_df = pd.read_csv(output_path, sep='\t')
            return dash_table.DataTable(
                targets_df.to_dict("records"),
                [{"name": i, "id": i} for i in targets_df.columns],
                filter_action="native",
                filter_options={"placeholder_text": "Filter column..."},
                page_size=20,
                    style_table={
                    "overflowX": "auto",   
                    "width": "100%"       
                }, 
                style_cell={
                    "minWidth": "120px",  
                    "whiteSpace": "normal",  
                    "textAlign": "center",
                },
                style_header={
                    "fontWeight": "bold"
                }
            )
