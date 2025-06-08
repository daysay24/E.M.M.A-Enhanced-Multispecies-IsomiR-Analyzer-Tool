import dash

import dash_bootstrap_components as dbc
from dash import dcc, html, Input, Output, State, callback_context 
from dash.dependencies import ALL
import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio
import dash_daq as daq

import numpy as np
import pandas as pd
import datetime
from datetime import datetime as dt
import pathlib
import os 

# Create ar Dash instance 
app = dash.Dash(
    __name__,
    meta_tags=[{"name": "viewport", "content": "width=device-width, initial-scale=1"}],
    external_stylesheets=[dbc.themes.BOOTSTRAP]
)
app.title = "E.M.M.A Enhanced Multispecies isomiR Analyzer Tool"
server = app.server
app.config.suppress_callback_exceptions = True

################
# PATH
################
# Base path
BASE_PATH = pathlib.Path(__file__).parent.resolve()
# Data path
DATA_PATH = BASE_PATH.joinpath("data").resolve()
# Output paths
OUTPUT_1_PATH = '../outputs/1_summarised_isomiRs'
OUTPUT_TARGET_PREDICTION = '../outputs/10_target_prediction'

#################
# IMPORT DATASETS
################# 
# Read the meta file to get details of experimental design
metadata_df = pd.read_csv(DATA_PATH.joinpath("metadata.csv"))

# species-groups dict from metadata, for dropdown 
groups_by_species = metadata_df.groupby('species').apply(lambda x: set(x['group'])).to_dict()

# Species-alias dict from metadata, for dropdown 
species_list = dict(zip(metadata_df['species'].unique(), metadata_df['alias'].unique()))

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
            html.H6("E.M.M.A Enhanced Multispecies isomiR Analyzer Tool"),
            html.H3("IsomiR Target Prediction")
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
            html.P("Select species"),
            dcc.Dropdown(
                id="species-select",
                options=[{"label": v, "value": k} for k,v in species_list.items()],
                multi=False,
            ),
            html.Br(),
            html.P("Select groups"),
            dcc.Dropdown(
                id="group-select",
                multi=False,
            ),  
            html.Br(),
            html.P("Select canonical miRNA"),
            dcc.Dropdown(
                id="canonical-select",
                multi=True,
            ), 
            html.Br(),
            html.P("Select isomiR type"),
            dcc.Dropdown(
                id="isomir-type-select",
                options=[{"label": v, "value": v} for v in ["3'isomiR", "5'isomiR", "Both end isomiR"]],
                multi=False,
            ), 
            html.Br(),

            html.Button('Download', id='export-btn'),

            dbc.Modal(
                [
                    dbc.ModalHeader(dbc.ModalTitle("Export")),
                    dbc.ModalBody("Figures are exported successfully. Please check the /outputs/graph folders !"),
                    dbc.ModalFooter(
                        dbc.Button("Close", id="close", className="ms-auto", n_clicks=0, color="success")
                    ),
                ],
                id="modal",
                is_open=False,
            )
        ],
    )

# MAIN
app.layout = html.Div(
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
            id="right-column",
            className="nine columns",
            style={
                "overflowY": "scroll"
            }
        ),
    ],
)

def predict_target(data, selected_species, selected_group, selected_canonical, selected_isomir_type):
    data = data[data['species'] == selected_species]
    data = data[data['group'] == selected_group]
    data = data[data['mirna_name'].isin(selected_canonical)]
    data = data[data['type'].isin(selected_isomir_type)]

    # Run target prediction here 
    print(data)


# TODO update groups option based on speicies 
@app.callback(
    Output("group-select", "options"),
    Input("species-select", "value")
)
def update_group_options(selected_species):
    if not selected_species:
        return []  # Empty options if no species selected

    groups = groups_by_species.get(selected_species, set())

    return [{"label": g, "value": g} for g in sorted(groups)]

# TODO update canonical option based on species + group 
@app.callback(
    [
        Output('data', 'data'),
        Output('canonical-select', 'options'),
    ],
    [
        Input('species-select', 'value'),
        Input("group-select", 'value')
    ]
) 
def load_data_and_canonical_list(selected_species, selected_group):
    # list of canonical
    mirnas = []
    # group 
    group_df = pd.DataFrame()

    # load data 
    group_df_list = []
    
    if selected_species and selected_group:
        for rep_file in os.listdir(f"{OUTPUT_1_PATH}/{selected_species}/{selected_group}"):
            rep_df = pd.read_csv(f'{OUTPUT_1_PATH}/{selected_species}/{selected_group}/{rep_file}', sep='\t')
            print(rep_df)
            # Get replicate name 
            rep_name = rep_file.split('.')[0]
            rep_df['rep_name'] = rep_name
            # Select a subset of important columns 
            rep_df = rep_df[['mirna_name', 'tag_sequence', '#count_tags']]
            group_df_list.append(rep_df)
        group_df = pd.concat(group_df_list, ignore_index=True) if group_df_list else pd.DataFrame()
        
        if not group_df.empty:
            mirnas = group_df['mirna_name'].unique() 

    return group_df.to_json(), mirnas


# TODO enable / disable download 
@app.callback(
    Output('export-btn', 'disabled'), 
    [
        Input('species-select', 'value'),
        Input('group-select', 'value'),
        Input('canonical-select', 'value'),
        Input('isomir-type-select', 'value')
    ]
) 
def disable_btn(selected_species, selected_group, selected_canonical, selected_isomir_type):
    if selected_species and selected_group and selected_canonical and selected_isomir_type:
        return False

# TODO nclick download input: output/1/speices/group1 output output/target/speices/group1 
@app.callback(
    Output("modal", "is_open"),
    [
        Input('export-btn', 'n_clicks'),
        Input("close", "n_clicks")
    ],
    [   
        State('data', 'data'),
        State('species-select', 'value'),
        Input('group-select', 'value'),
        Input('canonical-select', 'value'),
        Input('isomir-type-select', 'value'),
        State("modal", "is_open")
    ]
)
def export(n_clicks_open, n_clicks_close, data, selected_species, selected_group, selected_canonical, selected_isomir_type, is_open):
    ctx = callback_context  # or `ctx = ctx` if using Dash 2.4+

    if not ctx.triggered:
        return is_open

    triggered_id = ctx.triggered[0]['prop_id'].split('.')[0]

    if triggered_id == 'export-btn' and selected_species and selected_group and selected_canonical and selected_isomir_type:
        predict_target(data, selected_species, selected_group, selected_canonical, selected_isomir_type)
        return not is_open

    elif triggered_id == 'close':
        return not is_open

    return is_open

# Run the server
if __name__ == "__main__":
    app.run(debug=True)