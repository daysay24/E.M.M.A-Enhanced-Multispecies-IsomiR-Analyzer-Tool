import dash
import dash_bootstrap_components as dbc
from dash import dcc, html, Input, Output, State, callback_context, callback, register_page 
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

register_page(__name__, "/")

################
# PATH
################
# Base path
BASE_PATH = pathlib.Path(__file__).parent.resolve()
# Data path
DATA_PATH = BASE_PATH.joinpath("../data").resolve()
# Output path
OUTPUT_PATH = '../outputs/graphs'

#################
# IMPORT DATASETS
################# 
# Read the meta file to get details of experimental design
metadata_df = pd.read_csv(DATA_PATH.joinpath("metadata.csv"))

# species-groups dict from metadata, for dropdown 
groups_by_species = metadata_df.groupby('species').apply(lambda x: set(x['group'])).to_dict()

# Species-alias dict from metadata, for dropdown 
species_list = dict(zip(metadata_df['species'].unique(), metadata_df['alias'].unique()))

# Create lists to store DataFrames before concatenation
canonical_isomirs_list = []
isomir_types_list = []
isomir_types_nt_list = []
templated_nontemplated_extended_list = []
nt_extended_list = []
templated_nontemplated_all_list = []

for species in species_list.keys(): 
    # Get the list of outputs
    outputs = os.listdir(DATA_PATH.joinpath(f'{species}'))
    for output in outputs: 
        output_df = pd.read_csv(DATA_PATH.joinpath(f'{species}/{output}'))
        output_df['species'] = species  # Add species column
        
        # Collect DataFrames in lists instead of appending directly
        if output == 'graph_1_data.csv':
            canonical_isomirs_list += [output_df]
        elif output == 'graph_2_data.csv':
            isomir_types_list += [output_df]
        elif output == 'graph_3_data.csv':
            isomir_types_nt_list += [output_df]
        elif output == 'graph_4_data.csv':
            templated_nontemplated_extended_list += [output_df]
        elif output == 'graph_5_data.csv':
            nt_extended_list += [output_df]
        elif output == 'graph_6_data.csv':
            templated_nontemplated_all_list += [output_df]

# Efficient concatenation at the end
canonical_isomirs_df = pd.concat(canonical_isomirs_list, ignore_index=True) if canonical_isomirs_list else pd.DataFrame()
isomir_types_df = pd.concat(isomir_types_list, ignore_index=True) if isomir_types_list else pd.DataFrame()
isomir_types_nt_df = pd.concat(isomir_types_nt_list, ignore_index=True) if isomir_types_nt_list else pd.DataFrame()
templated_nontemplated_extended_df = pd.concat(templated_nontemplated_extended_list, ignore_index=True) if templated_nontemplated_extended_list else pd.DataFrame()
nt_extended_df = pd.concat(nt_extended_list, ignore_index=True) if nt_extended_list else pd.DataFrame()
templated_nontemplated_all_df = pd.concat(templated_nontemplated_all_list, ignore_index=True) if templated_nontemplated_all_list else pd.DataFrame()

# Map analysis type with dataframe
# Analysis type list, for dropdown 
# TODO 
analysis_type_list = {
    "Canonical miRNAs & isomiRs (all groups)": canonical_isomirs_df, 
    "IsomiR types (rpm)": isomir_types_df,
    "IsomiR types (unique tags)": isomir_types_df, 
    "All isomiR types (charactised by nt)": isomir_types_nt_df,
    "3'isomiR types (charactised by nt)": isomir_types_nt_df,
    "5'isomiR types (charactised by nt)": isomir_types_nt_df,
    "Templated vs Non-templated at extended positions (%)": templated_nontemplated_extended_df,
    "Templated vs Non-templated at extended positions (unique tags)": templated_nontemplated_extended_df,
    "Nt characterisation at extended positions (%)": nt_extended_df,
    "Nt characterisation at extended positions (unique tags)": nt_extended_df,
    "Templated vs Non-templated at all positions": templated_nontemplated_all_df
}

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
            html.H4("IsomiR statistics dashboard")
        ]
    )

# Side bar - filter 
def generate_control_card():
    """

    :return: A Div containing controls for graphs.
    """
    return html.Div(
        id="control-card",
        children=[
            # Analysis type
            html.P("Select an analysis type", className="select-title"),
            html.Span(" ⓘ", id="tooltip-analysis-type", style={"cursor": "pointer"}),
            dbc.Tooltip(
                "Select an analysis type",
                target="tooltip-analysis-type",
                placement="right",
                className="info-tooltip"
            ),
            html.Div(id='analysis-type'),
            dcc.Dropdown(
                id="analysis-type-select",
                options=[{"label": i, "value": i} for i in analysis_type_list.keys()],
            ),
            html.Br(),
            
            # Graph type
            html.Div(id='graph-type-container', children=[
                html.P("Select a graph type", className="select-title"),
                html.Span(" ⓘ", id="tooltip-graph-type", style={"cursor": "pointer"}),
                dcc.Dropdown(
                    id="graph-type-select",
                    options=[{"label": i, "value": i} for i in ['pie', 'bar']]
                ),
                dbc.Tooltip(
                    "Select a graph type. This is only available for analysis type IsomiR type.",
                    target="tooltip-graph-type",
                    placement="right",
                    className="info-tooltip"
                ),
                html.Br()
            ]),
            
            # Species
            html.P("Select species", className="select-title"),
            html.Span(" ⓘ", id="tooltip-species", style={"cursor": "pointer"}),
            dbc.Tooltip(
                "Select one or more species",
                target="tooltip-species",
                placement="right",
                className="info-tooltip"
            ),
            dcc.Dropdown(
                id="species-select",
                options=[{"label": v, "value": k} for k,v in species_list.items()],
                multi=True,
            ),
            html.Br(),

            # Groups
            html.P("Select groups", className="select-title"),
            html.Span(" ⓘ", id="tooltip-groups", style={"cursor": "pointer"}),
            dbc.Tooltip(
                "Select one or more groups. Analysis type and species must be selected first. Listed groups are the union of available groups of selected species.",
                target="tooltip-groups",
                placement="right",
                className="info-tooltip"
            ),
            dcc.Dropdown(
                id="group-select",
                multi=True,
            ),  
            html.Br(),

            html.P(id='legend-title'),
            dcc.Checklist(
                id = 'legend-checklist',
                labelStyle={"display": "flex"},
                value=[]
            ),
            html.Br(),

            html.Div(
                id="export-toggle",
                children=[
                    html.Div(
                        children=[
                            html.P('Export', className="select-title"),
                            html.Span(" ⓘ", id="tooltip-export", style={"cursor": "pointer"})
                        ]
                    ),
                    daq.BooleanSwitch(id='export-switch', on=False),
                    dbc.Tooltip(
                        "Only available after an analysis type, species and groups have been selected: (1) Switch on the toggle to enable export (2) Select graphs (3) Select a export file format (4) Download graphs",
                        target="tooltip-export",
                        placement="right",
                        className="info-tooltip"
                    ),
                ]            
            ),

            html.Div(id='export-container', children=[
                dcc.Dropdown(
                    id="select-export-format",
                    options=[{"label": v, "value": v} for v in ['svg', 'jpeg', 'pdf']],
                    value='svg'
                ),
                html.Div(id="export-folder-tree"),
                html.Button('Download', id='export-btn')
            ]),


            dbc.Modal(
                [
                    dbc.ModalHeader(dbc.ModalTitle("Export")),
                    dbc.ModalBody("Figures are exported successfully. Please check the /outputs/graph folders !"),
                    dbc.ModalFooter(
                        dbc.Button("Close", id="close", className="ms-auto", n_clicks=0, color="success")
                    ),
                ],
                id="modal-target",
                is_open=False,
            )
        ],
    )

################
# Style calculation
###############
# Number of graph calculation 
def calculate_number_graphs(selected_analysis_type, selected_species_len, selected_groups_len): 
    if selected_analysis_type in ["Canonical miRNAs & isomiRs (all groups)", "All isomiR types (charactised by nt)",  "3'isomiR types (charactised by nt)", "5'isomiR types (charactised by nt)"]:
        return selected_species_len
    else: 
        return selected_species_len * selected_groups_len

# Size calculation 
def calculate_sizes(selected_analysis_type, selected_species_len, selected_groups_len):
    """
    Calculate sizes of: 
    - species container div (height, width)
    - subplot (height, width)

    Return a dict {}
    """
    # Number of graph 
    graph_number = calculate_number_graphs(selected_analysis_type, selected_species_len, selected_groups_len)
    
    graph_container_width = ""
    graph_container_height = ""
    graph_subplot_width = ""
    graph_subplot_height = "100%"

    if selected_groups_len == 1 or selected_analysis_type in ["Canonical miRNAs & isomiRs (all groups)", "All isomiR types (charactised by nt)",  "3'isomiR types (charactised by nt)", "5'isomiR types (charactised by nt)", "IsomiR types (rpm)", "IsomiR types (unique tags)"]: 
        if graph_number == 1:  
            graph_container_width = "100%"
            graph_container_height = "100%"
        else: 
            graph_container_width = "50%"
            graph_container_height = "50%"
        graph_subplot_width = "100%"
        
    else: 
        if selected_species_len < 3 :
            graph_container_height = str((1 / selected_species_len * 100)) + "%" 
        else:
            graph_container_height = str((1 / 3 * 100)) + "%" 
        graph_container_width = "100%"
        graph_subplot_width = str(1 / 3 * 100) + "%"
    
    return {
        "graph_container_width": "calc(" + graph_container_width + " - 10px)",
        "graph_container_height": "calc(" + graph_container_height + " - 10px)",
        "graph_subplot_width": graph_subplot_width,
        "graph_subplot_height": graph_subplot_height
    }
    

#################
# Graphs generation
#################
# Graph 1 
def generate_individual_graph_1(selected_analysis_type, species, groups, sizes, selected_legend_items, legend_item_color, figures):
    # Load data
    data = analysis_type_list[selected_analysis_type]
    data = data[data['species'] == species]
    data = data[data['group'].isin(groups)]
    data = data[data['type'].isin(selected_legend_items)]

    # Create figure
    fig = go.Figure()

    # Add a line for each 'type'
    for trace_type in data['type'].unique():
        df_sub = data[data['type'] == trace_type]
        fig.add_trace(go.Scatter(
            x=df_sub['group'],
            y=df_sub['rpm'],
            mode='lines+markers',
            name=trace_type,
            line=dict(color=legend_item_color.get(trace_type, 'pink'))  # fallback to default Plotly blue
        ))

    # Update layout
    fig.update_layout(
        autosize=True,
        margin=dict(
            pad=10,
            t=20,
            b=0,
            l=2,
            r=20
        ),
        yaxis_title="Count",
        barmode='stack',
        plot_bgcolor='white',
        showlegend=False
    )

    # Update x axis 
    fig.update_xaxes(
        title="<b>RPM</b>",
        mirror=True,
        ticks='outside',
        showline=True,
        linecolor='black',
    )

    # Update y axis 
    fig.update_yaxes(
        title="<b>Group</b>",
        mirror=True,
        ticks='outside',
        showline=True,
        linecolor='black',
        gridcolor='lightgrey'
    )

    # Figure name 
    fig_name = f'{selected_analysis_type}:{species}:{"_".join(groups)}'
    figures[fig_name] = fig

    return  {
        'id': fig_name,
        'figure': dcc.Graph(
            figure=fig, 
            config={'displayModeBar': False}, 
            style={
                "width": "100%", 
                "height": "100%"
            })
        }

# Graph 2 
def generate_individual_graph_2_pie(selected_analysis_type, species, group, sizes, selected_legend_items, legend_item_color, figures):
    # Load data
    data = analysis_type_list[selected_analysis_type]
    data = data[data['species'] == species]
    data = data[data['group'] == group]
    data = data[data['grouped_type'].isin(selected_legend_items)]

    # Value type
    value_type = ''
    if selected_analysis_type == 'IsomiR types (rpm)': 
        value_type = 'rpm'
    else:
        value_type = 'unique_tag' 
    
    fig = px.pie(
        data,
        values=value_type,
        names='grouped_type',
        color='grouped_type',
        color_discrete_map=legend_item_color
    )

    fig.update_traces(
        textposition='outside',
        insidetextorientation='radial',
        pull=[0]*len(data),  # ensure no slice is pulled out
        showlegend=False
    )

    fig.update_layout(
        autosize=True,
        margin=dict(pad=10, t=20, b=0, l=2, r=20),
        uniformtext_minsize=10,
        uniformtext_mode='hide',  # prevents text overlap
        plot_bgcolor='white'
    )

    # Figure name 
    fig_name = f'{selected_analysis_type}:{species}:{group}'
    figures[fig_name] = fig

    return  {
        'id': fig_name,
        'figure': dcc.Graph(
            figure=fig, 
            config={'displayModeBar': False},
            style={
                "width": "100%", 
                "height": "100%"
            }
        )
    }

# Graph 2 
def generate_individual_graph_2_bar(selected_analysis_type, species, groups, sizes, selected_legend_items, legend_item_color, figures):
    # Load data
    data = analysis_type_list[selected_analysis_type]
    data = data[data['species'] == species]
    data = data[data['group'].isin(groups)]
    data = data[data['grouped_type'].isin(selected_legend_items)]

    # Value type
    value_type = ''
    if selected_analysis_type == 'IsomiR types (rpm)': 
        value_type = 'rpm'
        data['percentage'] = data.groupby('group')['rpm'].transform(lambda x: (x / x.sum()) * 100)
    else:
        value_type = 'unique_tag' 
        data['percentage'] = data.groupby('group')['unique_tag'].transform(lambda x: (x / x.sum()) * 100)


    # Create stacked bar chart
    fig = go.Figure()

    # Create traces 
    traces = []

    # Group records by grouped_type
    df_grouped = data.groupby('grouped_type')
    for grouped_type, group in df_grouped:
        traces.append(go.Bar(
            x=group['group'],
            y=group['percentage'],
            name=grouped_type,
            marker=dict(color=legend_item_color.get(grouped_type, '#636EFA')),
        ))

    # Create the figure
    fig = go.Figure(traces)

    # Update layout
    fig.update_layout(
        autosize=True,
        margin=dict(
            pad=10,
            t=20,
            b=0,
            l=2,
            r=20
        ),
        yaxis=dict(title='<b>Percentage</b>', ticksuffix='%'),
        xaxis=dict(title=''),
        barmode='stack',
        plot_bgcolor='white',
        showlegend=False
    )

    # Update x axis 
    fig.update_xaxes(
        mirror=True,
        ticks='outside',
        showline=True,
        linecolor='black',
    )

    # Update y axis
    fig.update_yaxes(
        mirror=True,
        ticks='outside',
        showline=True,
        linecolor='black',
        gridcolor='lightgrey'
    )
    
    # Figure name 
    fig_name = f'{selected_analysis_type}:{species}:{"_".join(groups)}'
    figures[fig_name] = fig

    return  {
        'id': fig_name,
        'figure': dcc.Graph(
            figure=fig, 
            config={'displayModeBar': False}, 
            style={
                "width": "100%", 
                "height": "100%"
            })
        }

# Graph 3 
def generate_individual_graph_3(selected_analysis_type, species, groups, sizes, selected_legend_items, legend_item_color, figures):
    # Load data 
    data = analysis_type_list[selected_analysis_type]
    data = data[data['species'] == species]
    data = data[data['group'].isin(groups)]
    data = data[data['type_nt'].isin(selected_legend_items)]

    # value type
    if selected_analysis_type == "3'isomiR types (charactised by nt)":
        data = data[data['grouped_type'] == "3'isomiR"]
    elif selected_analysis_type == "5'isomiR types (charactised by nt)":
        data = data[data['grouped_type'] == "5'isomiR"]

    data['percentage'] = data.groupby('group')['rpm'].transform(lambda x: (x / x.sum()) * 100)

    # Create stacked bar chart
    fig = go.Figure()

    # Create traces 
    traces = []

    # Group records by type_nt
    df_grouped = data.groupby('type_nt')
    for type_nt, group in df_grouped:
        traces.append(go.Bar(
            x=group['group'],
            y=group['percentage'],
            name=type_nt,
            marker=dict(color=legend_item_color.get(type_nt, '#636EFA')),
        ))

    # Create the figure
    fig = go.Figure(traces)

    # Update layout
    fig.update_layout(
        autosize=True,
        margin=dict(
            pad=10,
            t=20,
            b=0,
            l=2,
            r=20
        ),
        yaxis=dict(title='<b>Percentage</b>', ticksuffix='%'),
        xaxis=dict(title=''),
        barmode='stack',
        plot_bgcolor='white',
        showlegend=False
    )

    # Update x axis 
    fig.update_xaxes(
        mirror=True,
        ticks='outside',
        showline=True,
        linecolor='black',
    )

    # Update y axis
    fig.update_yaxes(
        mirror=True,
        ticks='outside',
        showline=True,
        linecolor='black',
        gridcolor='lightgrey'
    )
    
    # Figure name 
    fig_name = f'{selected_analysis_type}:{species}:{"_".join(groups)}'
    figures[fig_name] = fig

    return  {
        'id': fig_name,
        'figure': dcc.Graph(
            figure=fig, 
            config={'displayModeBar': False}, 
            style={
                "width": "100%", 
                "height": "100%"
            })
        }

# Graph 4 
def generate_individual_graph_4(selected_analysis_type, species, group, sizes, selected_legend_items, legend_item_color, figures):
    # Load data
    data = analysis_type_list[selected_analysis_type]
    data = data[data['species'] == species]
    data = data[data['group'] == group]
    data = data[data['templated'].isin(selected_legend_items)]

    # Value type 
    value_type = ''
    y_title = ''
    if selected_analysis_type == 'Templated vs Non-templated at extended positions (unique tags)':
        value_type = 'count'
        y_title = '# Unique tags'
    else: 
        value_type = 'percentage'
        y_title = '%'
        data['percentage'] = data.groupby('position')['count'].transform(lambda x: (x / x.sum()) * 100)

    # Create traces for each medal type
    traces = []
    templated_categories = data['templated'].unique()
    for templated_category in templated_categories:
        df_filtered = data[data['templated'] == templated_category]
        traces.append(go.Bar(
            x=df_filtered['position'],
            y=df_filtered[value_type],
            name=templated_category,
            marker=dict(color=legend_item_color.get(templated_category, "#636EFA"))  # Default to blue if not specified
        ))

    # Create the figure
    fig = go.Figure(traces)

    # Update layout
    fig.update_layout(
        autosize=True,
        title=f"<b>{group}</b>",
        margin=dict(
            pad=10,
            t=40,
            b=0,
            l=2,
            r=20
        ),
        yaxis_title=f'<b>{y_title}</b>',
        xaxis_title="<b>Position</b>",
        barmode='stack',
        plot_bgcolor='white',
        showlegend=False
    )

    # Update x axis 
    fig.update_xaxes(
        mirror=True,
        ticks='outside',
        showline=True,
        linecolor='black',
    )

    # Update y axis 
    fig.update_yaxes(
        mirror=True,
        ticks='outside',
        showline=True,
        linecolor='black',
        gridcolor='lightgrey'
    )

    # Figure name 
    fig_name = f'{selected_analysis_type}:{species}:{group}'
    figures[fig_name] = fig

    return  {
        'id': fig_name,
        'figure': dcc.Graph(
            figure=fig, 
            config={'displayModeBar': False}, 
            style={
                "width": "100%", 
                "height": "100%"
            })
        }

# Graph 5
def generate_individual_graph_5(selected_analysis_type, species, group, sizes, selected_legend_items, legend_item_color, figures):
    # Load data
    data = analysis_type_list[selected_analysis_type]
    data = data[data['species'] == species]
    data = data[data['group'] == group]
    data = data[data['nucleotide'].isin(selected_legend_items)]

    # Value type 
    value_type = ''
    y_title = ''
    if selected_analysis_type == 'Nt characterisation at extended positions (unique tags)':
        value_type = 'count'
        y_title = '# Unique tags'
    else: 
        value_type = 'percentage'
        y_title = '%'
        data['percentage'] = data.groupby('position')['count'].transform(lambda x: (x / x.sum()) * 100)

    # Create stacked bar chart
    fig = go.Figure()

    # Create traces 
    traces = []

    # Group records by type_nt
    df_grouped = data.groupby('nucleotide')
    for nucleotide, grouped in df_grouped:
        traces.append(go.Bar(
            x=grouped['position'],
            y=grouped[value_type],
            name=nucleotide,
            marker=dict(color=legend_item_color.get(nucleotide, '#636EFA')),
        ))

    # Create the figure
    fig = go.Figure(traces)

    # Update layout
    fig.update_layout(
        autosize=True,
        title=f"<b>{group}</b>",
        margin=dict(
            pad=10,
            t=40, 
            b=0,
            l=2,
            r=20
        ),
        yaxis=dict(title=f'<b>{y_title}</b>', ticksuffix='%'),
        xaxis=dict(title='<b>Position</b>'),
        legend_title="Type",
        barmode='stack',
        plot_bgcolor='white',
        showlegend=False
    )

    # Update x axis 
    fig.update_xaxes(
        mirror=True,
        ticks='outside',
        showline=True,
        linecolor='black',
    )

    # Update y axis 
    fig.update_yaxes(
        mirror=True,
        ticks='outside',
        showline=True,
        linecolor='black',
        gridcolor='lightgrey'
    )

    fig_name = f'{selected_analysis_type}:{species}:{group}'
    figures[fig_name] = fig

    return  {
        'id': fig_name,
        'figure': dcc.Graph(
            figure=fig, 
            config={'displayModeBar': False}, 
            style={
                "width": "100%", 
                "height": "100%"
            })
        }

# Graph 6
def generate_individual_graph_6(selected_analysis_type, species, group, sizes, selected_legend_items, legend_item_color, figures):
    # Load data
    data = analysis_type_list[selected_analysis_type]
    data = data[data['species'] == species]
    data = data[data['group'] == group]
    data = data[data['templated'].isin(selected_legend_items)]

    # Get the max position 
    position_count = data[data['position'].str.isnumeric()].groupby('position')['count'].sum() 
    position_count.index = position_count.index.astype(int)
    max_position = position_count[position_count > 0].index.max()

    is_numeric = data['position'].str.isnumeric()
    numeric_position_df = data[is_numeric]
    numeric_position_df = numeric_position_df[(numeric_position_df['position'].astype(float) <= max_position)]
    df = pd.concat([data[~is_numeric], numeric_position_df], ignore_index=True)

    # Create traces
    traces = []
    for templated_category, grouped in df.groupby('templated'):
        traces.append(go.Bar(
            x=grouped['position'],
            y=grouped['count'],
            name=templated_category,
            marker=dict(color=legend_item_color.get(templated_category, "#636EFA"))  # Default to blue if not specified
        ))

    # Create the figure
    fig = go.Figure(traces)

    # Update layout
    fig.update_layout(
        autosize=True,
        title=f'<b>{group}</b>',
        margin=dict(
            pad=10,
            t=40,
            b=0,
            l=2,
            r=20
        ),
        yaxis_title="<b># Unique tags</b>",
        xaxis_title="<b>Positions</b>",
        barmode='stack',
        plot_bgcolor='white',
        showlegend=False
    )

    # Update x axis
    fig.update_xaxes(
        mirror=True,
        ticks='outside',
        showline=True,
        linecolor='black',
    )

    # Update y axis
    fig.update_yaxes(
        mirror=True,
        ticks='outside',
        showline=True,
        linecolor='black',
        gridcolor='lightgrey'
    )
    
    fig_name = f'{selected_analysis_type}:{species}:{group}'
    figures[fig_name] = fig

    return  {
        'id': fig_name,
        'figure': dcc.Graph(
            figure=fig, 
            config={'displayModeBar': False}, 
            style={
                "width": "100%", 
                "height": "100%"
            })
        }


# Graphs for a species of a type
def generate_species_graphs(selected_analysis_type, selected_graph_type, species, selected_groups, sizes, selected_legend_items, legend_item_color, figures):
   
    species_graphs = []
    
    if selected_analysis_type == 'Canonical miRNAs & isomiRs (all groups)':
        species_graphs.append(generate_individual_graph_1(selected_analysis_type, species, selected_groups, sizes, selected_legend_items, legend_item_color, figures))
    elif selected_analysis_type in ['IsomiR types (rpm)', 'IsomiR types (unique tags)']:
        if selected_graph_type == "bar":
            species_graphs.append(generate_individual_graph_2_bar(selected_analysis_type, species, selected_groups, sizes, selected_legend_items, legend_item_color, figures))
        else: 
            for group in selected_groups: 
                species_graphs.append(generate_individual_graph_2_pie(selected_analysis_type, species, group, sizes, selected_legend_items, legend_item_color, figures))
    elif selected_analysis_type in ['All isomiR types (charactised by nt)', "3'isomiR types (charactised by nt)", "5'isomiR types (charactised by nt)"]:
        species_graphs.append(generate_individual_graph_3(selected_analysis_type, species, selected_groups, sizes, selected_legend_items, legend_item_color, figures))
    elif selected_analysis_type in ["Templated vs Non-templated at extended positions (%)", 'Templated vs Non-templated at extended positions (unique tags)']:
        for group in selected_groups: 
            species_graphs.append(generate_individual_graph_4(selected_analysis_type, species, group, sizes, selected_legend_items, legend_item_color, figures))
    elif selected_analysis_type in ['Nt characterisation at extended positions (%)', 'Nt characterisation at extended positions (unique tags)']:
        for group in selected_groups: 
            species_graphs.append(generate_individual_graph_5(selected_analysis_type, species, group, sizes, selected_legend_items, legend_item_color, figures))
    elif selected_analysis_type == 'Templated vs Non-templated at all positions':
        for group in selected_groups: 
            species_graphs.append(generate_individual_graph_6(selected_analysis_type, species, group, sizes, selected_legend_items, legend_item_color, figures))

    return species_graphs

# Generate container of subplots 
def generate_graph_subplots(species_graphs, species, sizes):
    options = [{
        "label": [
            species_graph['figure']
        ],
        "value": species_graph['id'],

    } for species_graph in species_graphs]

    return dcc.Checklist(
        id={'type': 'species-graph-checklist', 'index': species},
        className='species-graph-checklist',
        options=options,
        value=[],
        labelStyle={
            "width": sizes['graph_subplot_width'], 
            "height": sizes['graph_subplot_height'],
            "minWidth": sizes['graph_subplot_width']
        }
    )

# Graph container card 
def generate_graph_containers(selected_analysis_type, selected_graph_type, selected_species, selected_groups, selected_legend_items, legend_item_color, figures): 
    species_container_divs = []

    selected_species_len = len(selected_species)
    selected_groups_len = len(selected_groups)

    sizes = calculate_sizes(selected_analysis_type, selected_species_len, selected_groups_len)

    for species in selected_species: 
        species_graphs = generate_species_graphs(selected_analysis_type, selected_graph_type, species, selected_groups, sizes, selected_legend_items, legend_item_color, figures) 

        species_container_divs.append(
            html.Div(
                id=f'{species}_container',
                className='graph_container',
                children=[
                    html.B(species_list[species]),
                    html.Hr()
                ]
                + [generate_graph_subplots(species_graphs, species, sizes)],
                style={
                    "width": sizes["graph_container_width"], 
                    "height": sizes["graph_container_height"],
                    "overflowX": "scroll"
                }
            )        
        )
    return species_container_divs

def generate_colors(item_number): 
    if item_number <= 5: 
        base_colors = ["#A8FCD5","#30A0C5", "#FFA6A6", "#FFD678", "#A7B6FF"]
        return base_colors[:item_number]
    elif item_number <= 12:
        return px.colors.qualitative.Set3[:item_number]
    else: 
        return px.colors.qualitative.Alphabet[24-item_number:-1]

def get_legend_item_color(selected_analysis_type, selected_species, selected_groups):
    if selected_analysis_type == "Canonical miRNAs & isomiRs (all groups)":
        return {
            "Canonical": "#A8FCD5",
            "IsomiR": "#A7B6FF"
        }
    elif selected_analysis_type in ["IsomiR types (rpm)", "IsomiR types (unique tags)"]:
        return {
            "3'isomiR":'#A8FCD5',
            "5'isomiR":'#30A0C5',
            "Both end isomiR":'#FFA6A6',
            "Canonical":'#FFD678',
            "Others": '#A7B6FF'
        }
    elif selected_analysis_type in ["Templated vs Non-templated at extended positions (%)", "Templated vs Non-templated at extended positions (unique tags)", "Templated vs Non-templated at all positions"]:
        return {
            "Templated": "#A8FCD5",
            "Nontemplated": "#A7B6FF"
        }
    elif selected_analysis_type in ["Nt characterisation at extended positions (%)", "Nt characterisation at extended positions (unique tags)"]:
        return {
            "a": "#A8FCD5",
            "c": "#A7B6FF",
            "g": "#FFA6A6",
            "u":"#30A0C5"
        }
    else:
        data = analysis_type_list[selected_analysis_type]
        data = data[data['species'].isin(selected_species)]
        data = data[data['group'].isin(selected_groups)]
        if selected_analysis_type == "3'isomiR types (charactised by nt)":
            data = data[data['grouped_type'] == "3'isomiR"]
        elif selected_analysis_type == "5'isomiR types (charactised by nt)":
            data = data[data['grouped_type'] == "5'isomiR"]

        legend_items = list(data['type_nt'].unique())
        colors = generate_colors(len(legend_items))
        
        return dict(zip(legend_items, colors))

def get_legend_title(selected_analysis_type):
    if selected_analysis_type in ["Templated vs Non-templated at extended positions (%)", "Templated vs Non-templated at extended positions (unique tags)", "Templated vs Non-templated at all positions"]: 
        return 'Templated'
    elif selected_analysis_type in ["Nt characterisation at extended positions (%)", "Nt characterisation at extended positions (unique tags)"] :
        return 'Nucleotide'
    else:
        return 'Variation type'
    
def export_figures(selected_analysis_type, selected_figures, stored_figures, selected_format): 
    if not os.path.exists(OUTPUT_PATH):
        os.makedirs(OUTPUT_PATH)

    if stored_figures and selected_figures:    
        current_analysis_type_folders = [folder for folder in os.listdir(OUTPUT_PATH)
            if os.path.isdir(os.path.join(OUTPUT_PATH, folder)) and folder.startswith(selected_analysis_type)]
        
        analysis_type_folder_path = f'{OUTPUT_PATH}/{selected_analysis_type}' if len(current_analysis_type_folders) == 0 else f'{OUTPUT_PATH}/{selected_analysis_type} ({len(current_analysis_type_folders)})'
        os.makedirs(analysis_type_folder_path)

        for selected_figure_set in selected_figures: 
            for selected_figure in selected_figure_set:
                selected_figure_names = selected_figure.split(':')
                species = selected_figure_names[1]
                group = selected_figure_names[2]
                fig = stored_figures[selected_figure]
            
                figure_path = f'{analysis_type_folder_path}/{species}/{group}'
                os.makedirs(figure_path)
                pio.write_image(fig, f'{figure_path}/Figure.{selected_format}', format=selected_format)

def generate_folder_tree(selected_analysis_type, selected_figures, selected_graph_type):

    species_subfolders = []
    for selected_figure_set in selected_figures:
        if len(selected_figure_set) > 0:
            species = selected_figure_set[0].split(":")[1]
            group_subfolders = []
            
            for selected_figure in selected_figure_set: 
                group = selected_figure.split(":")[2]
                group_subfolders.append(
                    html.Div([html.Span(group, className="folder fa fa-folder-o"), html.Span(f"Figure.svg", className="file fa fa-file-excel-o")], className="foldercontainer")
                )
            species_subfolders.append(
                html.Div([html.Span(species, className="folder fa fa-folder-o")] + group_subfolders, className="foldercontainer")
            )
    
    return [
        html.Div(
            [
                html.Span(selected_analysis_type, className="folder fa fa-folder-o", **{"data-isexpanded": "true"}),
            ] + species_subfolders, 
            className="foldercontainer")
    ]
        
# MAIN
layout = html.Div(
        id="app-container",
        children=[
            # Storage 
            dcc.Store(id="legend-item-color"),
            dcc.Store(id="stored-figures"),
            # Side bar
            html.Div(
                id="left-column",
                className="three columns",
                children=[
                    description_card(), 
                    generate_control_card()
                ]
            ),
            # Right column
            html.Div(
                id="right-column",
                className="nine columns",
                style={
                    "overflowY": "scroll"
                }
            )
        ],
    )


@callback(
    Output("analysis-type", "children"),
    Input("analysis-type-select", "value"),
    prevent_initial_call=True
)
def update_analysis_type(selected_analysis_type):
    if not selected_analysis_type: 
        return ''
    return selected_analysis_type

@callback(
    Output("group-select", "options"),
    Input("species-select", "value"),
    prevent_initial_call=True
)
def update_group_options(selected_species):
    if not selected_species:
        return []  # Empty options if no species selected
    
    # Get the intersection of all groups for selected species
    available_groups = set()
    for species in selected_species: 
        groups = groups_by_species.get(species, set())
        if len(available_groups) == 0:
            available_groups = groups
        else:
            available_groups = available_groups & groups

    return [{"label": g, "value": g} for g in sorted(available_groups)]

@callback(
    [   Output('legend-checklist', 'options'),
        Output('legend-checklist', 'value'),
        Output('legend-item-color', 'data'),
        Output('legend-title', 'children')
    ],
    [
        Input('analysis-type-select', 'value'),
        Input('species-select', 'value'),
        Input('group-select', 'value'),
    ],
    prevent_initial_call=True    
)
def update_legend_checklist(selected_analysis_type, selected_species, selected_groups):
    
    if not selected_analysis_type or not selected_species or not selected_groups:
        return [], [], {}, ''
    
    legend_item_color = get_legend_item_color(selected_analysis_type, selected_species, selected_groups)
    legend_title = get_legend_title(selected_analysis_type)

    return [
        {
            "label": [
                html.Div(
                    className="legend-indicator",
                    style={
                        "minWidth": "15px",
                        "height": "15px",
                        "backgroundColor": color,
                        "borderRadius": "3px"
                    }
                ),
                html.Span(legend_item, style={"font-size": 15, "padding-left": 5}),
            ],
            "value": legend_item,
        } 
        for legend_item, color in legend_item_color.items()
    ], list(legend_item_color.keys()), legend_item_color, legend_title

@callback(
    Output('right-column', 'children'),
    Output('stored-figures', 'data'),
    [
        Input('analysis-type-select', 'value'),
        Input('graph-type-select', 'value'),
        Input('species-select', 'value'),
        Input('group-select', 'value'),
        Input('legend-checklist', 'value'),
        Input('legend-item-color', 'data')
    ],
    prevent_initial_call=True
)
def update_graphs(selected_analysis_type, selected_graph_type, selected_species, selected_groups, selected_legend_items, legend_item_color):
    figures = {}
    if not selected_species or not selected_groups or not selected_analysis_type:
        return [], []

    if selected_analysis_type in ["IsomiR types (rpm)", "IsomiR types (unique tags)"] and not selected_graph_type:
        return [], []
        
    results = generate_graph_containers(selected_analysis_type, selected_graph_type, selected_species, selected_groups, selected_legend_items, legend_item_color, figures)
    return results, figures

@callback(
    Output("modal-target", "is_open"),
    [
        Input('export-btn', 'n_clicks'),
        Input("close", "n_clicks")
    ],
    [   State('select-export-format', 'value'),
        State('analysis-type-select', 'value'),
        State({'type': 'species-graph-checklist', 'index': ALL}, 'value'),
        State('stored-figures', 'data'),
        State("modal-target", "is_open")
    ],
    prevent_initial_call=True
)
def export(n_clicks_open, n_clicks_close, selected_format, selected_analysis_type, selected_figures, stored_figures, is_open):
    ctx = callback_context  # or `ctx = ctx` if using Dash 2.4+

    if not ctx.triggered:
        return is_open

    triggered_id = ctx.triggered[0]['prop_id'].split('.')[0]

    if triggered_id == 'export-btn' and selected_format:
        export_figures(selected_analysis_type, selected_figures, stored_figures, selected_format)
        return not is_open

    elif triggered_id == 'close':
        return not is_open

    return is_open

@callback(
    Output('export-switch', 'disabled'),
    Output('export-switch', 'on'),
    Input('stored-figures', 'data'),
    State('export-switch', 'on'),
    prevent_initial_call=True
)
def disable_btn(stored_figures, export_on):
    if not stored_figures:
        return True, False
    else: 
        return False, export_on
        

@callback(
    Output({'type': 'species-graph-checklist', 'index': ALL}, 'className'),
    Output({'type': 'species-graph-checklist', 'index': ALL}, 'value'),
    Output('export-container', 'style'),
    [
        Input('export-switch', 'on'),
        Input('export-switch', 'disabled'),
        State('species-select', 'value')
    ]
)
def toggle_checklist(export_on, export_disabled,  selected_species):
    export_container_style = {'display': 'block' if export_on else 'none'}
    classes = []
    values = []

    if not export_disabled and selected_species:
        classes = ['species-graph-checklist' if export_on else 'species-graph-checklist export-off'] * len(selected_species)
        values = [[] for _ in selected_species]
    
    return classes, values, export_container_style

@callback(
    Output('export-folder-tree', 'children'),
    [
        Input('export-switch', 'on'),
        Input({'type': 'species-graph-checklist', 'index': ALL}, 'value'),
        Input('graph-type-select', 'value'),
        State('analysis-type-select', 'value')
    ]
)
def show_folder_tree(export_on, selected_figures, selected_graph_type, selected_analysis_type):
    if export_on == True: 
        return generate_folder_tree(selected_analysis_type, selected_figures, selected_graph_type)
    return []
        
@callback(
    Output('graph-type-container', 'style'),
    Input('analysis-type-select', 'value'),
)
def show_graph_type_select(selected_analysis_type):
    if selected_analysis_type in ["IsomiR types (rpm)", "IsomiR types (unique tags)"]:
        return {'display': 'block'}
    else: 
        return {'display': 'none'}








