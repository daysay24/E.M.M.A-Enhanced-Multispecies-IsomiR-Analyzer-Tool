import dash
from dash import Dash, dcc, html, Input, Output, callback, page_container
import dash_bootstrap_components as dbc

app = Dash(
    __name__, 
    use_pages=True,
    suppress_callback_exceptions=True,
    meta_tags=[{"name": "viewport", "content": "width=device-width, initial-scale=1"}],
    external_stylesheets=[dbc.themes.BOOTSTRAP]
)

app.layout = html.Div(
    [
        dcc.Location(id="url"),
        html.Div(
            id="navigator-bar",
            children =
            [
                html.H6("E.M.M.A Enhanced Multispecies isomiR Analyzer Tool"),
                html.Div(
                    [
                        dcc.Dropdown(
                            id="page-selector",
                            options=[
                                {"label": "Statistics Dashboard", "value": "/"},
                                {"label": "Target Prediction", "value": "target_prediction"},
                            ],
                            value="/",
                            persistence=True,
                            clearable=False,
                            searchable=False,
                            style={"minWidth": "200px"},
                        ),
                    ],
                    style={"marginLeft": "auto", "display": "flex", "gap": 10}
                )
            ],
        ),
        page_container,
    ],

) 

@callback(
    Output("url", "pathname"),
    Input("page-selector", "value"),
    prevent_initial_call=True
)
def navigate(value):
    return value

if __name__ == '__main__':
    app.run(debug=True)