import dash
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from dash import dcc, html, Input, Output
from pathlib import Path
from plotly._subplots import SubplotRef
from typing import List, Tuple, TypeAlias

grid_ref_: TypeAlias = List[List[Tuple[SubplotRef]]]


def realign_subplot_axes(
        fig: go.Figure,
        x_axes: bool | dict = False,
        y_axes: bool | dict = False
) -> None:
    """For a given plotly figure, allow all x-axes (column-wise) and/or y-axes (row-wise) to have
    their own ranges. The function operates in-place, modifying the figure without returning it.

    Args:
        fig: The plotly figure to be modified.
        x_axes, y_axes:
            If True, the respective axes will be realigned.
            If a dictionary, they will additionally be updated with the specified settings.

    Examples:

        realign_subplot_axes(fig, y_axes=True)
        realign_subplot_axes(fig, x_axes=dict(showticklabels=True))
    """
    if not (x_axes or y_axes):
        return
    subplot_rows: grid_ref_ = fig._grid_ref  # a list of lists, indicative of the grid dimensions
    n_cols: int = len(subplot_rows[0])  # needed in either case
    if x_axes:
        x_settings = x_axes if isinstance(x_axes, dict) else None
        for col in range(1, n_cols + 1):
            # set the layout.xaxis.matches property of all subplots in the same column to the
            # x-axis name of the first subplot in that column. 'x1' is accepted for 'x' (the very first).
            fig.update_xaxes(x_settings, col=col, matches=f"x{col}")
    if y_axes:
        y_settings = y_axes if isinstance(y_axes, dict) else None
        n_rows = len(subplot_rows)
        for row in range(1, n_rows + 1):
            # set the layout.yaxis.matches property of all subplots in the same row to the
            # y-axis name of the first subplot in that row. Y-axis names are numbered row-wise,
            # such that the first y-axis in the second row is 'y3' (if there are 2 columns).
            first_y = (row - 1) * n_cols + 1
            fig.update_yaxes(y_settings, row=row, matches=f"y{first_y}")


# get working directory path
cwd = Path.cwd()
print(cwd)
data_path = cwd.joinpath('csv', 'alpha_basis_data.csv')
print(cwd)

# read in data
data = pd.read_csv(data_path)
vl = ['D', 'T', 'Q', '5', '6', ]
app = dash.Dash(__name__, external_stylesheets=['https://codepen.io/chriddyp/pen/bWLwgP.css'])
server = app.server
app.layout = html.Div([
    html.H1('Basis set error'),
    dcc.RadioItems(
        id='iso_type',
        options=[{'label': i, 'value': i} for i in ['alphaE', 'gammaE']],
        value='alphaE'
    ),
    dcc.Dropdown(
        id='View',
        options=[{'label': i, 'value': i} for i in ['Static', 'Frequency']],
        value='Type'
    ),
    # Title the sliders
    dcc.RangeSlider(
        id='valence',
        min=0, max=4, step=1,
        marks={i: '{}'.format(vl[i]) for i in range(0, 5)},
        value=[0, 4],

    ),
    dcc.RangeSlider(
        id='omega',
        min=0, max=8, step=1,
        marks={i: '{}'.format(i) for i in range(0, 9)},
        value=[0, 8]
    ),

    dcc.Graph(id='alphaE_graph')
])


# Add a callback option to control the value of omega
# create a callback to update the valence level


@app.callback(

    Output('alphaE_graph', 'figure'),
    [
        Input('iso_type', 'value'),
        Input('omega', 'value'),
        Input('valence', 'value'),
        Input('View', 'value')
    ],

)
def update_figure(iso_type, omega, valence, View):
    filtered_data = data
    if iso_type == 'gammaE':
        filtered_data = filtered_data.query('gamma.abs() > 1e-3')
    filtered_data = filtered_data.query('omega>=@omega[0] and omega<=@omega[1]')
    # map the valence to ['D','T', 'Q', '5', '6', '7', '8']
    valence_map = {0: 'D', 1: 'T', 2: 'Q', 3: '5', 4: '6', 5: '7', 6: '8'}
    valence = [valence_map[i] for i in range(valence[0], valence[1] + 1)]
    filtered_data = filtered_data.query('valence in @valence')
    # set the order of color to ['First-row', 'Fluorine', 'Second-row']
    filtered_data = filtered_data.query('mol_system in ["First-row", "Fluorine", "Second-row"]')

    # remove unused valence level catergories

    # make a 3D scatter plot of the alphaE data
    if View=='Type':
        fig = px.strip(filtered_data, x='Type', y=iso_type, color='mol_system', facet_col_spacing=0.05,
                       facet_col='valence',
                       category_orders={'mol_system': ['First-row', 'Fluorine', 'Second-row']},
                       hover_data=['molecule', 'basis', 'valence', 'omega', 'Type', 'alpha',
                                   'alphaMRA',
                                   'gamma', 'gammaMRA',
                                   ],
                       color_discrete_sequence=px.colors.qualitative.D3,
                       title='Percent Error in Invariants',
                       )
    else:
        fig = px.strip(filtered_data, x='omega', y=iso_type, color='mol_system', facet_col_spacing=0.05,
                       facet_col='Type',
                       facet_row='valence',
                       category_orders={'mol_system': ['First-row', 'Fluorine', 'Second-row']},
                       hover_data=['molecule', 'basis', 'valence', 'omega', 'Type', 'alpha',
                                   'alphaMRA',
                                   'gamma', 'gammaMRA',
                                   ],
                       color_discrete_sequence=px.colors.qualitative.D3,
                       title='Percent Error in Invariants',
                       )
    # update ax titles

    # update titles

    return fig


if __name__ == '__main__':
    app.run_server(debug=True)
