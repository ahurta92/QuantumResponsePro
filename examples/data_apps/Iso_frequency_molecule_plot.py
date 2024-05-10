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
    dcc.Dropdown(
        id='molecule',
        options=[{'label': i, 'value': i} for i in data.molecule.unique()],
        multi=True,
        value=['H2O'],
    ),
    dcc.Checklist(
        id='valence',
        options=[{'label': i, 'value': i} for i in data.valence.unique()],
        value=['D', 'T', 'Q', '5', '6'],
        inline=True
    ),
    dcc.Checklist(
        id='Type',
        options=[{'label': i, 'value': i} for i in data.Type.unique()],
        value=['aug-cc-pVnZ', 'aug-cc-pCVnZ', 'd-aug-cc-pVnZ', 'd-aug-cc-pCVnZ'],
        inline=True
    ),
    dcc.RangeSlider(
        id='omega',
        min=0, max=8, step=1,
        marks={i: '{}'.format(i) for i in range(0, 9)},
        value=[0, 0]
    ),

    dcc.Graph(id='alphaE_graph')
])


# Add a callback option to control the value of omega
# create a callback to update the valence level


@app.callback(

    Output('alphaE_graph', 'figure'),
    [
        Input('molecule', 'value'),
        Input('omega', 'value'),
        Input('valence', 'value'),
        Input('Type', 'value')
    ],

)
def update_figure(molecule, omega, valence, Type):
    filtered_data = data
    iso_type = 'alphaE'
    # filter the data
    filtered_data = filtered_data.query('molecule.isin(@molecule)')

    filtered_data = filtered_data.query('omega>=@omega[0] and omega<=@omega[1]')
    filtered_data = filtered_data.query('valence.isin(@valence)')
    filtered_data = filtered_data.query('Type.isin(@Type)')
    # update ax titles
    fig = px.line(filtered_data, x='valence', y=iso_type, facet_col='ij', color='Type', facet_row='molecule',
                  color_discrete_sequence=px.colors.qualitative.D3,
                  hover_data=['molecule', 'basis', 'valence', 'omega', 'Type', 'alpha', 'alphaMRA'], log_y=False,
                  markers=True,
                  render_mode='svg',
                  line_shape='linear',
                  title='Polarizability Component Percent Error',
                  width=800,
                  height=250 * len(molecule),
                  )

    # fig.update_yaxes(showticklabels=True)

    realign_subplot_axes(fig, y_axes=True)

    fig.update_layout(legend_title_text='Basis Set')
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    # replace y-axis titles
    # fig.update_yaxes(title_text='Percent Error')
    # add horizontal lines at 0 and +/- .05
    fig.add_hline(y=0, line_dash="dot", line_color="black", opacity=.5, )
    # fig.add_hline(y=-.05, line_dash="dot", line_color="black",opacity=.5)

    # update titles

    return fig


if __name__ == '__main__':
    app.run_server(debug=True)
