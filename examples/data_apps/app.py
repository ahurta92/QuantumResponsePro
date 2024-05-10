from pathlib import Path

import dash
import pandas as pd
import plotly.express as px
from dash import dcc, html, Input, Output

data_path = Path('alpha_basis_data.csv')

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
        Input('valence', 'value')
    ],

)
def update_figure(iso_type, omega, valence):
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
    # update ax titles

    # update titles



    return fig


if __name__ == '__main__':
    app.run_server(debug=True)
