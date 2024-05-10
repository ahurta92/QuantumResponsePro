from pathlib import Path

import dash
import pandas as pd
import plotly.express as px
from dash import dcc, html, Input, Output



data_path = Path('cluster_data.csv')
# read in data
data = pd.read_csv(data_path)
data.cluster = data.cluster.astype('category')
data.dropna(inplace=True)
print(data)
vl = ['D', 'T', 'Q', '5']
app = dash.Dash(__name__, external_stylesheets=['https://codepen.io/chriddyp/pen/bWLwgP.css'])
server = app.server
app.layout = html.Div(
    [
        html.Div([
            html.H1('Basis set error'),
            dcc.RadioItems(
                id='iso_type',
                options=[{'label': i, 'value': i} for i in ['alphaE', 'gammaE']],
                value='alphaE'
            ),
            # Title the sliders
            dcc.RangeSlider(
                id='valence',
                min=0, max=3, step=1,
                marks={i: '{}'.format(vl[i]) for i in range(0, 4)},
                value=[0, 3],

            ),
            dcc.RangeSlider(
                id='omega',
                min=0, max=8, step=1,
                marks={i: '{}'.format(i) for i in range(0, 9)},
                value=[0,0]
            ),
            dcc.Dropdown(
                id='cluster',
                options=[{'label': i, 'value': i} for i in data.cluster.unique()],
                value='0'
            ),

            dcc.Graph(id='Cluster Basis Set Error Data'),
        ],
            #       style={'width': '49%', 'display': 'inline-block', 'padding': '0 20'}
        ),
    ])


@app.callback(
    Output('Cluster Basis Set Error Data', 'figure'),
    [
        Input('iso_type', 'value'),
        Input('omega', 'value'),
        Input('valence', 'value'),
        Input('cluster', 'value')
    ],

)
def update_figure(iso_type, omega, valence, cluster):
    filtered_data = data.copy()
    if iso_type == 'gammaE':
        filtered_data = filtered_data.query('gamma.abs() > 1e-3')
    filtered_data = filtered_data.query('omega>=@omega[0] and omega<=@omega[1]')
    filtered_data = filtered_data.query('cluster==@cluster')
    # map the valence to ['D','T', 'Q', '5', '6', '7', '8']
    valence_map = {0: 'D', 1: 'T', 2: 'Q', 3: '5', 4: '6', 5: '7', 6: '8'}
    valence = [valence_map[i] for i in range(valence[0], valence[1] + 1)]
    filtered_data = filtered_data.query('valence in @valence')
    filtered_data.dropna(inplace=True)
    fig = px.strip(filtered_data, x='valence', y=iso_type, color='Type', facet_col_spacing=0.05,
                   width=1000,
                   height=1000,
                   category_orders={'cluster': ['0', '1', '2', '3', '4', '5', '6', ],
                                    'Type':['aug-cc-pVnZ','aug-cc-pCnZ','d-aug-cc-pVnZ','d-aug-cc-pCnZ']
                                    },
                   hover_data=['molecule', 'basis', 'valence', 'omega', 'Type', 'alpha',
                               'alphaMRA',
                               'gamma', 'gammaMRA',
                               ],

                   )
    # add an extra legend with the cluster molecules in it
    # add a Table with the cluster data
    molecule = filtered_data.molecule.unique()
    # add the indiviual molecules to the legend
    # add legened groups





    return fig


if __name__ == '__main__':
    app.run_server(debug=True)
