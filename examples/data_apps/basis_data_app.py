from pathlib import Path

import dash
import pandas as pd
import plotly.express as px
from dash import dcc, html, Input, Output

data_path = Path('alpha_basis_data.csv')

# read in data
data = pd.read_csv(data_path)

app = dash.Dash(__name__, external_stylesheets=['https://codepen.io/chriddyp/pen/bWLwgP.css'])
app.layout = html.Div([
    html.H1('AlphaE Error'),
    dcc.Dropdown(
        id='molecule',
        options=[{'label': i, 'value': i} for i in data.molecule.unique()],
        value='H2O'
    ),
    dcc.RadioItems(
        id='iso_type',
        options=[{'label': i, 'value': i} for i in ['alphaE', 'gammaE', 'alpha', 'gamma']],
        value='alphaE'
    ),
    # Title the sliders
    dcc.RangeSlider(
        id='valence',
        min=0, max=5, step=1,
        marks={i: '{}'.format(i) for i in range(0, 5)},
        value=[0, 5],

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
    [Input('molecule', 'value'),
     Input('iso_type', 'value'),
     Input('omega', 'value'),
     Input('valence', 'value')
     ],

)
def update_figure(molecule, iso_type, omega, valence):
    filtered_data = data.query('molecule==@molecule')
    filtered_data = filtered_data.query('omega>=@omega[0] and omega<=@omega[1]')
    # map the valence to ['D','T', 'Q', '5', '6', '7', '8']
    valence_map = {0: 'D', 1: 'T', 2: 'Q', 3: '5', 4: '6', 5: '7', 6: '8'}
    valence = [valence_map[i] for i in range(valence[0], valence[1] + 1)]
    filtered_data = filtered_data.query('valence in @valence')

    # remove unused valence level catergories
    # make a 3D scatter plot of the alphaE data
    fig = px.scatter_3d(filtered_data, x='valence', y='omega', z=iso_type, color='Type',
                        hover_data=['molecule', 'basis', 'mol_system', 'alpha', 'alphaMRA', 'gamma',
                                    'gammaMRA'
                                    ],
                        color_discrete_sequence=px.colors.qualitative.Safe,
                        opacity=0.8,
                        title='AlphaE Error',
                        )
    fig.update_layout(width=1000, height=1000)  # Adjust these values as needed
    # zlabel to $\alpha\%$ Error
    fig.update_layout(scene_zaxis_title=r'Percent Error')
    return fig


if __name__ == '__main__':
    app.run_server(debug=True)
