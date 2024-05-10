from pathlib import Path

import dash
import pandas as pd
import plotly.express as px
from dash import dcc, html, Input, Output

data_path = Path('beta_hrs.csv')

data_path2 = Path('beta_parallel.csv')
# read in datakk
data = pd.read_csv(data_path)
data2 = pd.read_csv(data_path2)

print(data)
vl = ['D', 'T', 'Q', ]
app = dash.Dash(__name__, external_stylesheets=['https://codepen.io/chriddyp/pen/bWLwgP.css'])
server = app.server
app.layout = html.Div([
    html.H1('Basis set error'),
    # Title the sliders
    dcc.RangeSlider(
        id='valence',
        min=0, max=2, step=1,
        marks={i: '{}'.format(vl[i]) for i in range(0, 3)},
        value=[0, 2],
    ),
    dcc.Dropdown(
        id='Beta_Type',
        options=[{'label': i, 'value': i} for i in ['beta_hrs', 'beta_para']],
        value='beta'
    ),

    dcc.Graph(id='betaE_graph')
])


# Add a callback option to control the value of omega
# create a callback to update the valence level


@app.callback(

    Output('betaE_graph', 'figure'),
    [
        Input('valence', 'value'),
        Input('Beta_Type', 'value')
    ],

)
def update_figure( valence, Beta_Type):
    if Beta_Type == 'beta_hrs':
        filtered_data = data
    else:
        filtered_data = data2
    filtered_data = filtered_data.query('b==0 and c==0 ')
    # map the valence to ['D','T', 'Q', '5', '6', '7', '8']
    valence_map = {0: 'D', 1: 'T', 2: 'Q', 3: '5', 4: '6', 5: '7', 6: '8'}
    valence = [valence_map[i] for i in range(valence[0], valence[1] + 1)]
    filtered_data = filtered_data.query('valence in @valence')
    # set the order of color to ['First-row', 'Fluorine', 'Second-row']
    filtered_data = filtered_data.query('mol_system in ["First-row", "Fluorine", "Second-row"]')

    # remove unused valence level catergories
    # make a 3D scatter plot of the alphaE data
    fig = px.strip(filtered_data, x='Type', y='betaE', color='mol_system',
                     facet_col_spacing=0.05,
                   facet_col='valence',
                   category_orders={'mol_system': ['First-row', 'Fluorine', 'Second-row']},
                   hover_data=['molecule', 'basis', 'valence',  'Type', 'Beta',
                               'beta_MRA',
                               ],
                   color_discrete_sequence=px.colors.qualitative.D3,
                   title='Percent Error in Invariants',
                   )

    return fig


if __name__ == '__main__':
    app.run_server(debug=True)
