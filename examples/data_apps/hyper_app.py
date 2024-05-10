from pathlib import Path

import dash
import pandas as pd
import plotly.express as px
from dash import dcc, html, Input, Output

data_path_hrs = Path('beta_hrs.csv')
data_path_para = Path('beta_parallel.csv')

# read in data
data_hrs = pd.read_csv(data_path_hrs)
data_para = pd.read_csv(data_path_para)

# Start a Dash app
app = dash.Dash(__name__)
server = app.server

# App layout
app.layout = html.Div([
    dcc.Dropdown(
        id='data-source-selector',
        options=[
            {'label': 'Beta HRS', 'value': 'source1'},
            {'label': 'Beta Parallel', 'value': 'source2'}
        ],
        value='source1'
    ),
    dcc.Dropdown(
        id='molecule-selector',
        # Options will be set in the callback
        value=None
    ),
    dcc.Graph(id='3d-scatter-plot')
])


# Callback to update graph based on dropdown
# Callback to update molecule options based on selected data source
@app.callback(
    Output('molecule-selector', 'options'),
    Input('data-source-selector', 'value')
)
def set_molecule_options(selected_source):
    if selected_source == 'source1':
        df = data_hrs
    else:
        df = data_para
    return [{'label': molecule, 'value': molecule} for molecule in df['molecule'].unique()]


# Callback to update the plot based on selected data source and molecule
@app.callback(
    Output('3d-scatter-plot', 'figure'),
    [Input('data-source-selector', 'value'),
     Input('molecule-selector', 'value')]
)
def update_plot(selected_source, selected_molecule):
    if selected_source == 'source1':
        df = data_hrs
    else:
        df = data_para

    filtered_data = df[df['molecule'] == selected_molecule] if selected_molecule else df
    filtered_data = filtered_data.query('b>c')
    # and limit to b<=4
    filtered_data = filtered_data.query('b<=4')


    # Create a figure based on the selected data
    fig = px.scatter_3d(filtered_data, x='b', y='c', z='Beta', color='basis')
    fig.update_layout(width=1000, height=1000)  # Adjust these values as needed

    return fig


# Run the app
if __name__ == '__main__':
    app.run_server(debug=True)
