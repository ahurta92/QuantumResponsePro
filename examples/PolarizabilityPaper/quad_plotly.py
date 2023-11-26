from pathlib import Path
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

import dash
from dash import dcc, html, Input, Output
from quantumresponsepro.BasisMRADataAssembler import make_detailed_df

october_absolute_kain_path = Path('/mnt/data/madness_data/october_beta_absolute_kain/')
november_database_path = Path('/mnt/data/madness_data/beta_paper_database/')
paper_path = Path('/home/adrianhurtado/projects/writing/mra-tdhf-polarizability/Figures_v2')
paper_path = Path('/home/adrianhurtado/projects/writing/thesis2023/PresentationMaterials')

data_path = Path('/home/adrianhurtado/projects/writing/thesis2023/PresentationMaterials/beta_parallel.svg'
                 '.csv')

data_path = Path('/home/adrianhurtado/projects/writing/thesis2023/PresentationMaterials/beta_hrs'
                 '.csv')

# read in data
data = pd.read_csv(data_path)

#data=data.query('b>c')


# Start a Dash app
app = dash.Dash(__name__)

# App layout
app.layout = html.Div([
    dcc.Dropdown(
        id='dataset-dropdown',
        options=[{'label': i, 'value': i} for i in data['molecule'].unique()],
        value=data['molecule'].unique()[0]
    ),
    dcc.Graph(id='3d-scatter-plot')
])



# Callback to update graph based on dropdown
@app.callback(
    Output('3d-scatter-plot', 'figure'),
    Input('dataset-dropdown', 'value')
)
# i need to update the color bar to have more colors

def update_graph(selected_dataset):
    filtered_data = data[data['molecule'] == selected_dataset]
    fig = px.scatter_3d(filtered_data, x='b', y='c', z='Beta', color='basis', opacity=0.7,
                        )


    scale=1.0
    fig.update_layout(width=scale*1000, height=scale*800)  # Set the width and height as desired
    return fig

# Run the app
if __name__ == '__main__':
    app.run_server(debug=True)

