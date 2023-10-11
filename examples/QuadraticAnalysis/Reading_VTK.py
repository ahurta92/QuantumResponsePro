import vtk
from pathlib import Path
from quantumresponsepro import MadnessReader
from quantumresponsepro import MadnessResponse
from vtkmodules.util import numpy_support

# Read the .vts file
reader = vtk.vtkXMLStructuredGridReader()
rho0_path = '/mnt/data/madness_data/vtk_plots/rho0.vts'

mol_path = Path('/mnt/data/madness_data/august_no_symmetry/hf/CH3SH/')

vtk_mol_path = mol_path.joinpath('vtk')


# write a function which returns the folder path of the vtk files for each frequency

def get_density_paths(freq_vtk_path, index):
    freq_keys = list(freq_vtk_path.keys())
    density_paths = {}
    density_path = freq_vtk_path[freq_keys[index]].joinpath('density.vts')
    density_names = {}
    density_names[0] = 'ground'
    density_names[1] = 'response_0'
    density_names[2] = 'response_1'
    density_names[3] = 'response_2'
    return density_path, density_names


def read_vtk_data(path, name, min):
    reader = vtk.vtkXMLStructuredGridReader()
    reader.SetFileName(path)
    reader.Update()
    structured_grid = reader.GetOutput()
    dims = structured_grid.GetDimensions()
    bounds = structured_grid.GetBounds()
    print("bounds: ", bounds)
    print(structured_grid)
    structured_grid.GetPointData().GetArray(name)
    x, y, z = np.meshgrid(
        np.linspace(bounds[0], bounds[1], dims[0]),
        np.linspace(bounds[2], bounds[3], dims[1]),
        np.linspace(bounds[4], bounds[5], dims[2])
    )

    # Assuming scalar data is stored in a field named "YourField"
    vtk_data_array = structured_grid.GetPointData().GetArray(name)
    print("vtk_data_array: ", vtk_data_array)
    numpy_data = numpy_support.vtk_to_numpy(vtk_data_array)

    # Reshape the data
    reshaped_data = numpy_data.reshape(dims[2], dims[1], dims[0])
    values = np.abs(reshaped_data.flatten())
    values[np.abs(values) < min] = np.nan

    iso_data = go.Isosurface(
        x=x.flatten(),
        y=y.flatten(),
        z=z.flatten(),
        value=values,
        isomin=values.min(),
        isomax=values.max(),
        colorscale='Picnic',
        surface_count=15,  # number of isosurfaces, 2 by default: only min and max
        caps=dict(x_show=False, y_show=False),
        opacity=0.6,
    )
    return iso_data


def ball_and_stick(xc, op, database_path, molecule, xshift=1.0, yshift=0.0, zshift=0.0):
    mra_mol = MadnessResponse(molecule, xc, op, database_path)
    molecule_dict = mra_mol.ground_info['molecule']
    geometry = molecule_dict['geometry']
    symbols = molecule_dict['symbols']
    x = []
    y = []
    z = []
    for atom in geometry:
        x.append(atom[0] * 1.0 + xshift)
        y.append(atom[1] * 1.0 + yshift)
        z.append(atom[2] * 1.0 + zshift)

    # create a dictionary of colors for each atom
    colors = {'O': 'red', 'H': 'white', 'C': 'black', 'N': 'blue', 'Cl': 'green',
              'Na': 'purple',
              'F': 'orange', 'S': 'yellow', 'P': 'pink', 'I': 'brown', 'Br': 'cyan',
              'B': 'grey',
              'Ne': 'magenta', 'He': 'grey', 'Ar': 'grey', 'Kr': 'grey', 'Xe': 'grey',
              'Li': 'grey', 'Mg': 'grey', 'Al': 'grey', 'Si': 'grey', 'K': 'grey', 'Ca': 'grey',
              'Ti': 'grey', 'Be': 'grey', 'Fe': 'grey', 'Cu': 'grey', 'Zn': 'grey',
              'Ag': 'grey', }
    colors = [colors[s] for s in symbols]
    sizes = [25 for s in symbols]
    # Create scatter plot for atoms
    scatter = go.Scatter3d(x=x, y=y, z=z, mode='markers+text', text=symbols,
                           marker=dict(color=colors, size=sizes, opacity=1.0, ))

    return scatter


# Create the isosurface plot

mol = 'H2O'
xc = 'hf'
op = 'dipole'
database = Path('/mnt/data/madness_data/august_no_symmetry')
database = Path('/mnt/data/madness_data/development')

mad_reader = MadnessReader(database)
freq_vtk_path = mad_reader.get_frequency_vtk_path(mol, xc, op)
freq_keys = list(freq_vtk_path.keys())

density_path, names = get_density_paths(freq_vtk_path, 0)
index = 1
data = []

import dash
from dash import dcc
from dash import html
from dash.dependencies import Input, Output
import plotly.graph_objects as go
import numpy as np

# Initial figure
fig = go.Figure()

molecule = ball_and_stick(xc, op, database, mol)


def get_vtk_datas(path, names, min):
    reader = vtk.vtkXMLStructuredGridReader()
    reader.SetFileName(path)
    reader.Update()
    structured_grid = reader.GetOutput()
    dims = structured_grid.GetDimensions()
    bounds = structured_grid.GetBounds()
    print("bounds: ", bounds)

    x, y, z = np.meshgrid(
        np.linspace(bounds[0], bounds[1], dims[0]),
        np.linspace(bounds[2], bounds[3], dims[1]),
        np.linspace(bounds[4], bounds[5], dims[2])
    )
    # Assuming scalar data is stored in a field named "YourField"
    values = []
    for name in names:
        vtk_data_array = structured_grid.GetPointData().GetArray(name)
        numpy_data = numpy_support.vtk_to_numpy(vtk_data_array)
        # Reshape the data
        reshaped_data = numpy_data.reshape(dims[2], dims[1], dims[0])
        values_i = reshaped_data.flatten()
        values_i[np.abs(values_i) < min] = np.nan
        values.append(values_i)

    return x, y, z, values


X, Y, Z, densities = get_vtk_datas(density_path, names, .01)

# Initialize Dash app
app = dash.Dash(__name__)

# Dash app layout
app.layout = html.Div([
    dcc.Graph(id='iso-plot', figure=fig),
    html.Button('Ground', id='btn-ground'),
    html.Button('Density X', id='btn-x'),
    html.Button('Density Y', id='btn-y'),
    html.Button('Density Z', id='btn-z')
])


# Callback to update isosurface
@app.callback(
    Output('iso-plot', 'figure'),
    [Input('btn-ground', 'n_clicks'),
     Input('btn-x', 'n_clicks'),
     Input('btn-y', 'n_clicks'),
     Input('btn-z', 'n_clicks')]
)
def toggle_density(n1, n2, n3, n4):
    changed_id = [p['prop_id'] for p in dash.callback_context.triggered][0]
    print('changed_id: ', changed_id)
    density = np.zeros_like(X)

    if 'btn-ground' in changed_id:
        density = densities[0]
        print('density: ', density)
    elif 'btn-x' in changed_id:
        density = densities[1]
    elif 'btn-y' in changed_id:
        density = densities[2]
    elif 'btn-z' in changed_id:
        density = densities[3]

    iso_data = go.Isosurface(
        x=X.flatten(),
        y=Y.flatten(),
        z=Z.flatten(),
        value=density,
        isomin=density.min(),
        isomax=density.max(),
        colorscale='Picnic',
        surface_count=15,  # number of isosurfaces, 2 by default: only min and max
        caps=dict(x_show=False, y_show=False),
        opacity=0.3,
    )
    fig = go.Figure(data=[molecule, iso_data])

    return fig


if __name__ == '__main__':
    app.run_server(debug=True)
