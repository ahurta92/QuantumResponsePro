import math

import pandas as pd
import matplotlib as mpl
from quantumresponsepro import BasisMRAData, DaltonRunner

from quantumresponsepro import BasisMRADataAnalyzer
from quantumresponsepro import MadnessResponse

import seaborn as sns

from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# Create some example data
import numpy as np
from quantumresponsepro.BasisMRADataAssembler import make_detailed_df
from ipywidgets import interact
import plotly.graph_objects as go

fd_path = Path('/mnt/data/madness_data/fd_compare3')
thesis_path = Path('/home/adrianhurtado/projects/writing/thesis2023/Figures_v2')
paper_path = Path('/home/adrianhurtado/projects/writing/mra-tdhf-polarizability/Figures_v2')
tromso_path = Path('/home/adrianhurtado/projects/writing/tromso_poster/figures')

paper_path = thesis_path

database_path = Path('/mnt/data/madness_data/post_watoc/august_high_prec')
compare_path = fd_path.joinpath('high-high')

pal = sns.color_palette("seismic_r", 4).as_hex()
p1 = [pal[1], pal[0], pal[2], pal[3]]
pal = sns.color_palette(p1)

pal2 = sns.color_palette("coolwarm_r", 4).as_hex()
p2 = [pal2[1], pal2[0], pal2[2], pal2[3]]
p3 = [pal2[1], pal2[2], ]
light_pal = sns.color_palette(p2)
simple_pal = sns.color_palette(p3)


# drop MRA data less than 1e-4
def tensor_to_numpy(j):
    array = np.empty(j["size"])
    array[:] = j["vals"]
    return np.reshape(array, tuple(j["dims"]))


class QuadDataAnalyzer:
    def __init__(self, db_path, molecule, type=None, new=False):
        self.molecule = molecule
        self.quad_data = None
        self.energy_diff = None
        self.alpha_diff = None
        self.database = db_path
        mra_data = BasisMRAData(db_path, new=new)
        # first prune quad data
        self.quad_data = (mra_data.all_quad_data.copy()).query(
            'molecule==@molecule')
        print(self.quad_data)

        self.quad_data.reset_index(inplace=True)

        # also limit the frequencies to frequencies for which MRA data is available

        self.quad_data['Afreq'] = self.quad_data['Afreq'].apply(lambda x: math.trunc(x * 1e3) / 1e3)
        self.quad_data['Bfreq'] = self.quad_data['Bfreq'].apply(lambda x: math.trunc(x * 1e3) / 1e3)
        self.quad_data['Cfreq'] = self.quad_data['Cfreq'].apply(lambda x: math.trunc(x * 1e3) / 1e3)

        # truncate Beta to 5 decimal places
        self.quad_data['Beta'] = self.quad_data['Beta'].apply(lambda x: math.trunc(x * 1e4) / 1e4)

        # change the sign of Afreq for if basis is MRA
        self.quad_data['Afreq'] = self.quad_data.apply(
            lambda x: -x['Afreq'] if x['basis'] == 'MRA' else x['Afreq'], axis=1)

        avail_freqs = self.quad_data.query('basis=="MRA"').Afreq.unique()

        self.quad_data = self.quad_data.query('Afreq.isin(@avail_freqs)')
        self.quad_data.set_index('ijk', inplace=True)
        # I rename the indices so that MRA and dalton indices match
        # self.quad_data.rename(index={'XXY': 'XYX', 'XXZ': 'XZX', 'XYZ': 'XZY', 'YYZ': 'YZY'},
        #                      inplace=True)

        zijk = ['XXZ', 'XZX', 'XYZ', 'XZY', 'XZZ', 'YXZ', 'YZX', 'YYZ', 'YZY', 'YZZ', 'ZXX', 'ZXY',
                'ZXZ',
                'ZYX', 'ZYY', 'ZZX', 'ZZY', 'ZZZ']
        zijk = []
        # zijk = ['XYX', 'XZX', 'YYY', 'YZY', 'YZZ', 'ZZZ', ]

        # if the ijk contains a Z then the sign of the Beta is negative
        self.quad_data.loc[((self.quad_data['basis'] == 'MRA') & (self.quad_data.index.isin(
            zijk))), 'Beta'] = self.quad_data.loc[((self.quad_data['basis'] == 'MRA') & (
            self.quad_data.index.isin(zijk))), 'Beta'] * -1.0

        self.quad_data.reset_index(inplace=True)

        multi_index = ['ijk', 'Afreq', 'Bfreq', 'Cfreq', 'molecule']
        self.quad_data.set_index(multi_index, inplace=True)
        basis_data = self.quad_data.query('basis=="aug-cc-pVQZ"').copy(
        ).drop_duplicates()
        unique_index = basis_data.index.unique()
        # MRA and basis dimension do not mathc after the query
        # mad_data = self.quad_data.query('molecule==@molecule & basis=="MRA" ').copy()
        # mad_data = mad_data.query('index.isin(@unique_index)')
        # self.quad_data = self.quad_data.query('index.isin(@unique_index)')
        # self.quad_data.drop_duplicates(inplace=True)

        # rename madness data to mad_beta
        # mad_data.rename(columns={'Beta': 'mad_beta'}, inplace=True)
        # truncate frequency columns at 5 decimal places
        # basis_data = pd.concat([basis_data, mad_data], axis=1)
        # basis_data['sign'] = np.sign(basis_data['Beta'].values * basis_data['mad_beta'].values)
        # basis_data['Madness'] = basis_data['mad_beta'].values * basis_data['sign'].values

        # print('BASIS DATA:\n', basis_data)

        # set MRA Beta data in quad_data to basis_data
        # if basis  is MRA Beta data = basis_data
        # self.quad_data.loc[self.quad_data['basis'] == 'MRA', 'Beta'] = basis_data['Madness']

        # now compute the differences
        # quad_diff = self.compute_molecule_quad_diff(self.quad_data, molecule).reset_index()
        # quad_diff = make_detailed_df(quad_diff)

        self.alpha_data = mra_data.alpha_eigen.copy()

        alpha_diff = self.compute_molecule_alpha_diff(molecule).reset_index()
        alpha_diff = make_detailed_df(alpha_diff)

        energy_diff = mra_data.energy_diff.copy().query('molecule==@molecule').reset_index()
        energy_diff = make_detailed_df(energy_diff)
        if type is not None:
            print(type, "is not `None`")
            energy_diff = energy_diff.query('Type.isin(@type)').copy()
            alpha_diff = alpha_diff.query('Type.isin(@type)').copy()
            # quad_diff = quad_diff.query('Type.isin(@type)').copy()
            # drop unused categories
            energy_diff['Type'] = energy_diff['Type'].cat.remove_unused_categories()
            alpha_diff['Type'] = alpha_diff['Type'].cat.remove_unused_categories()
            # quad_diff['Type'] = quad_diff['Type'].cat.remove_unused_categories()
            print(energy_diff)

        self.energy_diff = energy_diff
        self.alpha_diff = alpha_diff
        # self.quad_diff = quad_diff

    def computed_molecule_quad_basis_diff(self, quad_data, molecule, basis):
        # get the MADNESS value for the molecule
        mad_data = quad_data.query('molecule==@molecule & basis=="MRA"').copy()
        # get the basis set value for the molecule
        #  negate the sign of Afreq for mad_data
        basis_data = quad_data.query('molecule==@molecule & basis==@basis').copy().drop_duplicates()

        mad_data.rename(columns={'Beta': 'mad_beta'}, inplace=True)
        basis_data = pd.concat([basis_data, mad_data.drop(columns='basis')], axis=1)
        basis_data['sign'] = np.sign(basis_data['Beta'].values * basis_data['mad_beta'].values)
        basis_data['Madness'] = basis_data['mad_beta'].values * basis_data['sign'].values

        basis_data['diff'] = (basis_data['Beta'] - basis_data['Madness']) / basis_data[
            'Madness'] * 100
        # drop the Beta column
        basis_data.drop(columns=['Beta'], inplace=True)
        basis_data.drop(columns=['Madness'], inplace=True)
        basis_data.drop(columns=['sign'], inplace=True)

        return basis_data

    def compute_molecule_quad_diff(self, quad_data, molecule):
        data = []
        b_data = quad_data.query('basis!="MRA"')
        for basis in b_data.basis.unique():
            data.append(self.computed_molecule_quad_basis_diff(quad_data, molecule, basis))

        return pd.concat(data, axis=0, )

    def compute_molecule_alpha_basis_diff(self, molecule, basis='MRA'):
        alpha_data = self.alpha_data.query('alpha.abs()>1e-3').copy()
        multi_index = ['ij', 'omega']
        # get the MADNESS value for the molecule
        mad_data = alpha_data.query('molecule==@molecule & basis=="MRA"').copy()
        basis_data = alpha_data.query('molecule==@molecule & basis==@basis').copy()
        ij_min = basis_data.ij.unique()
        # now remove the ijk values from the madness data that are not in the basis set data
        mad_data = mad_data.query('ij.isin(@ij_min)')

        basis_data.set_index(multi_index, inplace=True)
        mad_data.set_index(multi_index, inplace=True)
        # now compute the percent error in the basis-set data
        basis_data['diff'] = (basis_data['alpha'] - mad_data['alpha']) / mad_data['alpha'] * 100
        # drop the Beta column
        basis_data.drop(columns=['alpha'], inplace=True)
        return basis_data

    def compute_molecule_alpha_diff(self, molecule):
        data = []
        b_data = self.alpha_data.query('basis!="MRA"')
        for basis in b_data.basis.unique():
            data.append(self.compute_molecule_alpha_basis_diff(molecule, basis))

        return pd.concat(data, axis=0, )

    def unit_sphere_representation(self, basis, omega_1, omega_2, radius=1, num_points=1000,
                                   basis_error=False, colormap='Blues', sizeref=0.5):
        if basis_error:
            data = self.quad_diff
            # perhaps change diff to Beta
            data['Beta'] = data['diff']
        else:
            data = self.quad_data
        data = data.query('basis==@basis & Afreq==@omega_1 & Bfreq==@omega_2')
        # Initialize results
        results = []
        # Generate 1000 points on the unit sphere
        points = generate_points_on_sphere(radius, num_points)

        beta_tensor = beta_df_np(data)
        # Compute projection
        for p in points:
            bi = beta_proj(beta_tensor, p)
            results.append(bi)
        results = np.array(results)

        p1 = points.copy()
        p2 = points.copy() + results

        x = p1[:, 0]
        y = p1[:, 1]
        z = p1[:, 2]
        u = p2[:, 0]
        v = p2[:, 1]
        w = p2[:, 2]

        quiver = go.Cone(x=x, y=y, z=z, u=u, v=v, w=w, colorscale=colormap, sizemode="scaled",
                         sizeref=sizeref)
        # return the quiver
        return quiver

    def unit_sphere_representation_basis(self, basis_sets, omega_1, omega_2, radius=1,
                                         num_points=1000,
                                         colormap='Blues', sizeref=0.5, shift=1.0):
        data = self.quad_diff
        data['Beta'] = data['diff']

        x = np.zeros((len(basis_sets) * num_points))
        y = np.zeros((len(basis_sets) * num_points))
        z = np.zeros((len(basis_sets) * num_points))
        u = np.zeros((len(basis_sets) * num_points))
        v = np.zeros((len(basis_sets) * num_points))
        w = np.zeros((len(basis_sets) * num_points))

        # for each basis set generate the xyzuvw data
        for i, basis in enumerate(basis_sets):
            bdata = data.query('basis==@basis & Afreq==@omega_1 & Bfreq==@omega_2')
            # Initialize results
            results = []
            # Generate 1000 points on the unit sphere
            points = generate_points_on_sphere(radius, num_points)

            beta_tensor = beta_df_np(bdata)
            # Compute projection
            for p in points:
                bi = beta_proj(beta_tensor, p)
                results.append(bi)
            results = np.array(results)
            p1 = points.copy()
            p2 = points.copy() + results
            x[i * num_points:(i + 1) * num_points] = p1[:, 0] + shift * (i - len(basis_sets) / 2)
            y[i * num_points:(i + 1) * num_points] = p1[:, 1]
            z[i * num_points:(i + 1) * num_points] = p1[:, 2]

            u[i * num_points:(i + 1) * num_points] = p2[:, 0] + shift * (i - len(basis_sets) / 2)
            v[i * num_points:(i + 1) * num_points] = p2[:, 1]
            w[i * num_points:(i + 1) * num_points] = p2[:, 2]
            # append the data to the quiver

        quiver = go.Cone(x=x, y=y, z=z, u=u, v=v, w=w, colorscale=colormap, sizemode="scaled",
                         sizeref=sizeref, cmin=0)
        # return the quiver
        return quiver

    def unit_sphere_representation_MRA(self, basis_sets, omega_1, omega_2, radius=1,
                                       num_points=1000,
                                       colormap='Blues', sizeref=0.5, shift=1.0, sizemode='scaled'):
        data = self.quad_data

        x = np.zeros((len(basis_sets) * num_points))
        y = np.zeros((len(basis_sets) * num_points))
        z = np.zeros((len(basis_sets) * num_points))
        u = np.zeros((len(basis_sets) * num_points))
        v = np.zeros((len(basis_sets) * num_points))
        w = np.zeros((len(basis_sets) * num_points))

        text = []

        for i, basis in enumerate(basis_sets):
            for j in range(num_points):
                text.append(basis)

        #minimal_index_ijk = ['XXX', 'YYY', 'ZZZ', 'XYZ', 'XXY', 'YXY', 'ZXX', 'ZXZ', 'YYZ', 'YZZ']

        # for each basis set generate the xyzuvw data
        for i, basis in enumerate(basis_sets):
            bdata = data.query('basis==@basis & Afreq==@omega_1 & Bfreq==@omega_2')
            print(bdata)
            bdata.reset_index(inplace=True)

            bdata.set_index('ijk', inplace=True)
            bdata.sort_index(inplace=True)
            #bdata = bdata.loc[minimal_index_ijk]
            # Initialize results
            results = []
            # Generate 1000 points on the unit sphere
            points = generate_points_on_sphere(radius, num_points)

            beta_tensor = beta_df_np(bdata)
            # Compute projection
            for p in points:
                bi = beta_proj(beta_tensor, p)
                results.append(bi)
            results = np.array(results)
            p1 = points.copy()
            p2 = points.copy() + results
            x[i * num_points:(i + 1) * num_points] = p1[:, 0] + shift * (i - len(basis_sets) / 2)
            y[i * num_points:(i + 1) * num_points] = p1[:, 1]
            z[i * num_points:(i + 1) * num_points] = p1[:, 2]

            u[i * num_points:(i + 1) * num_points] = p2[:, 0] + shift * (i - len(basis_sets) / 2)
            v[i * num_points:(i + 1) * num_points] = p2[:, 1]
            w[i * num_points:(i + 1) * num_points] = p2[:, 2]

            # also generate text basis set name data to add to quiver
            # append the data to the quiver

        quiver = go.Cone(x=x, y=y, z=z, u=u, v=v, w=w, colorscale=colormap, sizemode=sizemode,
                         sizeref=sizeref, text=text)
        # return the quiver
        return quiver

    def basis_set_analysis_plot(self, molecule, valence, type=None):
        quad_data = self.quad_diff.query('Beta.abs() > 1e-2').query(
            'Afreq == 0 & Bfreq == 0 & Cfreq == 0')
        alpha_data = self.alpha_diff.query('omega==0')

        r_plot, axes = plt.subplots(1, 3, figsize=(9, 4), frameon=True, layout='constrained')
        if type is None:
            pal = light_pal
        else:
            pal = simple_pal

        plot_energy_diff(qda.energy_diff, valence, axes[0], color_pal=pal)
        alpha_plot = plot_component_diff(qda.alpha_diff.query("omega==0"), valence, 'ij', axes[1],
                                         color_pal=pal)
        beta_plot = plot_component_diff(qda.quad_diff.query("Afreq==0"), valence, 'ijk', axes[2],
                                        color_pal=pal)

        axes[0].set_title('Energy')

        axes[1].set_title(r'$\alpha_{ij}(0;0,0)$')
        axes[1].set_yscale('symlog', linthresh=1e-2, linscale=0.25)

        axes[2].set_title(r'$\beta_{ijk}(0;0,0)$')
        # axes[2].set_yscale('symlog', linthresh=1e-2, linscale=0.25)

        return r_plot, axes


def generate_points_on_sphere(radius, num_points):
    points = np.zeros((num_points, 3))
    phi = math.pi * (3. - math.sqrt(5.))  # Golden angle

    for i in range(num_points):
        y = radius - (i / float(num_points - 1)) * 2 * radius  # y goes from 1 to -1
        ri = math.sqrt(radius ** 2 - y * y)  # radius at y

        theta = phi * i  # Golden angle increment

        x = math.cos(theta) * ri
        z = math.sin(theta) * ri

        points[i, :] = [x, y, z]

    return points


def beta_df_np(beta_df):
    beta_df = beta_df.copy()
    xyz_to_012 = {'X': 0, 'Y': 1, 'Z': 2}
    beta_tensor = np.zeros((3, 3, 3))
    beta_df.reset_index(inplace=True)
    beta_df.set_index('ijk', inplace=True)
    beta_df.sort_index(inplace=True)
    print('beta_df_np: ', beta_df)
    # for each row in beta_zero
    for i, row in beta_df.iterrows():
        # get the indices of the row
        indices = [xyz_to_012[x] for x in i]
        # set the value of the tensor at the indices to the value of the row
        # by kleinman symmetry
        beta_tensor[indices[0], indices[1], indices[2]] = row['Beta']
        beta_tensor[indices[0], indices[2], indices[1]] = row['Beta']
        beta_tensor[indices[2], indices[0], indices[1]] = row['Beta']

    return beta_tensor


def beta_proj(beta, E):
    return np.tensordot(beta, np.outer(E, E), axes=([0, 1], [0, 1])).T


# write a function which takes in a molecule and returns the geometry and symbols
# from the MadnessResponse Class

def get_molecule_geometry(molecule, database_path):
    mra_mol = MadnessResponse(molecule, 'hf', 'dipole', database_path)
    molecule_dict = mra_mol.ground_info['molecule']
    geometry = molecule_dict['geometry']
    symbols = molecule_dict['symbols']
    return geometry, symbols


# we will use plotly for this
def ball_and_stick(molecule, database_path):
    mra_mol = MadnessResponse(molecule, 'hf', 'dipole', database_path)
    molecule_dict = mra_mol.ground_info['molecule']
    geometry = molecule_dict['geometry']
    symbols = molecule_dict['symbols']

    x = []
    y = []
    z = []
    for atom in geometry:
        x.append(atom[0])
        y.append(atom[1])
        z.append(atom[2])

    # create a dictionary of colors for each atom
    colors = {'O': 'red', 'H': 'white', 'C': 'black', 'N': 'blue', 'Cl': 'green', 'Na': 'purple',
              'F': 'orange', 'S': 'yellow', 'P': 'pink', 'I': 'brown', 'Br': 'cyan', 'B': 'grey',
              'Ne': 'magenta', 'He': 'grey', 'Ar': 'grey', 'Kr': 'grey', 'Xe': 'grey',
              'Li': 'grey', 'Mg': 'grey', 'Al': 'grey', 'Si': 'grey', 'K': 'grey', 'Ca': 'grey',
              'Ti': 'grey', 'Be': 'grey', 'Fe': 'grey', 'Cu': 'grey', 'Zn': 'grey', 'Ag': 'grey', }
    colors = [colors[s] for s in symbols]
    sizes = [50 for s in symbols]
    # Create scatter plot for atoms
    scatter = go.Scatter3d(x=x, y=y, z=z, mode='markers', text=symbols,
                           marker=dict(color=colors, size=sizes, opacity=1.0, ))

    return scatter


def plot_component_diff(data, valence, component, ax, color_pal=None):
    data = data.copy()
    data['valence'] = pd.Categorical(data['valence'], valence)
    # make the polarization order ['V','CV']
    # make the polarization order ['V','CV']
    try:
        data['polarization'] = pd.Categorical(data['polarization'], ['V', 'CV'])
    except ValueError as v:
        print(v)
        print(data)
        pass
    # font scale *2
    g = sns.lineplot(data=data, x='valence', y='diff', hue='Type',
                     style=component,
                     markers=True, ax=ax, markersize=10, legend='full', dashes=False,
                     palette=color_pal)

    handles, labels = g.get_legend_handles_labels()
    N = math.ceil(len(handles) / 2)
    filtered_handles = handles[N:]
    filtered_labels = labels[N:]
    g.legend(filtered_handles, filtered_labels)

    ax.get_xaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
    ax.get_yaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
    ax.minorticks_on()
    ax.grid(which='major', color='w', linewidth=1.0)
    ax.grid(which='minor', color='w', linewidth=0.5)
    ax.set_xlabel('Valence')
    ax.set_ylabel('Percent Error')
    ax.axhline(y=0, color='k')

    # for the line plots, only keep the "style" legend not the hue legend
    handles, labels = ax.get_legend_handles_labels()

    # add the legends for the components only

    return g


def plot_energy_diff(data, valence, ax, color_pal=None):
    data['valence'] = pd.Categorical(data['valence'], valence)
    # make the polarization order ['V','CV']
    try:
        data['polarization'] = pd.Categorical(data['polarization'], ['V', 'CV'])
    except ValueError as v:
        print(v)
        print(data)
        pass

    # show minor ticks and grid lines  using with
    g = sns.lineplot(data=data, x='valence', y='error', hue='Type',
                     markers=True, ax=ax, palette=color_pal)

    ax.get_xaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
    ax.get_yaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
    ax.minorticks_on()
    ax.grid(which='major', color='w', linewidth=1.0)
    ax.grid(which='minor', color='w', linewidth=0.5)

    ax.set_xlabel('Valence')
    ax.set_ylabel('Error (a.u.)')
    ax.axhline(y=0, color='k')

    return g


compare_path = Path('/mnt/data/madness_data/post_watoc/august_high_prec')

molecule = 'LiBH4'
qda = QuadDataAnalyzer(database_path, molecule, new=False)
# r_plot, axes = qda.basis_set_analysis_plot(compare_path, 'CH3SH', ['D', 'T', 'Q', ], type=None)
#########r_plot.savefig(paper_path.joinpath("nacl_plot.png"), dpi=1000)
# plt.show()

pd.set_option('display.max_columns', None)
# CH3NH2 CH3OH CH2BH C6H6 CH3SH
database_path = Path('/mnt/data/madness_data/fd_compare3/high-high')
database_path = Path('/mnt/data/madness_data/post_watoc/august_high_prec')

scatter = ball_and_stick(molecule, database_path)
quiver = qda.unit_sphere_representation_MRA(['aug-cc-pVQZ', 'MRA', ],
                                            0.0, 0.0,
                                            sizeref=3.5, shift=2.5, sizemode='scaled')

# Create the figure and plot
fig = go.Figure(data=[scatter, quiver, ])
fig.update_layout(scene=dict(aspectmode='cube', xaxis=dict(range=[-5, 5]),
                             yaxis=dict(range=[-5, 5]),
                             zaxis=dict(range=[-5, 5])))
fig.show()

# get the dipole moments
