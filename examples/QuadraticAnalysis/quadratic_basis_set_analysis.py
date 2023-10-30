import math

import matplotlib as mpl
import matplotlib.pyplot as plt
# Create some example data
import numpy as np
import pandas as pd
import seaborn as sns
from pathlib import Path
from quantumresponsepro import BasisMRAData
from quantumresponsepro.BasisMRADataAssembler import make_detailed_df

fd_path = Path('/mnt/data/madness_data/fd_compare3')
thesis_path = Path('/home/adrianhurtado/projects/writing/thesis2023/Figures_v2')
paper_path = Path('/home/adrianhurtado/projects/writing/mra-tdhf-polarizability/Figures_v2')
tromso_path = Path('/home/adrianhurtado/projects/writing/tromso_poster/figures')

paper_path = thesis_path

database_path = Path('/mnt/data/madness_data/august_no_symmetry')

database = BasisMRAData(database_path, new=True)
print(database.all_polar_data)

# show the entire dataframe by setting the max rows to None
pd.set_option('display.max_columns', None)

# get a detailed view of the data
# now I need to take the difference of each basis set value with the MADNESS value
# I'll compute the percent error of each value
# drop Afreq Bfreq and Cfreq greater than 0

alpha_data = database.all_polar_data.copy()
alpha_data = alpha_data.query('omega==0')


# drop MRA data less than 1e-4
def tensor_to_numpy(j):
    array = np.empty(j["size"])
    array[:] = j["vals"]
    return np.reshape(array, tuple(j["dims"]))


class QuadDataAnalyzer:
    def __init__(self, db_path, molecule, type=None, ):
        self.quad_data = None
        self.energy_diff = None
        self.alpha_diff = None
        self.database = db_path
        mra_data = BasisMRAData(db_path, new=False)
        self.quad_data = mra_data.all_quad_data.copy().query('Afreq==0 &Bfreq==0 &Cfreq==0 & '
                                                             'Beta.abs()>1e-3')
        print(self.quad_data)
        self.quad_data.sort_values(by='ijk', inplace=True)

        quad_diff = self.compute_molecule_quad_diff(self.quad_data, molecule).reset_index()
        quad_diff = make_detailed_df(quad_diff)

        alpha_data = mra_data.all_polar_data.copy().query('omega==0')

        alpha_diff = compute_molecule_alpha_diff(alpha_data, molecule).reset_index()
        alpha_diff = make_detailed_df(alpha_diff)

        energy_diff = mra_data.energy_diff.copy().query('molecule==@molecule').reset_index()
        print(energy_diff)
        energy_diff = make_detailed_df(energy_diff)
        if type is not None:
            print(type, "is not `None`")
            energy_diff = energy_diff.query('Type.isin(@type)').copy()
            alpha_diff = alpha_diff.query('Type.isin(@type)').copy()
            quad_diff = quad_diff.query('Type.isin(@type)').copy()
            # drop unused categories
            energy_diff['Type'] = energy_diff['Type'].cat.remove_unused_categories()
            alpha_diff['Type'] = alpha_diff['Type'].cat.remove_unused_categories()
            quad_diff['Type'] = quad_diff['Type'].cat.remove_unused_categories()
            print(energy_diff)

        self.energy_diff = energy_diff
        self.alpha_diff = alpha_diff
        self.quad_diff = quad_diff

    def computed_molecule_quad_basis_diff(self, quad_data, molecule, basis):
        multi_index = ['ijk', 'Afreq', 'Bfreq', 'Cfreq', 'molecule']
        # get the MADNESS value for the molecule
        mad_data = quad_data.query('molecule==@molecule & basis=="MRA"').copy()
        # get the basis set value for the molecule
        basis_data = quad_data.query('molecule==@molecule & basis==@basis').copy().drop_duplicates()
        # from the unique ijk from the basis set data
        print(basis_data)
        print(mad_data)

        basis_data.set_index(multi_index, inplace=True)
        mad_data.set_index(multi_index, inplace=True)

        basis_data['Beta'].values[:] = np.sign(basis_data['Beta'].values * mad_data[
            'Beta'].values) * basis_data['Beta'].values
        basis_data['diff'] = (basis_data['Beta'] - mad_data['Beta']) / mad_data['Beta'] * 100
        # drop the Beta column
        basis_data.drop(columns=['Beta'], inplace=True)
        return basis_data

    def compute_molecule_quad_diff(self, quad_data, molecule):
        data = []
        b_data = quad_data.query('basis!="MRA"').copy().reset_index()
        print(b_data)
        for basis in b_data.basis.unique():
            data.append(self.computed_molecule_quad_basis_diff(quad_data, molecule, basis))

        return pd.concat(data, axis=0, )


def plot_component_diff(data, valence, component, ax, color_pal=None):
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
    print(N)
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


def compute_molecule_alpha_basis_diff(alpha_data, molecule, basis='MRA'):
    alpha_data = alpha_data.query('alpha.abs()>1e-3').copy()
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


def compute_molecule_alpha_diff(alpha_data, molecule):
    data = []
    b_data = alpha_data.query('basis!="MRA"')
    for basis in b_data.basis.unique():
        print(basis)
        data.append(compute_molecule_alpha_basis_diff(alpha_data, molecule, basis))

    return pd.concat(data, axis=0, )


pal = sns.color_palette("seismic_r", 4).as_hex()
p1 = [pal[1], pal[0], pal[2], pal[3]]
pal = sns.color_palette(p1)

pal2 = sns.color_palette("coolwarm_r", 4).as_hex()
p2 = [pal2[1], pal2[0], pal2[2], pal2[3]]
p3 = [pal2[1], pal2[2], ]
light_pal = sns.color_palette(p2)
simple_pal = sns.color_palette(p3)

sns.set_theme('paper', 'darkgrid', palette=light_pal, font='sans-serif', font_scale=1.0)


def basis_set_analysis_plot(database_path, molecule, valence, type=None):
    qda = QuadDataAnalyzer(database_path, molecule, type)

    quad_data = qda.quad_data.query('Beta.abs() > 1e-2').query(
        'Afreq == 0 & Bfreq == 0 & Cfreq == 0')
    alpha_data = qda.alpha_diff.query('omega==0')

    r_plot, axes = plt.subplots(1, 3, figsize=(9, 4), frameon=True, layout='constrained')
    if type is None:
        pal = light_pal
    else:
        pal = simple_pal

    plot_energy_diff(qda.energy_diff, valence, axes[0], color_pal=pal)
    alpha_plot = plot_component_diff(qda.alpha_diff, valence, 'ij', axes[1], color_pal=pal)
    beta_plot = plot_component_diff(qda.quad_diff, valence, 'ijk', axes[2], color_pal=pal)

    axes[0].set_title('Energy')

    axes[1].set_title(r'$\alpha_{ij}(0;0,0)$')
    axes[1].set_yscale('symlog', linthresh=1e-2, linscale=0.25)

    axes[2].set_title(r'$\beta_{ijk}(0;0,0)$')
    # axes[2].set_yscale('symlog', linthresh=1e-2, linscale=0.25)

    return r_plot, axes


mol = 'H2O'
r_plot, axes = basis_set_analysis_plot(database_path, mol, ['D', 'T', 'Q', ], type=None)
r_plot.savefig(paper_path.joinpath("{}.png".format(mol)), dpi=1000)
