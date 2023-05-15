import shutil

from quantumresponsepro import BasisMRADataCollection
from quantumresponsepro import BasisMRADataAnalyzer
import seaborn as sns

from pathlib import Path
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.pyplot as plt
import pandas as pd
import json
from quantumresponsepro.BasisMRADataAssembler import partition_molecule_list
from quantumresponsepro import Tabler
from sklearn.cluster import AgglomerativeClustering
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import StandardScaler
import numpy as np


def set_ax_inset(g: sns.FacetGrid, mol, loc='upper right', yb=.01, iso_type='alpha',
                 omega=[0, 4, 8],
                 vlevel=['T', 'Q', '5'],
                 width='80%', height='50%'):
    b = 0
    ax = g.ax
    sns.move_legend(g, 'lower right')
    ylims = [3, .8, .18]
    ylims = dict(zip(['D', 'T', 'Q'], ylims))
    print(ylims)
    ax_inset = inset_axes(ax, width=2, height=2, bbox_to_anchor=(1.55, 1.0),
                          bbox_transform=ax.transAxes, borderpad=0,
                          )

    # Plot the same data in the inset_axes but zoom in on the specified range
    zoom_range = (0.75, 3.50)
    # ax_inset.set_xlim(*zoom_range)
    analyzer.plot_iso_valence_convergence_inset(ax_inset, mol=mol, iso_type='alpha', valence=vlevel,
                                                omega=omega, )
    # Optionally, draw a rectangle in the main plot to indicate the zoomed area
    rect = plt.Rectangle((zoom_range[0], -yb), zoom_range[1] - zoom_range[0], 2 * yb, lw=1,
                         edgecolor='red',
                         linestyle='--', facecolor='none')
    ax.add_patch(rect)
    ax_inset.set_ylim(-yb, yb)
    b += 1


august = Path('/mnt/data/madness_data/post_watoc/august')
paper_path = Path('/home/adrianhurtado/projects/writing/mra-tdhf-polarizability/Figures_v2')
thesis_path = Path('/home/adrianhurtado/projects/writing/thesis2023/Figures_v2')
paper_path = thesis_path
database = BasisMRADataCollection(august)
analyzer = BasisMRADataAnalyzer(database, .05, font_scale=1.5)
tabler = Tabler(database)


def mol_valence_matrix(mol):
    df = pd.DataFrame()
    for V in ['D', 'T', 'Q']:
        mol_valence = database.detailed_iso_diff.query(
            'molecule==@mol & valence==@V & omega==8').copy()
        mol_valence['abs_alpha'] = mol_valence.alpha.abs()
        mol_valence = mol_valence.sort_values('abs_alpha').drop('abs_alpha', axis=1)
        vcol = mol_valence['augmentation'].astype('str') + mol_valence['polarization'].astype(
            'str').map({'CV': '+core', 'V': ''})
        df[V] = vcol.reset_index(drop=True)
    return df


def get_convergence_data(mol_list):
    my_dict = {}
    for mol in mol_list:
        mol_key = mol_valence_matrix(mol).loc[0][1:].to_json()
        if mol_key in my_dict:
            my_dict[mol_key].append(mol)
        else:
            my_dict[mol_key] = [mol]

    df = pd.DataFrame()
    b = 0
    for key, val in my_dict.items():
        mol_series = pd.Series(json.loads(key))
        mol_series['molecules'] = val
        # mol_series=pd.DataFrame(key)
        df[b] = mol_series
        b += 1
    return df


# Assuming your data is a 2D numpy array where each row is a vector for a molecule
# data = np.array([...])

DZ = ['aug-cc-pVDZ', 'aug-cc-pCVDZ', 'd-aug-cc-pVDZ', 'd-aug-cc-pCVDZ', 'aug-cc-pVTZ',
      'd-aug-cc-pVTZ']
single = ['aug-cc-pVTZ', 'aug-cc-pCVTZ', 'aug-cc-pVQZ', 'aug-cc-pCVQZ']
double = ['d-aug-cc-pVTZ', 'd-aug-cc-pCVTZ', 'd-aug-cc-pVQZ', 'd-aug-cc-pCVQZ']


def cluster_basis_data(data, n_clusters=8):
    scaler = StandardScaler()
    X = data.to_numpy()
    data_scaled = scaler.fit_transform(X)

    # Create an AgglomerativeClustering instance with n_clusters
    agglo = AgglomerativeClustering(n_clusters=n_clusters)

    # Fit the model to your data and get the cluster assignments in one step
    labels = agglo.fit_predict(data_scaled)
    data['cluster'] = labels

    return data


first_row_path = paper_path.joinpath('molecules/first')
second_row_path = paper_path.joinpath('molecules/second')
fluorine_path = paper_path.joinpath('molecules/fluorine')
molecules_path = paper_path.joinpath('molecules')

from quantumresponsepro.BasisMRADataAssembler import partition_molecule_list

first, second, fluorine = partition_molecule_list(database.molecules)

omega = [0, 4, 8]

width = '30%'
height = '20%'

i = 2
mapping = {
    'aug': 1,
    'aug+core': 2,
    'd-aug': 3,
    'd-aug+core': 4,
}
df = pd.DataFrame()
for mol in database.molecules:
    df[mol] = mol_valence_matrix(mol).stack().map(mapping)
df.reset_index(drop=True, inplace=True)

# data = tabler.get_basis_data()[0][DZ].abs().dropna()
# data_first = data.query('molecule in @second')
# data_first = data_first.query('molecule != "Ne"')
X = cluster_basis_data(df.T, n_clusters=8)
print(X)

for cluster in X['cluster'].unique():
    mol_list = X.query('cluster==@cluster').index
    cluster_path = molecules_path.joinpath(f'cluster_{cluster}')
    if not cluster_path.exists():
        cluster_path.mkdir()
    else:
        shutil.rmtree(cluster_path)
        cluster_path.mkdir()

    for mol in mol_list:
        g = analyzer.plot_iso_valence_convergence_v2(mol, 'alpha', ['D', 'T', 'Q', '5'], omega)
        set_ax_inset(g, mol, loc='lower right', yb=.10, iso_type='alpha', omega=omega, width=width,
                     height=height)
        g.fig.savefig(cluster_path.joinpath(f'{mol}_alpha_converge.svg'), dpi=300)

# g = analyzer.plot_iso_valence_convergence_v2(mol, 'alpha', ['D', 'T', 'Q', '5'], omega, sharey=True)
# g.fig.show()
#
# g = analyzer.plot_iso_valence_convergence(mol, 'gamma', ['D', 'T', 'Q', '5'], omega, sharey=True)
# g.fig.show()
# # g = analyzer.plot_iso_valence_convergence('NH3', 'gamma', ['D', 'T', 'Q', '5'], omega)
# # g.fig.show()

#########g = analyzer.plot_alpha_eigen('NH3', 'gamma', ['D', 'T', 'Q', '5'], omega)
#########g.fig.show()
