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
from sklearn.cluster import DBSCAN
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture
from sklearn.mixture import BayesianGaussianMixture
from sklearn.metrics import silhouette_score

from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import StandardScaler
import numpy as np


def set_face_color(ax_dict):
    pal2 = sns.color_palette("coolwarm_r", 4).as_hex()
    p2 = [pal2[1], pal2[0], pal2[2], pal2[3]]
    light_pal = sns.color_palette(p2)
    b = 0
    for i, ax in ax_dict.items():
        back_color = light_pal[b % 4]
        back_color = back_color + (0.3,)
        ax.set_facecolor(back_color)

        b += 1


def freq_inset(g: sns.FacetGrid, data, loc='upper right', ylims=[3, .8, .18],
               iso_type='alpha',
               omega=[0, 4, 8],
               width='80%', height='50%'):
    b = 0
    ax_dict = g.axes_dict
    ylims = [3, .8, .18]
    ylims = dict(zip(['D', 'T', 'Q'], ylims))
    print(ylims)

    for i, ax in ax_dict.items():
        vlevel = i[0]
        btype = i[1]
        print(vlevel)
        yb = ylims[vlevel]
        if i[1] == 'd-aug-cc-pCVnZ' or i[1] == 'd-aug-cc-pVnZ' or (
                (vlevel == 'Q') and (i[1] == 'aug-cc-pCVnZ' or btype == 'aug-cc-pVnZ')):
            # Create the inset_axes instance
            ax_inset = inset_axes(ax, width=width, height=height, loc=loc)

            # Plot the same data in the inset_axes but zoom in on the specified range
            zoom_range = (1.50, 2.50)
            # ax_inset.set_xlim(*zoom_range)
            analyzer.cluster_iso_plot_ax(data, ax_inset, v_level=vlevel, b_type=btype,
                                         iso_type=iso_type,
                                         omegas=omega,
                                         pal='colorblind', )
            # ax_inset.set_title('Zoomed Inset')

            # Optionally, draw a rectangle in the main plot to indicate the zoomed area
            rect = plt.Rectangle((zoom_range[0], -yb), zoom_range[1] - zoom_range[0], 2 * yb, lw=1,
                                 edgecolor='red',
                                 linestyle='--', facecolor='none')
            ax.add_patch(rect)
            ax_inset.set_ylim(-yb, yb)
            b += 1
        else:
            ax.axhline(yb, color='k', linestyle='--', linewidth=1)
            ax.axhline(yb, color='k', linestyle='--', linewidth=1)


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


def inv_symlog(y, linthresh):
    """Inverse of symmetric log transformation."""

    return np.sign(y) * linthresh * (np.exp(np.abs(y)) - 1)


def symlog(x, linthresh):
    """Symmetric log transformation."""

    return np.sign(x) * np.log1p(np.abs(x / linthresh))


def cluster_gaussian_mixture_with_bic(data):
    """
    Bayesian Information Criterioon (BIC) for model selction amoung a finite set of models
    Based on likelhood function and introduces a penalty term for complexity to avoid overfitting
    In this context used to select the number of clusters
    """
    scaler = StandardScaler()
    X = data.to_numpy()
    threshold = 5e-2

    X = symlog(X, threshold)

    X = scaler.fit_transform(X)
    # print(data_scaled)
    print(X)

    n_clusters_range = range(1, 10)
    bic_scores = []
    convergence_type = 'full'
    tol = 1e-3

    for n_clusters in n_clusters_range:
        # Create and fit a Gaussian Mixture Model
        gmm = GaussianMixture(n_components=n_clusters, covariance_type=convergence_type,
                              random_state=42,
                              tol=tol, n_init=10, )
        gmm.fit(X)

        # Compute BIC for the current clustering
        bic_scores.append(gmm.bic(X))

    # Find the number of clusters that gives the minimum BIC
    optimal_clusters = n_clusters_range[np.argmin(bic_scores)]
    # optimal_clusters = 8
    print(f"Optimal number of clusters: {optimal_clusters}")

    plt.plot(n_clusters_range, bic_scores, marker='o')
    plt.xlabel('Number of clusters')
    plt.ylabel('BIC score')
    plt.title('BIC score per number of clusters')
    plt.show()

    # Create an AgglomerativeClustering instance with n_clusters
    # gmm = BayesianGaussianMixture(n_components=optimal_clusters,
    #                               random_state=None,
    #                               covariance_type='full',
    #                               tol=1e-6, )
    gmm = GaussianMixture(n_components=optimal_clusters,
                          covariance_type=convergence_type,
                          random_state=42,
                          warm_start=True,
                          reg_covar=1e-4,
                          n_init=10,
                          tol=tol, )
    gmm.fit(X)

    # Fit the model to your data and get the cluster assignments in one step
    labels = gmm.fit_predict(X)
    score = silhouette_score(X, labels, metric='euclidean')
    data['cluster'] = labels
    print("silhouette score", score)

    average_vectors = []
    std_vectors = []
    for i in range(optimal_clusters):
        cluster_points = X[labels == i]
        average_vector = np.mean(cluster_points, axis=0)
        std_vector = np.std(cluster_points, axis=0)
        average_vectors.append(average_vector)
        std_vectors.append(std_vector)

    avg_df = pd.DataFrame(average_vectors, columns=data.columns[:-1])
    avg_df['mean'] = avg_df.mean(axis=1)
    avg_df = avg_df.sort_values('mean', ascending=False)
    avg_df.drop('mean', axis=1, inplace=True)

    std_df = pd.DataFrame(std_vectors, columns=data.columns[:-1])

    sorted_index = avg_df.index
    std_df = std_df.reindex(sorted_index)
    std_df.sort_index(inplace=True)

    cluster_map = {sorted_index[i]: i for i in range(len(sorted_index))}
    avg_df = avg_df.reset_index(drop=True)
    # avg_df = avg_df.apply(lambda x: inv_symlog(x, threshold))
    avg_df = pd.DataFrame(scaler.inverse_transform(avg_df), columns=avg_df.columns,
                          index=avg_df.index)
    print(cluster_map)
    data['cluster'] = data['cluster'].map(cluster_map)
    # avg_df = np.exp(avg_df) + xmin

    # Generate a line plot
    plt.figure(figsize=(10, 6))  # Adjust the size of the plot as needed
    for i, row in avg_df.iterrows():
        plt.plot(row, label=f'Cluster {i}', marker='o')
    # plt.yscale('log')
    plt.axhline(y=0, color='black', linestyle='--', linewidth=1)  # Add a horizontal line at 0
    plt.xlabel('Basis set')  # Label for the x-axis
    plt.ylabel('Average basis set error')  # Label for the y-axis
    plt.title('Average basis set errors for each cluster')  # Title of the plot
    plt.legend()  # Show the legend
    plt.xticks(rotation=90)  # Rotate the x-axis labels for better visibility if needed
    plt.tight_layout()  # Adjust the layout for better visibility if needed
    plt.show()

    return data, avg_df


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
    'd-aug': -1,
    'aug+core': 1.25,
    'd-aug+core': -1.25,
}
df = pd.DataFrame()
for mol in database.molecules:
    df[mol] = mol_valence_matrix(mol).stack().map(mapping)
df.reset_index(drop=True, inplace=True)

DZ = [
    'aug-cc-pVDZ', 'd-aug-cc-pVDZ',
    'aug-cc-pCVTZ', 'd-aug-cc-pCVTZ',
    'aug-cc-pVQZ', 'd-aug-cc-pVQZ',
    'aug-cc-pCVQZ', 'd-aug-cc-pCVQZ',
]

data = tabler.get_basis_data()[0].dropna()
# data_first = data.query('molecule in @second')
# data_first = data_first.query('molecule != "Ne"')
# X = cluster_basis_data(df.T, n_clusters=4)
# X = cluster_basis_data_DBSCAN(df.T, eps=0.5, min_samples=5)
X, avg_vectors = cluster_gaussian_mixture_with_bic(data)

avg_vectors.plot()

cluster_path = molecules_path.joinpath(f'clusters')
if not cluster_path.exists():
    cluster_path.mkdir()
else:
    shutil.rmtree(cluster_path)
    cluster_path.mkdir()

grouped = X.groupby('cluster')

for cluster_id, group in grouped:
    cluster_mol_path = cluster_path.joinpath(f'cluster_{cluster_id}_molecules.txt')
    with open(cluster_mol_path, 'w') as f:
        for molecule in group.index:
            f.write(f"{molecule}\n")

iso_diff = database.detailed_iso_diff.copy()
iso_diff['cluster'] = iso_diff['molecule'].map(X['cluster'])
g = analyzer.freq_iso_plot_cluster(iso_diff, ['D', 'T', 'Q'], 'alpha', 'all', 'row', omegas=omega,
                                   pal='colorblind',
                                   )
freq_inset(g, iso_diff, omega=[8], loc='lower right', iso_type='alpha', width='50%', height='30%')
set_face_color(g.axes_dict)
g.fig.show()
# Generate a line plot
basis_order = [
    'aug-cc-pVDZ',
    'aug-cc-pCVDZ',
    'd-aug-cc-pVDZ',
    'd-aug-cc-pCVDZ',
    'aug-cc-pVTZ',
    'aug-cc-pCVTZ',
    'd-aug-cc-pVTZ',
    'd-aug-cc-pCVTZ',
    'aug-cc-pVQZ',
    'aug-cc-pCVQZ',
]
rdata = iso_diff.query('basis==@basis_order')
sns.relplot(data=rdata, x=rdata.basis, y='alpha', hue='cluster', kind='line', errorbar='sd',
            palette='colorblind', height=10, aspect=1.5)
# plt.yscale('log')
plt.axhline(y=0, color='black', linestyle='--', linewidth=1)  # Add a horizontal line at 0
plt.xlabel('Basis set')  # Label for the x-axis
plt.ylabel('Average basis set error')  # Label for the y-axis
plt.yscale('symlog')  # Set the y-axis to be symlog
plt.title('Average basis set errors for each cluster')  # Title of the plot
plt.legend()  # Show the legend
plt.xticks(rotation=90)  # Rotate the x-axis labels for better visibility if needed
plt.tight_layout()  # Adjust the layout for better visibility if needed
plt.show()

for cluster in X['cluster'].unique():
    mol_list = X.query('cluster==@cluster').index
    cluster_path_i = cluster_path.joinpath(f'cluster_{cluster}')
    if not cluster_path_i.exists():
        cluster_path_i.mkdir()
    else:
        shutil.rmtree(cluster_path)
        cluster_path_i.mkdir()

    for mol in mol_list:
        g = analyzer.plot_iso_valence_convergence_v2(mol, 'alpha', ['D', 'T', 'Q', '5'], omega)
        set_ax_inset(g, mol, loc='lower right', yb=.10, iso_type='alpha', width=width,
                     height=height)
        g.fig.savefig(cluster_path_i.joinpath(f'{mol}_alpha_converge.svg'), dpi=300)

# g.fig.show()
#
# g = analyzer.plot_iso_valence_convergence(mol, 'gamma', ['D', 'T', 'Q', '5'], omega, sharey=True)
# g.fig.show()
# # g = analyzer.plot_iso_valence_convergence('NH3', 'gamma', ['D', 'T', 'Q', '5'], omega)
# # g.fig.show()

#########g = analyzer.plot_alpha_eigen('NH3', 'gamma', ['D', 'T', 'Q', '5'], omega)
#########g.fig.show()
