import json
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import shutil
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from pathlib import Path
from quantumresponsepro import BasisMRADataAnalyzer
from quantumresponsepro import BasisMRADataCollection
from quantumresponsepro import Tabler
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import silhouette_score
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler

august = Path('/mnt/data/madness_data/post_watoc/august')
paper_path = Path('/home/adrianhurtado/projects/writing/mra-tdhf-polarizability/Figures_v2')
thesis_path = Path('/home/adrianhurtado/projects/writing/thesis2023/Figures_v2')
tromso_poster_path = Path('/home/adrianhurtado/projects/writing/tromso_poster/figures')
tromso_poster_path_supple = Path(
    '/home/adrianhurtado/projects/writing/tromso_poster/supplementary/figures')
paper_path = paper_path
database = BasisMRADataCollection(august)
analyzer = BasisMRADataAnalyzer(database, .02)
tabler = Tabler(database)

linthresh = .01
subticks = [2, 3, 4, 5, 6, 7, 8, 9]
linscale = 0.10

sns.set_context('paper')
sns.set_theme(style="whitegrid", font_scale=1.0, rc={'ytick.left': True, 'xtick.bottom': False, })


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


def set_face_color_cluster(g: sns.FacetGrid):
    pal2 = sns.color_palette("coolwarm_r", 4).as_hex()
    p2 = [pal2[1], pal2[0], pal2[2], pal2[3]]
    light_pal = sns.color_palette(p2)
    b = 0
    axes = g.axes
    for i, axi in enumerate(axes):
        back_color = light_pal[i]
        back_color = back_color + (0.2,)
        print(i, back_color)
        for j, axij in enumerate(axi):
            axij.set_facecolor(back_color)


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
        if (btype == 'd-aug-cc-pCVnZ') and (vlevel == 'Q' or vlevel == 'T'):
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
    data_matrix = data.to_numpy()
    threshold = analyzer.mra_ref

    data_matrix = symlog(data_matrix, threshold)

    data_matrix = scaler.fit_transform(data_matrix)

    # pca = PCA(n_components=4)
    # principal_components = pca.fit_transform(data_matrix)

    # print(data_scaled)
    # data_matrix = principal_components

    n_clusters_range = range(1, 16)
    bic_scores = []
    convergence_type = 'full'
    tol = 1e-3

    for n_clusters in n_clusters_range:
        # Create and fit a Gaussian Mixture Model
        gmm = GaussianMixture(n_components=n_clusters, covariance_type=convergence_type,
                              random_state=42,
                              tol=tol, n_init=10, )
        gmm.fit(data_matrix)

        # Compute BIC for the current clustering
        bic_scores.append(gmm.bic(data_matrix))

    # Find the number of clusters that gives the minimum BIC
    optimal_clusters = n_clusters_range[np.argmin(bic_scores)]
    # optimal_clusters = 4
    print(f"Optimal number of clusters: {optimal_clusters}")

    plt.figure()  # Adjust the size of the plot as needed
    plt.plot(n_clusters_range, bic_scores, marker='o')
    plt.xlabel('Number of clusters')
    plt.ylabel('BIC score')
    plt.title('BIC score per number of clusters')
    plt.tight_layout()
    plt.savefig(tromso_poster_path.joinpath('BIC.svg'))
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
                          n_init=40,
                          tol=tol, )
    gmm.fit(data_matrix)

    # Fit the model to your data and get the cluster assignments in one step
    labels = gmm.fit_predict(data_matrix)
    score = silhouette_score(data_matrix, labels, metric='euclidean')
    data['cluster'] = labels
    print("silhouette score", score)

    average_vectors = []
    std_vectors = []
    for i in range(optimal_clusters):
        cluster_points = data_matrix[labels == i]
        average_vector = np.mean(cluster_points, axis=0)
        std_vector = np.std(cluster_points, axis=0)
        average_vectors.append(average_vector)
        std_vectors.append(std_vector)

    avg_df = pd.DataFrame(average_vectors, )
    avg_df['mean'] = avg_df.mean(axis=1)
    avg_df = avg_df.sort_values('mean', ascending=False)
    avg_df.drop('mean', axis=1, inplace=True)

    std_df = pd.DataFrame(std_vectors, )

    sorted_index = avg_df.index
    std_df = std_df.reindex(sorted_index)
    std_df.sort_index(inplace=True)

    cluster_map = {sorted_index[i]: i for i in range(len(sorted_index))}
    avg_df = avg_df.reset_index(drop=True)
    # avg_df = avg_df.apply(lambda x: inv_symlog(x, threshold))
    # avg_df = pd.DataFrame(scaler.inverse_transform(avg_df), columns=avg_df.columns,
    #                      index=avg_df.index)
    print(cluster_map)
    data['cluster'] = data['cluster'].map(cluster_map)
    # avg_df = np.exp(avg_df) + xmin

    # Generate a line plot
    fig = plt.figure()  # Adjust the size of the plot as needed
    for i, row in avg_df.iterrows():
        plt.plot(row, label=f'Cluster {i}', marker='o')
    # plt.yscale('log')
    plt.axhline(y=0, color='black', linestyle='--', linewidth=1)  # Add a horizontal line at 0
    plt.xlabel('Basis set')  # Label for the x-axis
    plt.ylabel('Average basis set error')  # Label for the y-axis
    plt.title('Average basis set errors for each cluster')  # Title of the plot
    plt.legend()  # Show the legend
    plt.xticks(rotation=60)  # Rotate the x-axis labels for better visibility if needed
    plt.tight_layout()  # Adjust the layout for better visibility if needed
    plt.show()
    fig.savefig(tromso_poster_path.joinpath('GMM_averages.svg'))

    return data, avg_df


def plot_iso_valence_cluster_frequecy(data, iso_type, valence, omega, sharey=False):
    data = data.query('omega.isin(@omega) & valence.isin(@valence)')
    facet_kws = {"sharey": sharey, 'despine': True, }
    data['valence'].cat.remove_unused_categories()

    f = sns.relplot(data=data,
                    x='omega',
                    kind='scatter',
                    size=4.25,
                    style='Type',
                    hue='Type',
                    col='valence',
                    row='cluster',
                    facet_kws=facet_kws,
                    palette='colorblind',
                    y=iso_type
                    )
    f.map(plt.axhline, y=0, color='k', dashes=(2, 1), zorder=0, alpha=0.8)
    f.map(plt.axhline, y=1, color='red', dashes=(2, 1), zorder=0, alpha=0.8)
    f.map(plt.axhline, y=-1, color='red', dashes=(2, 1), zorder=0, alpha=0.8)
    f.set_titles('{col_name}')
    f.set_xlabels('Valence [n]')
    f.set_ylabels('Percent Error')

    f.tight_layout()
    for ax in f.axes.flatten():  # If you have multiple axes due to faceting
        ax.grid(visible=True, which="both", axis='y', alpha=0.5)
        ax.minorticks_on()
        ax.tick_params(axis='y', which="both", top="off", left="on", right="on", bottom="off", )
    # 'g' now is a FacetGrid object, and we can access the underlying axes
    return f


def plot_iso_valence_cluster_convergence(data, iso_type, valence, omega, sharey=False):
    data = data.query('omega.isin(@omega) & valence.isin(@valence)').copy()
    facet_kws = {"sharey": sharey, 'despine': True, 'sharex': True,
                 "margin_titles": True, }

    mapping = {i: i + 1 for i in range(len(data.cluster.unique()))}
    data['cluster'] = data['cluster'].map(mapping)
    data['cluster'] = data['cluster'].astype('category')

    with sns.axes_style("darkgrid") and sns.plotting_context("paper", font_scale=2.3):
        f = sns.relplot(data=data,
                        x='valence',
                        y=iso_type,
                        col='cluster',
                        height=5,
                        aspect=0.5,
                        row='mol_system',
                        kind='line',
                        hue='Type',
                        facet_kws=facet_kws,
                        palette='colorblind',
                        markers=True,
                        legend='brief',
                        errorbar=('ci', 95),
                        )
        f.map(plt.axhline, y=0, color='k', dashes=(1, 1), zorder=0, alpha=0.8, ls='--')
        f.set_titles(row_template="{row_name}", col_template="Cluster {col_name}", template="{"
                                                                                            "row_name}\
                                                                                         | {col_name}")

        f.set_xlabels('Valence [n]')
        f.set_ylabels('Percent Error')
        sns.move_legend(
            f, "upper center",
            ncol=4, title=None, frameon=False)
        f.figure.subplots_adjust(top=0.93, right=0.97, left=0.08, bottom=0.08, wspace=0.05)

    return f


def plot_iso_valence_cluster_convergence_2(data, iso_type, valence, omega, sharey=False):
    data = data.query('omega.isin(@omega) & valence.isin(@valence)').copy()
    facet_kws = {"sharey": sharey, 'despine': True, 'sharex': True, }
    mapping = {i: i + 1 for i in range(len(data.cluster.unique()))}
    data['cluster'] = data['cluster'].map(mapping)
    data['cluster'] = data['cluster'].astype('category')

    with sns.axes_style("darkgrid"):
        f = sns.relplot(data=data,
                        x='valence',
                        kind='line',
                        aspect=0.50,
                        col='Type',
                        style='cluster',
                        hue='cluster',
                        row='mol_system',
                        facet_kws=facet_kws,
                        palette='colorblind',
                        markers=True,
                        legend='brief',
                        errorbar=('ci', 95),
                        y=iso_type
                        )
        f.map(plt.axhline, y=0, color='k', dashes=(1, 1), zorder=0, alpha=0.8, ls='--')
        f.set_titles('Cluster {col_name}')
        f.set_xlabels('Valence [n]')
        f.set_ylabels('Percent Error')
        # sns.move_legend(
        #    f, "upper center",
        #    ncol=4, title=None, frameon=False, fontsize=12)
        f.figure.subplots_adjust(top=0.87, right=0.99, left=0.08, bottom=0.15, wspace=0.05)

    return f


first_row_path = paper_path.joinpath('molecules/first')
second_row_path = paper_path.joinpath('molecules/second')
fluorine_path = paper_path.joinpath('molecules/fluorine')
molecules_path = paper_path.joinpath('molecules')

omega = [0, 4, 8]
aspect = 0.9
g = analyzer.freq_iso_plot_v2(['D', 'T', 'Q'], 'alpha', 'all', True, omegas=omega,
                              pal='colorblind',
                              )
for ax in g.axes_dict.values():
    ax.axhline(-analyzer.mra_ref, color='green', linestyle='--')
    ax.axhline(analyzer.mra_ref, color='green', linestyle='--')
    ax.set_yscale('symlog', linthresh=linthresh, base=10, linscale=linscale,
                  subs=subticks)
    ax.xaxis.grid(True, "minor", linewidth=0.25)
    ax.yaxis.grid(True, "minor", linewidth=0.25)
# e_fig.despine(left=True, bottom=True)`
set_face_color(g.axes_dict)
# set_ax_inset(g, loc='lower right', iso_type='alpha', width='50%', height='50%')
g.fig.savefig(paper_path.joinpath('alpha_freq_DTQ.svg'), dpi=300)

from quantumresponsepro.BasisMRADataAssembler import partition_molecule_list

first, second, fluorine = partition_molecule_list(database.molecules)

omega = [0, 8]

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

cluster_path = molecules_path.joinpath(f'alpha_clusters')
if not cluster_path.exists():
    cluster_path.mkdir()
else:
    shutil.rmtree(cluster_path)
    cluster_path.mkdir()

grouped = X.groupby('cluster')

cluster_mol_path = cluster_path.joinpath(f'cluster_molecules.txt')
with open(cluster_mol_path, 'w') as f:
    for cluster_id, group in grouped:
        f.write(r'\textbf{' + f'Cluster {cluster_id + 1} ' + f'({len(group.index)})' + ':} ')
        for molecule in group.index:
            f.write("\ce{" + f"{molecule}" + "}, ")
        f.write('\n\n')

for cluster_id, group in grouped:
    cluster_mol_path = cluster_path.joinpath(f'cluster_{cluster_id + 1}_molecules.txt')

    with open(cluster_mol_path, 'w') as f:
        f.write(r'\textbf{' + f'Cluster {cluster_id + 1} ' + f'({len(group.index)})' + ':} ')
        for molecule in group.index:
            f.write("\ce{" + f"{molecule}" + "}, ")

iso_diff = database.detailed_iso_diff.copy()
iso_diff['cluster'] = iso_diff['molecule'].map(X['cluster'])
iso_diff['cluster'] = iso_diff['cluster'].astype('category')
g = plot_iso_valence_cluster_convergence(iso_diff, 'alpha', ['D', 'T', 'Q', '5'], omega,
                                         sharey='row')


# Define a custom formatting function for the axis ticks
def custom_format(value, pos):
    # Apply your desired formatting logic here
    formatted_value = "{:.2f}".format(value)  # Format the value as needed
    return formatted_value


for ax in g.axes_dict.values():
    ax.axhline(-analyzer.mra_ref, color='green', linestyle='--')
    ax.axhline(analyzer.mra_ref, color='green', linestyle='--')
    ax.set_yscale('symlog', linthresh=linthresh, base=10, linscale=linscale,
                  subs=subticks)
    ax.xaxis.grid(True, "minor", linewidth=0.25)
    ax.yaxis.grid(True, "minor", linewidth=0.25)
# e_fig.despine(left=True, bottom=True)
# Add 'H2' text to the first subplot
axs = g.axes  # Get the axes object from the FacetGrid
print(axs)
axs[0, 4].text(0.05, 0.95, 'O$_3$', ha='left', va='top', transform=axs[0, 4].transAxes,
               bbox=dict(facecolor='red', alpha=0.5, edgecolor='black'), fontsize=20)
axs[0, 5].text(0.05, 0.95, 'Be', ha='left', va='top', transform=axs[0, 5].transAxes,
               bbox=dict(facecolor='red', alpha=0.5, edgecolor='black'), fontsize=20)

axs[1, 4].text(0.05, 0.95, 'BF', ha='left', va='top', transform=axs[1, 4].transAxes,
               bbox=dict(facecolor='red', alpha=0.5, edgecolor='black'), fontsize=20)

axs[2, 2].text(0.05, 0.95, 'SiH$_3$Cl,  PH$_3$', ha='left', va='top', transform=axs[2, 2].transAxes,
               bbox=dict(facecolor='red', alpha=0.5, edgecolor='black'), fontsize=20)

# Add 'N2' text to the second subplot
g.fig.savefig(cluster_path.joinpath('alpha_convergence.svg'))
g.fig.show()

g = analyzer.freq_iso_plot_cluster(iso_diff, ['D', 'T', 'Q'], 'alpha',
                                   omegas=[0, 4, 8],
                                   sharey=False, aspect=1)
for ax in g.axes_dict.values():
    ax.axhline(-analyzer.mra_ref, color='green', linestyle='--')
    ax.axhline(analyzer.mra_ref, color='green', linestyle='--')
    ax.set_yscale('symlog', linthresh=linthresh, base=10, linscale=linscale,
                  subs=subticks)
    ax.xaxis.grid(True, "minor", linewidth=0.25)
    ax.yaxis.grid(True, "minor", linewidth=0.25)
# for ax in g.axes_dict.values():
#    ax.set_yscale('symlog', linthresh=linthresh, base=10, linscale=linscale,
#                  subs=subticks)
set_face_color_cluster(g)
# freq_inset(g, iso_diff, loc='lower right', ylims=[3, .8, .18], width='30%', height='40%',
# omega=[8])
g.fig.savefig(cluster_path.joinpath('alpha_frequency_convergence.svg'))
g.fig.show()

plt.figure()

with sns.axes_style('darkgrid') and sns.plotting_context('paper', font_scale=1.0):
    p = plt.figure(figsize=(4.5, 3))
    ax = p.add_subplot(111)
    data = iso_diff.copy()
    mapping = {i: i + 1 for i in range(len(data.cluster.unique()))}
    data['cluster'] = data['cluster'].map(mapping)
    data['cluster'] = data['cluster'].astype('category')
    sns.histplot(data.query('omega==0 & basis =="aug-cc-pVQZ"'), x="cluster",
                 hue="mol_system",
                 multiple="stack",
                 ax=ax,
                 shrink=.8)
    plt.xticks(range(1, len(data.cluster.unique())))
    sns.move_legend(
        ax, "lower center",
        bbox_to_anchor=(.5, 1), ncol=3, title=None, frameon=False,
    )
    plt.xlabel('')
    # p.subplots_adjust(bottom=0.13, left=0.15, right=0.95, top=0.90)
    plt.show()
    p.savefig(cluster_path.joinpath('cluster_molecule_count.svg'))

if False:

    for cluster in X['cluster'].unique():
        mol_list = X.query('cluster==@cluster').index
        cluster_path_i = cluster_path.joinpath(f'cluster_{cluster + 1}')
        if not cluster_path_i.exists():
            cluster_path_i.mkdir()
        else:
            shutil.rmtree(cluster_path)
            cluster_path_i.mkdir()

        for mol in mol_list:
            g = analyzer.plot_iso_valence_convergence_v2(mol, 'alpha', ['D', 'T', 'Q', '5', '6'],
                                                         omega)
            g.fig.savefig(cluster_path_i.joinpath(f'{mol}_alpha_converge.svg'), dpi=300)

# g.fig.show()
#
# g = analyzer.plot_iso_valence_convergence(mol, 'gamma', ['D', 'T', 'Q', '5'], omega, sharey=True)
# g.fig.show()
# # g = analyzer.plot_iso_valence_convergence('NH3', 'gamma', ['D', 'T', 'Q', '5'], omega)
# # g.fig.show()

#########g = analyzer.plot_alpha_eigen('NH3', 'gamma', ['D', 'T', 'Q', '5'], omega)
#########g.fig.show()
