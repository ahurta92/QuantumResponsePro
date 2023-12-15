import json
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import shutil
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import silhouette_score
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler
from pathlib import Path


sns.set_context('paper', )
sns.set_style('darkgrid')

paper_path = Path('/home/adrianhurtado/projects/writing/thesis2023/Figure_v2/')
def mol_valence_matrix(database):
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


# Assuming your data is a 2D numpy array where each row is a vector for a molecule
# data = np.array([...])

DZ = ['aug-cc-pVDZ', 'aug-cc-pCVDZ', 'd-aug-cc-pVDZ', 'd-aug-cc-pCVDZ', 'aug-cc-pVTZ',
      'd-aug-cc-pVTZ']
single = ['aug-cc-pVTZ', 'aug-cc-pCVTZ', 'aug-cc-pVQZ', 'aug-cc-pCVQZ']
double = ['d-aug-cc-pVTZ', 'd-aug-cc-pCVTZ', 'd-aug-cc-pVQZ', 'd-aug-cc-pCVQZ']

def plot_iso_valence_cluster_convergence(data, iso_type, valence, omega, sharey=False):
    data = data.query('omega.isin(@omega) & valence.isin(@valence)')
    facet_kws = {"sharey": sharey, 'despine': True, 'sharex': True, }
    mapping = {i: i + 1 for i in range(len(data.cluster.unique()))}
    data['cluster'] = data['cluster'].map(mapping)
    data['cluster'] = data['cluster'].astype('category')

    with sns.axes_style("darkgrid"):
        f = sns.relplot(data=data,
                        x='valence',
                        kind='line',
                        aspect=0.5,
                        style='Type',
                        hue='Type',
                        col='cluster',
                        facet_kws=facet_kws,
                        palette='colorblind',
                        markers=True,
                        legend=True,
                        errorbar=('ci', 95),
                        y=iso_type
                        )
        f.map(plt.axhline, y=0, color='k', dashes=(1, 1), zorder=0, alpha=0.8, ls='--')
        f.set_titles('{col_name}')
        f.set_xlabels('Valence [n]')
        f.set_ylabels('Percent Error')
        sns.move_legend(f, loc='upper center', ncol=4, title=None,
                        frameon=False)
        # f.axes[0].legend(shadow=True, fancybox=True,fontsize=12)

        f.set(yscale='symlog')

    plt.tight_layout()
    return f








#
# g = analyzer.plot_iso_valence_convergence(mol, 'gamma', ['D', 'T', 'Q', '5'], omega, sharey=True)
# g.fig.show()
# # g = analyzer.plot_iso_valence_convergence('NH3', 'gamma', ['D', 'T', 'Q', '5'], omega)
# # g.fig.show()

#########g = analyzer.plot_alpha_eigen('NH3', 'gamma', ['D', 'T', 'Q', '5'], omega)
#########g.fig.show()
