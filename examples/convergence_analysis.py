import json
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from pathlib import Path
from quantumresponsepro import BasisMRADataAnalyzer
from quantumresponsepro import BasisMRADataCollection
from quantumresponsepro.BasisMRADataAssembler import partition_molecule_list


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
analyzer = BasisMRADataAnalyzer(database, .05, )


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


first, second, fluorine = partition_molecule_list(database.molecules)
conv_df = get_convergence_data(first)


omega = [0, 4, 8]

width = '30%'
height = '20%'

i = 2

print(conv_df.iloc[:-1, i])
for mol in conv_df.loc['molecules'][i]:
    g = analyzer.plot_iso_valence_convergence_v2(mol, 'alpha', ['D', 'T', 'Q', '5'], omega)
    set_ax_inset(g, mol, loc='lower right', yb=.10, iso_type='alpha', omega=omega, width=width,
                 height=height)
    g.fig.show()
    # g.fig.savefig(first_row_path.joinpath(f'{mol}_alpha_converge.svg'), dpi=300)

# g = analyzer.plot_iso_valence_convergence_v2(mol, 'alpha', ['D', 'T', 'Q', '5'], omega, sharey=True)
# g.fig.show()
#
# g = analyzer.plot_iso_valence_convergence(mol, 'gamma', ['D', 'T', 'Q', '5'], omega, sharey=True)
# g.fig.show()
# # g = analyzer.plot_iso_valence_convergence('NH3', 'gamma', ['D', 'T', 'Q', '5'], omega)
# # g.fig.show()

#########g = analyzer.plot_alpha_eigen('NH3', 'gamma', ['D', 'T', 'Q', '5'], omega)
#########g.fig.show()
