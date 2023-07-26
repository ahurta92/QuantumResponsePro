import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

from quantumresponsepro import DaltonRunner, BasisMRADataCollection, BasisMRADataAnalyzer
from quantumresponsepro.BasisMRADataAssembler import make_detailed_df
from pathlib import Path

mol = "H2"

xc = "hf"
op = 'dipole'

database_path = Path("/mnt/data/madness_data/post_watoc/alrich_test_set")

# create a DaltonRunner object and set run_dalton to True
dr = DaltonRunner(database_path, True)

all_basis = []
zetas = ['D', 'T', 'Q', '5', '6']

for zeta in zetas:
    basis = ["cc-pV{}Z".format(zeta), "aug-cc-pV{}Z".format(zeta), 'aug-cc-pCV{}Z'.format(zeta),
             'd-aug-cc-pV{}Z'.format(zeta),
             'd-aug-cc-pCV{}Z'.format(zeta), 't-aug-cc-pV{}Z'.format(zeta), ]
    all_basis = all_basis + basis
# run Dalton for the given molecule, basis, and xc
triple = ['t-aug-cc-pVDZ', 't-aug-cc-pVTZ', 't-aug-cc-pVQZ', 't-aug-cc-pCVDZ',
          't-aug-cc-pCVTZ',
          't-aug-cc-pCVQZ', 't-aug-cc-pV5Z', 't-aug-cc-pV6Z']

results = {}

for b in ['def2-QZVPPD']:
    try:
        results[b] = dr.get_polar_json(mol, xc, op, b)
    except (FileNotFoundError, KeyError):
        print(f"File not found for {b}")
        continue

print(results)

# read the results from the database
database = BasisMRADataCollection(database_path, new=True)
analyzer = BasisMRADataAnalyzer(database, .02, )

zeta = 'D'
zetas = []

for zeta in zetas:

    basis = ["aug-cc-pV{}Z".format(zeta), 'aug-cc-pCV{}Z'.format(zeta),
             'd-aug-cc-pV{}Z'.format(zeta),
             'd-aug-cc-pCV{}Z'.format(zeta)]
    data = analyzer.get_basis_energy_data(basis)
    a_data, g_data = analyzer.get_basis_iso_data(basis)

    mols = ['H2', "Ne", "H2O", "SO2"]
    for mol in mols:
        mol_data = data.query("molecule==@mol")
        alpha_data = a_data.query("molecule==@mol")
        datas = [mol_data, alpha_data]
        # create a 1 x 2 subplot

        fig, axs = plt.subplots(nrows=1, ncols=2, )
        fig.suptitle('Basis set convergence for {}'.format(mol))
        titles = ["Ground state energy", "Isotropic polarizability"]
        for i, ax in enumerate(axs):
            pdata = datas[i]
            print(pdata)
            pdata.T.plot(ax=ax, )
            ax.set_title(titles[i])

            ax.axhline(y=0, color='k', linestyle='--') if i == 1 else None
            ax.set_xlabel('Basis set', )
            axs[i].tick_params(axis='x', rotation=60)
            ax.set_ylabel('Absolute Error (a.u.)' if i == 0 else 'Percent Error')
        fig.tight_layout()
        plt.show()

linthresh = 1e-2
mols = ['H2', "Ne", "H2O", "SO2", "NaLi", "HF"]
mols = []
for mol in mols:
    # create a 1 x 2 subplot
    fig, axs = plt.subplots(nrows=1, ncols=2, )
    fig.suptitle('Basis set convergence for {}'.format(mol))
    analyzer.plot_energy_convergence_v2(mol, ['D', 'T', 'Q', '5'], ax=axs[0])
    analyzer.plot_iso_valence_convergence_v2(mol, 'alpha', ['D', 'T', 'Q', '5'],
                                             omega=[0], ax=axs[1])

    # g.fig.suptitle('Basis set convergence for {}'.format(mol))
    fig.tight_layout()
    plt.show()

zetas = ['D', 'T', 'Q']
for mol in mols:
    fig, axes = plt.subplots(nrows=len(zetas), ncols=2, figsize=(10, 10), sharex=True, sharey=False)
    fig.suptitle('Basis set convergence for {}'.format(mol))

    for i, zeta in enumerate(zetas):
        axs = axes[i]
        # create a 1 x 2 subplot
        analyzer.plot_energy_convergence_v2(mol, [zeta], ax=axs[0], reverse=True)
        analyzer.plot_iso_valence_convergence_v2(mol, 'alpha', [zeta],
                                                 omega=[0], ax=axs[1], reverse=True)

        if i == 0:
            axs[0].set_title('Ground state energy')
            axs[1].set_title('Isotropic polarizability')
        else:
            axs[0].set_title('')
            axs[1].set_title('')

        # g.fig.suptitle('Basis set convergence for {}'.format(mol))
    fig.tight_layout()
    plt.show()

zetas = ['D', 'T', 'Q', '5', '6']
all_basis = []

for zeta in zetas:
    basis = ["aug-cc-pV{}Z".format(zeta), 'aug-cc-pCV{}Z'.format(zeta),
             'd-aug-cc-pV{}Z'.format(zeta),
             'd-aug-cc-pCV{}Z'.format(zeta)]
    all_basis = all_basis + basis
# read detailed_all data from feather file if it exist else create it
try:
    # read the detailed_all data from feather file the working directory+ detailed_all_error.feather
    detailed_all_error = pd.read_feather(database_path + 'detailed_all_error.feather')
except:
    detailed_all = pd.DataFrame()
    all_error = pd.DataFrame()
    for freq in range(0, 8):
        print(freq)

        e_data = analyzer.get_basis_energy_data(all_basis)
        a_data, g_data = analyzer.get_basis_iso_data(all_basis, freq)
        # now do the same for the other molecules
        for mol in e_data.index:
            mol_energy_error = e_data.query("molecule==@mol").T
            mol_energy_error['molecule'] = mol
            mol_energy_error['omega'] = freq
            mol_energy_error.rename(columns={mol: 'energy'}, inplace=True)
            # reorder mol_energy_error columns
            mol_energy_error = mol_energy_error[['molecule', 'omega', 'energy']]
            mol_alpha_error = a_data.query("molecule==@mol").T
            mol_gamma_error = g_data.query("molecule==@mol").T
            mol_alpha_error.rename(columns={mol: 'alpha'}, inplace=True)
            mol_gamma_error.rename(columns={mol: 'gamma'}, inplace=True)
            mol_error = pd.concat([mol_energy_error, mol_alpha_error, mol_gamma_error], axis=1)
            all_error = pd.concat([all_error, mol_error], axis=0)

    detailed_all_error = make_detailed_df(
        all_error.reset_index().rename(columns={'index': 'basis'}))

    # save detailed_all_error to feather
    detailed_all_error.to_feather(database_path.joinpath('detailed_all_error.feather'))

sns.set_theme(style="white", color_codes=True, font_scale=2.0)
#

plt.show()

# limit to valence  = D T Q and omega =0
detailed_all_error = detailed_all_error.query("valence in ['D', 'T', 'Q'] ")
detailed_all_error.valence = detailed_all_error.valence.cat.remove_unused_categories()

# sharex=False, sharey=False
facet_kwd = dict(sharex='row', sharey='row', margin_titles=True, despine=False, legend_out=True, )
# Plot miles per gallon against horsepower with other semantics
g = sns.relplot(row='valence', x="alpha", y="energy", col='Type',
                alpha=.8, palette="colorblind", hue='cluster',
                height=6, data=detailed_all_error, facet_kws=facet_kwd, kind='scatter',
                )
# set y lable to energy
g.set_ylabels('Energy error (a.u)')
g.set_xlabels('Percent Error')

# g = sns.FacetGrid(data=detailed_all_error, col='valence', row='Type', margin_titles=True,
#                  sharex='row', sharey='row', despine=False, legend_out=True,
#                  palette='colorblind', height=6)
#
# g.map_dataframe(sns.stripplot, x="energy", y='alpha', hue="cluster", dodge=True,
#                jitter=True, legend="full",
#                )
#
g.set_titles(row_template="{row_name}", col_template="{col_name}")
# move the row titles to the right and the column titles to the top


# create vertical lines for alpha reference values and zero

for ax in g.axes.flat:
    ax.axhline(0, ls='--', color='k')
    ax.axvline(0, ls='--', color='k')
    # y scale is symlog
    # ax.set_yscale('symlog', linthresh=1e-2)
    # x scale is log
    # ax.set_xscale('log')
plt.show()
