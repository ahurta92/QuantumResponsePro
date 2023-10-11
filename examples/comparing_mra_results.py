import seaborn as sns
from matplotlib import pyplot as plt
from pathlib import Path
from quantumresponsepro import MadnessResponse, BasisMRADataAnalyzer
from quantumresponsepro.BasisMRADataCollection import BasisMRADataCollection

mol = "H2"
xc = "hf"
operator = "dipole"

data_dir_high = Path("/mnt/data/madness_data/post_watoc/alrich_test_set")

zetas = ['D', 'T', 'Q', '5', '6']
all_basis = []

for zeta in zetas:
    basis = ["cc-pV{}Z".format(zeta), "aug-cc-pV{}Z".format(zeta), 'aug-cc-pCV{}Z'.format(zeta),
             'd-aug-cc-pV{}Z'.format(zeta),
             'd-aug-cc-pCV{}Z'.format(zeta)]
    all_basis = all_basis + basis + ['']

molecules = ['H2']
op = 'dipole'

data_dir_high = Path("/mnt/data/madness_data/post_watoc/alrich_test_set")
database_high = BasisMRADataCollection(data_dir_high, new=True)

iso_diff = database_high.detailed_iso_diff
facet_kws = {'sharey': True, 'sharex': False}

zeta = ['D', 'T', 'Q', '5', '6']

g = sns.relplot(data=iso_diff.query("valence.isin(@zeta)"), x='valence', y='alpha',
                hue='Type',
                kind='line',
                style='Type',
                facet_kws=facet_kws,
                col='omega',
                col_wrap=3,
                markers=True,
                )

for ax in g.axes.flat:
    ax.axhline(0, ls='--', color='k')
    ax.axhline(.02, ls='--', color='green')
    ax.axhline(-.02, ls='--', color='green')
    # ax.axvline(0, ls='--', color='k')
    ax.set_yscale('symlog', linthresh=1e-2)
plt.show()

mols = ['H2']

data_dir_1 = Path("/mnt/data/madness_data/post_watoc/august")
database2 = BasisMRADataCollection(data_dir_1)

analyzer1 = BasisMRADataAnalyzer(database_high, .02)
analyzer2 = BasisMRADataAnalyzer(database2, .02)
fig, ax = plt.subplots(2, 2, figsize=(10, 10))

axs = ax.flatten()

freqs = [0, 8]

for i, omega in enumerate(freqs):
    g = analyzer1.plot_iso_valence_convergence_v2('H2', 'alpha', ['D', 'T,', 'Q', '5', '6'],
                                                  [omega],
                                                  axs[i])
for i, omega in enumerate(freqs):
    g = analyzer2.plot_iso_valence_convergence_v2('H2', 'alpha', ['D', 'T,', 'Q', '5', '6'],
                                                  [omega],
                                                  axs[i + 2])
plt.show()

g = analyzer1.plot_alpha_component_convergence('H2', ['xx', 'yy', 'zz'], ['D', 'T,', 'Q', '5',
                                                                          '6'], sharey=True)
for ax in g.axes.flat:
    ax.axhline(0, ls='--', color='k')
    ax.axhline(.02, ls='--', color='green')
    ax.axhline(-.02, ls='--', color='green')
    ax.set_yscale('symlog', linthresh=.02, linscale=1)

plt.show()

mols = ['H2', 'H2O', 'NaCl', 'NaCN', 'HF', 'F2', 'SO2', 'HCl']
for mol in mols:
    g = analyzer2.plot_alpha_component_convergence(mol, ['xx', 'yy', 'zz'], ['D', 'T,', 'Q', '5',
                                                                             '6'], sharey=True)
    for ax in g.axes.flat:
        ax.axhline(0, ls='--', color='k')
        ax.axhline(.02, ls='--', color='green')
        ax.axhline(-.02, ls='--', color='green')
        ax.set_yscale('symlog', linthresh=.02, linscale=1)

    plt.show()

zetas = ['Q', '5', '6']
all_basis = []

for zeta in zetas:
    basis = ["aug-cc-pV{}Z".format(zeta), 'aug-cc-pCV{}Z'.format(zeta),
             'd-aug-cc-pV{}Z'.format(zeta),
             'd-aug-cc-pCV{}Z'.format(zeta), 't-aug-cc-pV{}Z'.format(zeta), ]
    all_basis = all_basis + basis + ['def2-TZVPPD', 'def2-QZVPPD']

mad_high = MadnessResponse('H2', 'hf', 'dipole', data_dir_high)
comp_xx = mad_high.compare_to_dalton(data_dir_1, ['xx', ],
                                     all_basis)
comp_yy = mad_high.compare_to_dalton(data_dir_1, ['yy', ],
                                     all_basis)
comp_zz = mad_high.compare_to_dalton(data_dir_1, ['zz', ],
                                     all_basis)

print(comp_xx)
print(comp_yy)
print(comp_zz)
comp_ = mad_high.compare_to_dalton(data_dir_1, ['xx', 'yy', 'zz', ],
                                   all_basis)
comp_def2_QZVPPD = mad_high.compare_to_dalton(data_dir_high, ['xx', 'yy', 'zz', ], ['def2-QZVPPD'])

iso_data = database_high.iso_data.query('molecule == "H2"  ')
