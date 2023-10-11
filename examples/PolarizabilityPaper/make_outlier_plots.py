import seaborn as sns
from pathlib import Path
from quantumresponsepro import BasisMRADataAnalyzer
from quantumresponsepro import BasisMRADataCollection

august = Path('/mnt/data/madness_data/post_watoc/august')
mra_paper_path = Path('/home/adrianhurtado/projects/writing/mra-tdhf-polarizability/Figures_v2')
thesis_path = Path('/home/adrianhurtado/projects/writing/thesis2023/Figures_v2')
paper_path = Path('/home/adrianhurtado/projects/writing/mra-tdhf-polarizability/Figures_v2')
# paper_path = thesis_path
database = BasisMRADataCollection(august)
analyzer = BasisMRADataAnalyzer(database, .02)
sns.set_context('paper')
sns.set_theme(style="darkgrid", font_scale=1.0, rc={'ytick.left': True, 'xtick.bottom': False, })

percent = 0.85

sharey = True

omegas = [8]
size = (13, 11)
fig, axes = analyzer.plot_valence_outliers_large('energy', [0, 8], ['D', 'T', 'Q'], percent, 'row',
                                                 size)
fig.savefig(paper_path.joinpath('energy_outliers.svg'), dpi=300)
sns.despine(fig)
fig.show()
#
fig, axes = analyzer.plot_valence_outliers_large('alpha', omegas, ['D', 'T', 'Q'], percent, sharey,
                                                 size,
                                                 .5)
for ax in fig.axes:
    ax.set_xscale('symlog', linthresh=1e-2, linscale=.50, base=10)
    # ax.xaxis.grid(True, "minor", linewidth=0.25)
    # ax.yaxis.grid(True, "minor", linewidth=0.25)
sns.despine(fig)
fig.show()
fig.savefig(paper_path.joinpath('alpha_outliers.svg'), dpi=300)

size = (13, 11)
fig, axes = analyzer.plot_valence_outliers_large('gamma', omegas, ['D', 'T', 'Q'], percent, sharey,
                                                 size,
                                                 1)
for ax in fig.axes:
    ax.set_xscale('symlog', linthresh=1e-2, linscale=.50, base=10)
sns.despine(fig)
fig.show()
fig.savefig(paper_path.joinpath('gamma_outliers.svg'), dpi=300)
