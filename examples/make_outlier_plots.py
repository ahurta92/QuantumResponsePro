from quantumresponsepro import BasisMRADataCollection
from quantumresponsepro import BasisMRADataAnalyzer
import seaborn as sns

from pathlib import Path

august = Path('/mnt/data/madness_data/post_watoc/august')
mra_paper_path = Path('/home/adrianhurtado/projects/writing/mra-tdhf-polarizability/Figures_v2')
thesis_path = Path('/home/adrianhurtado/projects/writing/thesis2023/Figures_v2')
paper_path = Path('/home/adrianhurtado/projects/writing/mra-tdhf-polarizability/Figures_v2')
#paper_path = thesis_path
database = BasisMRADataCollection(august)
analyzer = BasisMRADataAnalyzer(database, .05, font_scale=1.5)

percent = 0.85

sharey = False

omegas = [8]
size = (13, 11)
fig = analyzer.plot_valence_outliers_large('energy', [0, 8], ['D', 'T', 'Q'], percent, 'row', size)
fig.savefig(paper_path.joinpath('energy_outliers.svg'), dpi=300)
sns.despine(fig)
fig.show()
#
fig = analyzer.plot_valence_outliers_large('alpha', omegas, ['D', 'T', 'Q'], percent, sharey, size,
                                           .5)
sns.despine(fig)
fig.show()
fig.savefig(paper_path.joinpath('alpha_outliers.svg'), dpi=300)

size = (13, 11)
fig = analyzer.plot_valence_outliers_large('gamma', omegas, ['D', 'T', 'Q'], percent, sharey, size,
                                           1)
sns.despine(fig)
fig.show()
fig.savefig(paper_path.joinpath('gamma_outliers.svg'), dpi=300)
