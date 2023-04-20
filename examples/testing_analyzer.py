from quantumresponsepro import BasisMRADataCollection
from quantumresponsepro import BasisMRADataAnalyzer
import seaborn as sns

from pathlib import Path

august = Path('/mnt/data/madness_data/post_watoc/august')
database = BasisMRADataCollection(august)
analyzer = BasisMRADataAnalyzer(database, .05)

e_fig = analyzer.plot_violin_strip('energy', 'error', ['D', 'T', 'Q'])
e_fig.fig.show()

size = (13, 11)
fig = analyzer.plot_valence_outliers_large_energy(['D', 'T', 'Q'], .85, False, size)
sns.despine(fig)
fig.show()

a_fig = analyzer.plot_violin_strip('response', 'alpha', ['D', 'T', 'Q'])
a_fig.fig.show()

a_fig = analyzer.plot_violin_strip('response', 'alpha', ['T', 'Q'])
a_fig.fig.show()

size = (13, 11)
fig = analyzer.plot_valence_outliers_large('alpha', 0, ['D', 'T', 'Q'], .85, False, size)


# g_fig = analyzer.plot_violin_strip('response', 'gamma', ['D', 'T', 'Q'])
# g_fig.fig.show()

sns.despine(fig)
fig.show()
