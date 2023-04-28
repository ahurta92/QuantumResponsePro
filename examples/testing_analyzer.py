from quantumresponsepro import BasisMRADataCollection
from quantumresponsepro import BasisMRADataAnalyzer
import seaborn as sns

from pathlib import Path

august = Path('/mnt/data/madness_data/post_watoc/august')
paper_path = Path('/home/adrianhurtado/projects/writing/mra-tdhf-polarizability/Figures_v2')
database = BasisMRADataCollection(august)
analyzer = BasisMRADataAnalyzer(database, .02)

e_fig = analyzer.plot_violin_strip('energy', 'error', ['D', 'T', 'Q'])
e_fig.fig.show()
e_fig.fig.savefig(paper_path.joinpath('energy.svg'), dpi=300)

size = (13, 11)
fig = analyzer.plot_valence_outliers_large('energy', 0, ['D', 'T', 'Q'], .85, 'row', size)
fig.savefig(paper_path.joinpath('energy_outliers.svg'), dpi=300)
sns.despine(fig)
fig.show()

a_fig = analyzer.plot_violin_strip('response', 'alpha', ['D', 'T', 'Q'])
a_fig.fig.show()
a_fig.savefig(paper_path.joinpath('alpha.svg'), dpi=300)

a_fig = analyzer.plot_violin_strip('response', 'alpha', ['T', 'Q', '5'])
for i, ax in a_fig.axes_dict.items():
    ax.set_ylim(-1.5, 1.1)
a_fig.fig.show()
a_fig.savefig(paper_path.joinpath('alpha_zoom.svg'), dpi=300)

size = (13, 11)
fig = analyzer.plot_valence_outliers_large('alpha', 0, ['D', 'T', 'Q'], .85, False, size, .5)
sns.despine(fig)
fig.show()
fig.savefig(paper_path.joinpath('alpha_outliers.svg'), dpi=300)

g = analyzer.freq_iso_plot('Q', 'alpha', 'all', False, border=0.5)
g.fig.show()
g.fig.savefig(paper_path.joinpath('alpha_freq.svg'), dpi=300)

g = analyzer.freq_iso_plot('Q', "alpha", ['First-row', 'Fluorine'], False, thresh=.1, border=.50)
sns.move_legend(g, "center right", title='Molecule', fancybox=True)
g.fig.show()
g.fig.savefig(paper_path.joinpath('alpha_freq_first_row.svg'), dpi=300)

g = analyzer.freq_iso_plot('Q', "alpha", ['Second-row'], 'row', thresh=.1, border=.50)
sns.move_legend(g, "center right", title='Molecule', fancybox=True)
g.fig.show()
g.fig.savefig(paper_path.joinpath('alpha_freq_second_row.svg'), dpi=300)

e_fig = analyzer.plot_violin_strip('response', 'gamma', ['D', 'T', 'Q'])
e_fig.fig.show()

size = (13, 11)
fig = analyzer.plot_valence_outliers_large('gamma', 0, ['D', 'T', 'Q'], .85, 'row', size)
sns.despine(fig)
fig.show()

g = analyzer.freq_iso_plot('Q', 'gamma', 'all', False, border=0.5)
g.fig.show()

g = analyzer.freq_iso_plot('Q', "gamma", ['First-row', 'Fluorine'], False, thresh=.1, border=.50)
sns.move_legend(g, "center right", title='Molecule', fancybox=True)
g.fig.show()

g = analyzer.freq_iso_plot('Q', "gamma", ['Second-row'], False, thresh=.1, border=.50)
sns.move_legend(g, "center right", title='Molecule', fancybox=True)
g.fig.show()
