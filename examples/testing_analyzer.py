from quantumresponsepro import BasisMRADataCollection
from quantumresponsepro import BasisMRADataAnalyzer
import seaborn as sns

from pathlib import Path

august = Path('/mnt/data/madness_data/post_watoc/august')
paper_path = Path('/home/adrianhurtado/projects/writing/mra-tdhf-polarizability/Figures_v2')
database = BasisMRADataCollection(august)
analyzer = BasisMRADataAnalyzer(database, .02)

# e_fig = analyzer.plot_violin_strip('energy', 'error', ['D', 'T', 'Q'])
# e_fig.fig.show()
# e_fig.fig.savefig(paper_path.joinpath('energy.svg'), dpi=300)

size = (13, 11)
fig = analyzer.plot_valence_outliers_large('energy', 0, ['D', 'T', 'Q'], .85, 'row', size)
fig.savefig(paper_path.joinpath('energy_outliers.svg'), dpi=300)
sns.despine(fig)
fig.show()
#

ax_lim_1 = [-20, -20, -2.1, -2.1]
ax_lim_2 = [2, 2, 1.2, 1.2]

a_fig = analyzer.plot_violin_strip('response', 'alpha', ['D', 'T', 'Q', '5'], sharey=False)
b = 0
for i, ax in a_fig.axes_dict.items():
    ax.set_ylim(ax_lim_1[b], ax_lim_2[b])
    ax.axhline(ax_lim_1[2], color='black', linewidth=1, alpha=0.9, zorder=0, dashes=(2, 1))
    ax.axhline(ax_lim_2[2], color='black', linewidth=1, alpha=0.9, zorder=0, dashes=(2, 1))
    b += 1

a_fig.fig.show()
a_fig.savefig(paper_path.joinpath('alpha.svg'), dpi=300)

#
# a_fig = analyzer.plot_violin_strip('response', 'alpha', ['T', 'Q', '5'])
# for i, ax in a_fig.axes_dict.items():
#     ax.set_ylim(-1.5, 1.1)
# a_fig.fig.show()
# a_fig.savefig(paper_path.joinpath('alpha_zoom.svg'), dpi=300)
#
# size = (13, 11)
# fig = analyzer.plot_valence_outliers_large('alpha', 0, ['D', 'T', 'Q'], .85, False, size, .5)
# sns.despine(fig)
# fig.show()
# fig.savefig(paper_path.joinpath('alpha_outliers.svg'), dpi=300)
#
g = analyzer.freq_iso_plot('T', 'alpha', 'all', False, border=0.25)
g.fig.show()
g.fig.savefig(paper_path.joinpath('alpha_freq_T.svg'), dpi=300)

g = analyzer.freq_iso_plot('Q', 'alpha', 'all', False, border=0.25)
g.fig.show()
g.fig.savefig(paper_path.joinpath('alpha_freq_Q.svg'), dpi=300)
#
# g = analyzer.freq_iso_plot('Q', "alpha", ['First-row', 'Fluorine'], False, thresh=.1, border=.50)
# sns.move_legend(g, "center right", title='Molecule', fancybox=True)
# g.fig.show()
# g.fig.savefig(paper_path.joinpath('alpha_freq_first_row.svg'), dpi=300)
#
# g = analyzer.freq_iso_plot('Q', "alpha", ['Second-row'], 'row', thresh=.1, border=.50)
# sns.move_legend(g, "center right", title='Molecule', fancybox=True)
# g.fig.show()
# g.fig.savefig(paper_path.joinpath('alpha_freq_second_row.svg'), dpi=300)

e_fig = analyzer.plot_violin_strip('response', 'gamma', ['D', 'T', 'Q', '5'], sharey=False)

ax_lim_1 = [-20, -20, -5.0, -5.0]
ax_lim_2 = [60, 60, 5, 5]

b = 0
for i, ax in e_fig.axes_dict.items():
    ax.set_ylim(ax_lim_1[b], ax_lim_2[b])
    ax.axhline(ax_lim_1[2], color='black', linewidth=1, alpha=0.9, zorder=0, dashes=(2, 1))
    ax.axhline(ax_lim_2[2], color='black', linewidth=1, alpha=0.9, zorder=0, dashes=(2, 1))
    b += 1

e_fig.fig.show()
e_fig.fig.savefig(paper_path.joinpath('gamma.svg'), dpi=300)

size = (13, 11)
fig = analyzer.plot_valence_outliers_large('gamma', 0, ['D', 'T', 'Q'], .85, False, size, 1)
sns.despine(fig)
fig.show()
fig.savefig(paper_path.joinpath('gamma_outliers.svg'), dpi=300)

g = analyzer.freq_iso_plot('T', 'gamma', 'all', False, border=1)

ax_lim_1 = [-20, -20, -20, -5]
ax_lim_2 = [35, 35, 35, 5]

b = 0
for i, ax in g.axes_dict.items():
    ax.set_ylim(ax_lim_1[b], ax_lim_2[b])
    ax.axhline(ax_lim_1[3], color='black', linewidth=1, alpha=0.9, zorder=0, dashes=(2, 1))
    ax.axhline(ax_lim_2[3], color='black', linewidth=1, alpha=0.9, zorder=0, dashes=(2, 1))
    b += 1
g.fig.show()
g.fig.savefig(paper_path.joinpath('gamma_freq_T.svg'), dpi=300)

g = analyzer.freq_iso_plot('Q', 'gamma', 'all', False, border=1)

ax_lim_1 = [-12.7, -12.7, -12.7, -5]
ax_lim_2 = [7.7, 7.7, 7.7, 5]

b = 0
for i, ax in g.axes_dict.items():
    ax.set_ylim(ax_lim_1[b], ax_lim_2[b])
    ax.axhline(ax_lim_1[3], color='black', linewidth=1, alpha=0.9, zorder=0, dashes=(2, 1))
    ax.axhline(ax_lim_2[3], color='black', linewidth=1, alpha=0.9, zorder=0, dashes=(2, 1))
    b += 1
g.fig.show()
g.fig.savefig(paper_path.joinpath('gamma_freq_Q.svg'), dpi=300)
