from quantumresponsepro import BasisMRADataCollection
from quantumresponsepro import BasisMRADataAnalyzer
import seaborn as sns

from pathlib import Path
import matplotlib
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.pyplot as plt

august = Path('/mnt/data/madness_data/post_watoc/august')
paper_path = Path('/home/adrianhurtado/projects/writing/mra-tdhf-polarizability/Figures_v2')
thesis_path = Path('/home/adrianhurtado/projects/writing/thesis2023/Figures_v2')
paper_path = Path('/home/adrianhurtado/projects/writing/mra-tdhf-polarizability/Figures_v2')
#paper_path = thesis_path
database = BasisMRADataCollection(august)
analyzer = BasisMRADataAnalyzer(database, .05, font_scale=3)


def set_ax_inset(g: sns.FacetGrid, loc='upper right', iso_type='alpha', TZ_lim=None,
                 calc_type='response',
                 width='80%', height='50%'):
    b = 0
    ax_dict = g.axes_dict

    for i, ax in ax_dict.items():
        if i == 'Q':
            if TZ_lim:
                height = TZ_lim[2]
                yb = TZ_lim[0]
            # Create the inset_axes instance
            ax_inset = inset_axes(ax, width=width, height=height, loc=loc)
            vlevel = i[0]
            analyzer.plot_violin_strip_ax(ax_inset, v_level=vlevel, calculation_type=calc_type,
                                          iso_type=iso_type, )
            if TZ_lim:
                ax_inset.set_ylim(-yb, yb)
            for spine in ax_inset.spines.values():
                spine.set_linewidth(1)
                spine.set_color('black')
            ax.legend().remove()
            b += 1
        if i == 'T' and TZ_lim:
            yb = TZ_lim[0]
            t_width = TZ_lim[1]
            t_height = TZ_lim[2]
            ax_inset = inset_axes(ax, width=t_width, height=t_height, loc=loc)
            ax_inset.set_ylim(-yb, yb)
            ax_inset.set_xlim(0, 1.5)
            vlevel = i[0]
            analyzer.plot_violin_strip_ax(ax_inset, v_level=vlevel, calculation_type=calc_type,
                                          iso_type=iso_type, )

            for spine in ax_inset.spines.values():
                spine.set_linewidth(1)
                spine.set_color('black')
            ax.legend().remove()
            b += 1


def set_lim(ax_dict, lim_1, lim_2):
    b = 0
    for i, ax in ax_dict.items():
        ax.set_ylim(lim_1[b], lim_2[b])
        ax.axhline(lim_1[3], color='black', linewidth=1, alpha=0.9, zorder=0, dashes=(5, 0))
        ax.axhline(lim_2[3], color='black', linewidth=1, alpha=0.9, zorder=0, dashes=(5, 0))
        b += 1


e_fig = analyzer.plot_violin_strip('energy', 'error', ['D', 'T', 'Q'], sharey=True)
set_ax_inset(e_fig, loc='upper right', iso_type='error', calc_type='energy', width='70%')
e_fig.fig.show()
e_fig.fig.savefig(paper_path.joinpath('energy.svg'), dpi=300)

a_fig = analyzer.plot_violin_strip('response', 'alpha', ['D', 'T', 'Q'], sharey=True)
set_ax_inset(a_fig, loc='lower right', TZ_lim=(2, '40%', '70%'), iso_type='alpha',
             calc_type='response',
             height='55%',
             width='70%')

a_fig.fig.show()
a_fig.savefig(paper_path.joinpath('alpha.svg'), dpi=300)

e_fig = analyzer.plot_violin_strip('response', 'gamma', ['D', 'T', 'Q'], sharey=True)
set_ax_inset(e_fig, loc='upper right', iso_type='gamma', calc_type='response',
             height='80%', width='50%', )
e_fig.fig.show()
e_fig.fig.savefig(paper_path.joinpath('gamma.svg'), dpi=300)
