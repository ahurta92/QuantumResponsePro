import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from pathlib import Path

from quantumresponsepro import BasisMRADataAnalyzer
from quantumresponsepro import BasisMRADataCollection

# Create some example data

august = Path('/mnt/data/madness_data/post_watoc/august')
paper_path = Path('/home/adrianhurtado/projects/writing/mra-tdhf-polarizability/Figures_v2')
thesis_path = Path('/home/adrianhurtado/projects/writing/thesis2023/Figures_v2')
paper_path = Path('/home/adrianhurtado/projects/writing/mra-tdhf-polarizability/Figures_v2')

tromso_path = Path('/home/adrianhurtado/projects/writing/tromso_poster/figures')

paper_path = paper_path

database = BasisMRADataCollection(august)
analyzer = BasisMRADataAnalyzer(database, .02)

sns.set_context('paper')
sns.set_theme(style="whitegrid", font_scale=2.0, rc={'ytick.left': True, 'xtick.bottom': True, })


def set_face_color(ax_dict):
    pal2 = sns.color_palette("coolwarm_r", 4).as_hex()
    p2 = [pal2[1], pal2[0], pal2[2], pal2[3]]
    light_pal = sns.color_palette(p2)
    b = 0
    for i, ax in ax_dict.items():
        back_color = light_pal[b % 4]
        back_color = back_color + (0.3,)
        ax.set_facecolor(back_color)

        b += 1


def set_ax_inset(g: sns.FacetGrid, loc='upper right', ylims=[3, .8, .18], iso_type='alpha',
                 width='80%', height='50%'):
    b = 0
    ax_dict = g.axes_dict
    ylims = [5, 1, .10]
    ylims = dict(zip(['D', 'T', 'Q'], ylims))
    print(ylims)

    for i, ax in ax_dict.items():
        vlevel = i[0]
        btype = i[1]
        yb = ylims[vlevel]
        if i[1] == 'd-aug-cc-pCVnZ' or i[1] == 'd-aug-cc-pVnZ' or (
                (vlevel == 'Q') and (i[1] == 'aug-cc-pCVnZ' or btype == 'aug-cc-pVnZ')):
            # Create the inset_axes instance
            ax_inset = inset_axes(ax, width=width, height=height, loc=loc)

            # Plot the same data in the inset_axes but zoom in on the specified range
            zoom_range = (1.50, 2.50)
            # ax_inset.set_xlim(*zoom_range)
            analyzer.iso_plot_ax(ax_inset, v_level=vlevel, b_type=btype, iso_type=iso_type,
                                 omegas=[8],
                                 pal='colorblind', )
            # ax_inset.set_title('Zoomed Inset')

            # Optionally, draw a rectangle in the main plot to indicate the zoomed area
            rect = plt.Rectangle((zoom_range[0], -yb), zoom_range[1] - zoom_range[0], 2 * yb, lw=1,
                                 edgecolor='red',
                                 linestyle='--', facecolor='none')
            ax.add_patch(rect)
            ax_inset.set_ylim(-yb, yb)
            b += 1
        else:
            ax.axhline(yb, color='k', linestyle='--', linewidth=1)
            ax.axhline(yb, color='k', linestyle='--', linewidth=1)


linthresh = .01

linscale = .5
subticks = [2, 3, 4, 5, 6, 7, 8, 9]
omega = [0, 8]
aspect = 0.9
g = analyzer.freq_iso_plot_v2(['D', 'T', 'Q'], 'alpha', 'all', False, omegas=omega,
                              pal='colorblind', aspect=1.3
                              )
g.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.5)
# g.add_legend(title=None, ncol=3, loc='lower center', frameon=False, borderaxespad=0.0, )

for ax in g.axes_dict.values():
    ax.axhline(-analyzer.mra_ref, color='green', linestyle='--')
    ax.axhline(analyzer.mra_ref, color='green', linestyle='--')
    # set a linear scale for the D and T levels
    # ax.xaxis.grid(True, "minor", linewidth=0.25)
    # ax.yaxis.grid(True, "minor", linewidth=0.25)
    positions = ax.get_xticks()
    new_labels = ['0', r'$\omega_{max}$']
    ax.set_xticklabels(new_labels)

    # Set the labels for the xticks
# e_fig.despine(left=True, bottom=True)`
set_face_color(g.axes_dict)
# set_ax_inset(g, loc='lower right', iso_type='alpha', width='30%', height='40%')
g.fig.savefig(paper_path.joinpath('alpha_freq_DTQ_2.svg'), dpi=300)
# set_ax_inset(g, loc='lower right', iso_type='alpha', width='50%', height='50%')
g.fig.show()
molecules = ['Be', 'H2O', 'HF', 'HCl', 'Na2', 'SiH4',
             'NaCl', 'SiH3Cl', 'BeH2', 'BH3', 'Ne', 'F2']

sns.set_theme(style="darkgrid")

g = analyzer.freq_iso_plot_v2_molecules(v_level=['D', 'T', 'Q'], iso_type='alpha',
                                        mol_set=molecules,
                                        sharey=False,
                                        omegas=[0, 1, 2, 3, 4, 5, 6, 7, 8],
                                        height=4,
                                        pal='colorblind', aspect=1.3,
                                        )

g.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.5)
g.add_legend(title='Molecules', loc='center right', frameon=True, )

# set_ax_inset(g, loc='lower right', iso_type='alpha', width='30%', height='40%')
g.fig.savefig(paper_path.joinpath('mol_freq_DTQ.svg'), dpi=300)
# set_ax_inset(g, loc='lower right', iso_type='alpha', width='50%', height='50%')
g.fig.show()
