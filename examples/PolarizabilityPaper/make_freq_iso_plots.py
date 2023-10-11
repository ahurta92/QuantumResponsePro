import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import FancyArrow
from pathlib import Path
from quantumresponsepro import BasisMRADataAnalyzer
from quantumresponsepro import BasisMRADataCollection

august = Path('/mnt/data/madness_data/post_watoc/august')
paper_path = Path('/home/adrianhurtado/projects/writing/mra-tdhf-polarizability/Figures_v2')
thesis_path = Path('/home/adrianhurtado/projects/writing/thesis2023/Figures_v2')
paper_path = Path('/home/adrianhurtado/projects/writing/mra-tdhf-polarizability/Figures_v2')
paper_path = thesis_path
database = BasisMRADataCollection(august)
analyzer = BasisMRADataAnalyzer(database, .02)
sns.set(font_scale=3)


# e_fig = analyzer.plot_violin_strip('energy', 'error', ['D', 'T', 'Q'])
# e_fig.fig.show()
# e_fig.fig.savefig(paper_path.joinpath('energy.svg'), dpi=300)

def set_lim(ax_dict, lim_1, lim_2):
    pal2 = sns.color_palette("coolwarm_r", 4).as_hex()
    p2 = [pal2[1], pal2[2], pal2[0], pal2[3]]
    light_pal = sns.color_palette(p2)
    b = 0
    for i, ax in ax_dict.items():
        ax.set_ylim(lim_1[b], lim_2[b])
        if b != 3:
            ax.axhline(lim_1[3], color='black', linewidth=1.5, zorder=0, )
            ax.axhline(lim_2[3], color='black', linewidth=1.5, zorder=0, )
        else:
            ax.axhline(lim_1[3] + .05, color='black', linewidth=2, zorder=0, )
            ax.axhline(lim_2[3] - .03, color='black', linewidth=2, zorder=0, )
        back_color = light_pal[b]
        back_color = back_color + (0.3,)
        ax.set_facecolor(back_color)

        b += 1


def add_arrows(gfig):
    # Get the figure and axes from the FacetGrid
    # plt.tight_layout()
    fig = gfig.fig
    # Add horizontal arrow
    h_arrow_x1 = 0.35
    h_arrow_x2 = 0.80
    h_arrow_y = 0.950
    arrow_text_y = h_arrow_y - .005
    # Add caption for the horizontal arrow
    fig.text(0.5, arrow_text_y, 'Double-Augmentation', ha='center', )
    h_arrow_dy = 0.05
    arrow_color = 'magenta'
    arrow_props = dict(facecolor=arrow_color, edgecolor=arrow_color, width=0.030, )

    # Create a custom arrow that spans the entire width of the figure
    horizontal_arrow = FancyArrow(h_arrow_x1, h_arrow_y, h_arrow_x2 - h_arrow_x1, 0,
                                  transform=fig.transFigure, figure=fig,
                                  **arrow_props)

    # Add the custom arrow to the figure
    fig.add_artist(horizontal_arrow)

    # Add vertical arrow and text
    v_arrow_x = .020
    v_arrow_y1 = 0.85
    v_arrow_y2 = 0.25
    fig.text(v_arrow_x - 0.009, 0.5, 'Polarized-Core', va='center', rotation=270, )

    # Create a custom arrow that spans the entire height of the figure
    arrow_props = dict(facecolor=arrow_color, edgecolor=arrow_color, width=0.015, )
    vertical_arrow = FancyArrow(v_arrow_x, v_arrow_y1, 0, v_arrow_y2 - v_arrow_y1,
                                transform=fig.transFigure, figure=fig,
                                **arrow_props)

    # Add the custom arrow to the figure
    fig.add_artist(vertical_arrow)

    # Add caption for the vertical arrow

    # Adjust the spacing around the plot
    # Display the plot


def set_new_legend(g):
    # Access the legend object
    g.legend.remove()
    g.add_legend(title=None, loc='upper left', ncol=3,
                 frameon=False)
    plt.tight_layout()
    # plt.subplots_adjust(left=0.1, right=0.90, top=0.90, bottom=0.10, hspace=0.1, wspace=0.1)


omega = [0, 2, 4, 6, 8]
aspect = 0.9
g = analyzer.freq_iso_plot('D', 'alpha', 'all', False, omegas=omega, pal='colorblind',
                           aspect=aspect)
ax_lim_1 = [-30, -30, -30, -3]
ax_lim_2 = [5, 5, 5, 3]
g.fig.subplots_adjust(top=0.85)
set_lim(g.axes_dict, ax_lim_1, ax_lim_2)
set_new_legend(g)
g.fig.show()
g.fig.savefig(paper_path.joinpath('alpha_freq_D.svg'), dpi=300)
#
g = analyzer.freq_iso_plot('T', 'alpha', 'all', False, border=0.25, omegas=omega, aspect=aspect)
ax_lim_1 = [-10, -10, -10, -1]
ax_lim_2 = [3, 3, 3, 1]
set_lim(g.axes_dict, ax_lim_1, ax_lim_2)
set_new_legend(g)
plt.tight_layout()
g.fig.show()
g.fig.savefig(paper_path.joinpath('alpha_freq_T.svg'), dpi=300)

g = analyzer.freq_iso_plot('Q', 'alpha', 'all', False, omegas=omega, border=.05, aspect=aspect)
ax_lim_1 = [-5, -5, -5, -.2]
ax_lim_2 = [2, 2, 2, .2]
set_lim(g.axes_dict, ax_lim_1, ax_lim_2)
set_new_legend(g)

g.fig.show()
g.fig.savefig(paper_path.joinpath('alpha_freq_Q.svg'), dpi=300)

g = analyzer.freq_iso_plot('D', 'gamma', 'all', False, border=1.0, omegas=omega)
ax_lim_1 = [-80, -80, -80, -20]
ax_lim_2 = [80, 80, 80, 20]
set_lim(g.axes_dict, ax_lim_1, ax_lim_2)
set_new_legend(g)
g.fig.show()
g.fig.savefig(paper_path.joinpath('gamma_freq_D.svg'), dpi=300)

g = analyzer.freq_iso_plot('T', 'gamma', 'all', False, border=1.0, omegas=omega)
ax_lim_1 = [-20, -20, -20, -4]
ax_lim_2 = [35, 35, 35, 4]
set_lim(g.axes_dict, ax_lim_1, ax_lim_2)
set_new_legend(g)
g.fig.show()
g.fig.savefig(paper_path.joinpath('gamma_freq_T.svg'), dpi=300)

g = analyzer.freq_iso_plot('Q', 'gamma', 'all', False, border=0.25, omegas=omega)
ax_lim_1 = [-12.7, -12.7, -12.7, -1]
ax_lim_2 = [7.7, 7.7, 7.7, 1]
set_lim(g.axes_dict, ax_lim_1, ax_lim_2)
set_new_legend(g)
g.fig.show()
g.fig.savefig(paper_path.joinpath('gamma_freq_Q.svg'), dpi=300)


def add_arrows_5z(gfig):
    # Get the figure and axes from the FacetGrid
    # plt.tight_layout()
    fig = gfig.fig
    # Add horizontal arrow
    h_arrow_x1 = 0.35
    h_arrow_x2 = 0.65
    h_arrow_y = 0.10
    h_text = h_arrow_y - 0.05
    # Add caption for the horizontal arrow
    fig.text(0.5, h_text, 'Double-Augmentation', ha='center')
    h_arrow_dy = 0.05
    arrow_props = dict(facecolor='black', edgecolor='black', width=0.005, )

    # Create a custom arrow that spans the entire width of the figure
    horizontal_arrow = FancyArrow(h_arrow_x1, h_arrow_y, h_arrow_x2 - h_arrow_x1, 0,
                                  transform=fig.transFigure, figure=fig,
                                  **arrow_props)

    # Add the custom arrow to the figure
    fig.add_artist(horizontal_arrow)

    # Adjust the spacing around the plot
    plt.subplots_adjust(left=0.05, right=0.95, top=0.90, bottom=0.25, hspace=0.1, wspace=0.10)
    # Display the plot


g = analyzer.freq_iso_plot('5', 'alpha', 'all', False, border=1.0, omegas=omega, pC=False)
ax_lim_1 = [-1.5, -.1, -1, -.1]
ax_lim_2 = [.1, 0.1, 7.7, .1]
set_lim(g.axes_dict, ax_lim_1, ax_lim_2)
set_new_legend(g)
# g.axes[1, 0].remove()
# g.axes[1, 1].remove()
# Remove the second row
g.fig.show()
g.fig.savefig(paper_path.joinpath('alpha_freq_5.svg'), dpi=300, bbox_inches='tight')

g = analyzer.freq_iso_plot('5', 'gamma', 'all', False, border=0.1, omegas=omega, pC=False)
ax_lim_1 = [-2, -.2, -1, -1]
ax_lim_2 = [2, 0.2, 7.7, 1]
set_lim(g.axes_dict, ax_lim_1, ax_lim_2)
set_new_legend(g)
# g.axes[1, 0].remove()
# g.axes[1, 1].remove()
# Remove the second row
g.fig.show()
g.fig.savefig(paper_path.joinpath('gamma_freq_5.svg'), dpi=300, bbox_inches='tight')
