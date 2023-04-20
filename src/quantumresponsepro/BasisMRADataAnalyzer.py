from quantumresponsepro import BasisMRADataCollection
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


class BasisMRADataAnalyzer:

    def __init__(self, data_collection: BasisMRADataCollection, MRA_ref):
        self.data_collection = data_collection
        self.mra_ref = MRA_ref
        self.all_basis = ['aug-cc-pVDZ', 'aug-cc-pVTZ', 'aug-cc-pVQZ',
                          'd-aug-cc-pVDZ', 'd-aug-cc-pVTZ', 'd-aug-cc-pVQZ',
                          'aug-cc-pCVDZ', 'aug-cc-pCVTZ', 'aug-cc-pCVQZ',
                          'd-aug-cc-pCVDZ', 'd-aug-cc-pCVTZ', 'd-aug-cc-pCVQZ']

    def get_basis_iso_data(self, basis_list, frequency=0):
        """Collect basis view of data for a given frequency"""
        df = self.data_collection.iso_diff_data.query('omega==@frequency')
        alphab = {}
        gammab = {}
        for b in basis_list:
            bdata = df.query('basis==@b')
            alphab[b] = bdata.set_index('molecule').alpha
            gammab[b] = bdata.set_index('molecule').gamma
        alpha_df = pd.DataFrame(alphab)
        gamma_df = pd.DataFrame(gammab)
        return alpha_df, gamma_df

    def get_basis_energy_data(self, basis_list):
        """
        Collect basis view of energy data
        """
        df = self.data_collection.energy_df
        e = {}
        for b in basis_list:
            bdata = df.query('basis==@b')
            e[b] = bdata.set_index('molecule').error
        e_df = pd.DataFrame(e)
        return e_df

    def plot_violin_strip(self, calculation_type, iso_type, valence_level, frequency=0):
        if calculation_type == 'response':
            iso_data_type = self.data_collection.detailed_iso_diff.query('omega==@frequency')
        elif calculation_type == 'energy':
            iso_data_type = self.data_collection.detailed_energy_diff

        height = 6
        aspect = 0.8

        pal = sns.color_palette("seismic_r", 4).as_hex()
        p1 = [pal[1], pal[0], pal[2], pal[3]]
        pal = sns.color_palette(p1)

        pal2 = sns.color_palette("coolwarm_r", 4).as_hex()
        p2 = [pal2[1], pal2[0], pal2[2], pal2[3]]
        light_pal = sns.color_palette(p2)

        g = sns.catplot(
            data=iso_data_type,
            height=height,
            col="valence",
            col_order=valence_level,
            y=iso_type,
            hue="Type",
            kind='strip',
            # col='mol_system',
            sharey=True,
            dodge=True,
            aspect=aspect,
            # jitter=True,
            alpha=.7,
            palette=pal,
            legend_out=False,
            legend=False
        )
        g.map_dataframe(sns.violinplot,
                        col="valence",
                        col_order=valence_level,
                        y=iso_type,
                        hue="Type",
                        palette=light_pal,
                        scale='count',
                        bw=.20,
                        inner='quartile',
                        alpha=0.5,
                        saturation=.8,
                        legend_out=False,
                        dodge=True)
        j = 0;
        for i, ax in g.axes_dict.items():
            ax.grid(visible=True, which="both", axis='y')
            ax.minorticks_on()
            for l in ax.lines:
                l.set_linestyle('--')
                l.set_linewidth(1.2)
                l.set_color('green')
                l.set_alpha(1.0)
            for l in ax.lines[1::3]:
                l.set_linestyle('-')
                l.set_linewidth(1.5)
                l.set_color('black')
                l.set_alpha(1.0)
            if j == 0:
                ax.set_ylabel('Percentage Error')
            if j == 2:
                # ax.set_ylim(-.5, .5)
                pass
            j += 1
            g.map(plt.axhline, y=0, color='k', linewidth=1.5, zorder=0). \
                set_titles("{col_name}Z")
        return g

    def iso_quartile_plot(self, iso_type, frequency, basis, ax, Q):
        iso_data = self.data_collection.detailed_iso_diff.query('omega==@frequency')
        a1, g1 = self.get_basis_iso_data(self.all_basis)

        a75 = a1.abs().quantile(q=Q)[basis]
        inner = set(iso_data.query('basis==@basis & alpha.abs() < @a75').molecule.unique())
        outer = set(iso_data.query('basis==@basis & alpha.abs() > @a75').molecule.unique())
        q4 = inner.union(outer)
        data = iso_data.query('basis==@basis').copy()
        data["dist_region"] = np.nan
        data.loc[data["molecule"].isin(outer), "dist_region"] = 'Outlier'
        data['molecule'] = data['molecule'].apply(format_formula)
        data = data.sort_values('dist_region').sort_values(iso_type)
        for i, m in enumerate(['Outlier']):
            dm = data.query('dist_region==@m')
            ax.set_title(basis)
            ax.set_ylabel('')
            ax.tick_params(
                axis='x',  # changes apply to the x-axis
                which='both',  # both major and minor ticks are affected
                bottom=True,  # ticks along the bottom edge are off
                top=False,  # ticks along the top edge are off
                labelbottom=True)  # labels
            ax.set_ylabel('')
            ax.set(xlabel=None)
            ax.set(ylabel=None)
            ax.set(xlabel='Percent Error')

            ax.axvline(x=0.0, color='k', ls='--', alpha=.5)
            ax.axvline(x=self.mra_ref, color='k', ls='--', alpha=.5)
            ax.axvline(x=-self.mra_ref, color='k', ls='--', alpha=.5)
            ax.axvline(x=.5, color='k', ls='--', alpha=.5)
            ax.axvline(x=-.5, color='k', ls='--', alpha=.5)

            sns.histplot(data=dm, y='molecule', hue='mol_system', weights=iso_type, ax=ax,
                         discrete=True,
                         binwidth=1, fill=True)
            l = ax.get_legend_handles_labels()
            # Put a legend below current axis
            ax.set(ylabel=None)

    def plot_valence_outliers_large(self, iso_type, frequency, nl, quartile, share_x, size: tuple):
        width = size[0]
        height = size[1]

        sns.set_theme(context='paper',
                      style='darkgrid', font='sans-serif',
                      font_scale=1.5, color_codes=True, rc=None)
        fig, axes = plt.subplots(nrows=3, ncols=4, figsize=(width, height),
                                 squeeze=False, layout='constrained', sharey=False, sharex=share_x)
        z = 0
        for i, n in enumerate(nl):
            DZ = ['aug-cc-pV{}Z'.format(n), 'aug-cc-pCV{}Z'.format(n),
                  'd-aug-cc-pV{}Z'.format(n), 'd-aug-cc-pCV{}Z'.format(n), ]
            for j, basis in enumerate(DZ):
                axij = axes[i, j]
                self.iso_quartile_plot(iso_type, frequency, basis, axij, quartile)

                if i == 2:
                    axij.set(xlabel='Percent Error')
                else:
                    axij.set(xlabel=None)

                if z != 0:
                    axij.get_legend().remove()
                else:
                    sns.move_legend(axij, "lower center", bbox_to_anchor=(.2, 1.05), ncol=1,
                                    title=None,
                                    frameon=False)
                z += 1
        return fig

    def plot_upperQ_energy(self, basis, ax, Q):
        a1 = self.get_basis_energy_data(self.all_basis)
        a75 = a1.quantile(q=Q)[basis]
        mdata = self.data_collection.detailed_energy_diff.copy()
        outer = set(mdata.query('basis==@basis & error > @a75').molecule.unique())
        data = mdata.query('basis==@basis').copy()
        data["dist_region"] = np.nan
        data.loc[data["molecule"].isin(outer), "dist_region"] = 'Outlier'
        data['molecule'] = data['molecule'].apply(format_formula)
        data = data.sort_values('dist_region').sort_values('error')
        for i, m in enumerate(['Outlier']):
            print(m)
            dm = data.query('dist_region==@m').copy()
            dm.loc[:, 'molecule'] = dm['molecule'].cat.remove_unused_categories()
            dm['molecule'] = dm['molecule'].astype('object')
            dm = dm.sort_values('error', ascending=True)
            ax.set_title(basis)
            ax.set_ylabel('')
            ax.tick_params(
                axis='x',  # changes apply to the x-axis
                which='both',  # both major and minor ticks are affected
                bottom=True,  # ticks along the bottom edge are off
                top=False,  # ticks along the top edge are off
                labelbottom=True)  # labels
            ax.set_ylabel('')
            ax.set(xlabel=None)
            ax.set(ylabel=None)
            ax.set(xlabel='Error')

            sns.histplot(data=dm.sort_values('error'), y=dm.molecule, hue=dm.mol_system,
                         weights=dm.error,
                         ax=ax, discrete=True, binwidth=1, fill=True)
            l = ax.get_legend_handles_labels()
            # Put a legend below current axis
            ax.set(ylabel=None)

    def plot_valence_outliers_large_energy(self, nl, Q, sharex, size: tuple):
        width = size[0]
        height = size[1]

        sns.set_theme(context='paper',
                      style='darkgrid', font='sans-serif',
                      font_scale=1.5, color_codes=True, rc=None)
        fig, axes = plt.subplots(nrows=3, ncols=4, figsize=(width, height),
                                 squeeze=False, layout='constrained', sharey=False, sharex=sharex)
        z = 0
        for i, n in enumerate(nl):
            DZ = ['aug-cc-pV{}Z'.format(n), 'aug-cc-pCV{}Z'.format(n),
                  'd-aug-cc-pV{}Z'.format(n), 'd-aug-cc-pCV{}Z'.format(n), ]
            for j, basis in enumerate(DZ):
                axij = axes[i, j]
                self.plot_upperQ_energy(basis, axij, Q)

                if i == 2:
                    axij.set(xlabel='Error')
                else:
                    axij.set(xlabel=None)

                if z != 0:
                    axij.get_legend().remove()
                else:
                    sns.move_legend(axij, "lower center", bbox_to_anchor=(.2, 1.05), ncol=1,
                                    title=None,
                                    frameon=False)
                z += 1
        return fig


b1 = [
    'aug-cc-pVDZ', 'aug-cc-pVTZ', 'aug-cc-pVQZ',
    'aug-cc-pCVDZ', 'aug-cc-pCVTZ', 'aug-cc-pCVQZ',
    'd-aug-cc-pVDZ', 'd-aug-cc-pVTZ', 'd-aug-cc-pVQZ',
    'd-aug-cc-pCVDZ', 'd-aug-cc-pCVTZ', 'd-aug-cc-pCVQZ', ]




EE = 5
cmap = 'BrBG'


def make_pretty_gamma_summary(styler):
    EE = 10

    styler.format("{:.2e}", precision=0)
    styler.background_gradient(
        axis=1, vmin=-EE, vmax=EE, cmap=cmap, low=0, high=0
    )
    # styler.apply(highlight_positive, axis=0)
    return styler


def make_pretty_summary(styler):
    styler.format("{:.2e}", precision=0)
    styler.background_gradient(
        axis=1, vmin=-EE, vmax=EE, cmap=cmap, low=0, high=0
    )
    # styler.apply(highlight_positive, axis=0)
    return styler


def make_e_pretty_summary(styler):
    EE = .05
    styler.format('{:.2e}', precision=0)
    styler.background_gradient(
        axis=1, vmin=-EE, vmax=EE, cmap=cmap, low=0, high=0
    )
    # styler.apply(highlight_positive, axis=0)
    return styler


def max_abs(row):
    max_val = max(abs(row['alpha']), abs(row['B']))
    orig_val = row['A'] if abs(row['A']) == max_val else row['B']
    return pd.Series({'max_abs_val': max_val, 'orig_val': orig_val})


# create a sample dataframe with a column of chemical formulas

# define a function to format the chemical formulas
def format_formula(formula):
    formatted_formula = ''
    for c in formula:
        if c.isdigit():
            # use Unicode subscript characters for digits
            formatted_formula += chr(8320 + int(c))
        else:
            formatted_formula += c
    return formatted_formula
