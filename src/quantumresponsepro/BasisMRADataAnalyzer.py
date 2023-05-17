from quantumresponsepro import BasisMRADataCollection
from quantumresponsepro.BasisMRADataAssembler import get_frequency_compare
from quantumresponsepro.BasisMRADataAssembler import partition_molecule_list
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import pandas as pd
import matplotlib.colors as mcolors


class BasisMRADataAnalyzer:

    def __init__(self, data_collection: BasisMRADataCollection, MRA_ref, font_scale=3):
        self.data_collection = data_collection
        self.mra_ref = MRA_ref
        self.all_basis = ['aug-cc-pVDZ', 'aug-cc-pVTZ', 'aug-cc-pVQZ',
                          'd-aug-cc-pVDZ', 'd-aug-cc-pVTZ', 'd-aug-cc-pVQZ',
                          'aug-cc-pCVDZ', 'aug-cc-pCVTZ', 'aug-cc-pCVQZ',
                          'd-aug-cc-pCVDZ', 'd-aug-cc-pCVTZ', 'd-aug-cc-pCVQZ']

        sns.set_theme(context='paper',
                      font='sans-serif',
                      font_scale=font_scale, color_codes=True, rc=None)

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
        df = self.data_collection.energy_diff
        e = {}
        for b in basis_list:
            bdata = df.query('basis==@b')
            e[b] = bdata.set_index('molecule').error
        e_df = pd.DataFrame(e)
        return e_df

    def plot_violin_strip_ax(self, ax, calculation_type, iso_type, v_level, frequency=0,
                             sharey=True,
                             sharex=False, border=0.0):
        if calculation_type == 'response':
            iso_data_type = self.data_collection.detailed_iso_diff.query(
                'omega==@frequency & valence==@v_level').copy()

        elif calculation_type == 'energy':
            iso_data_type = self.data_collection.detailed_energy_diff.query(
                'valence==@v_level').copy()
        iso_data_type['valence'] = iso_data_type['valence'].cat.remove_unused_categories()

        height = 6
        aspect = 0.8

        pal = sns.color_palette("seismic_r", 4).as_hex()
        p1 = [pal[1], pal[0], pal[2], pal[3]]
        pal = sns.color_palette(p1)

        pal2 = sns.color_palette("coolwarm_r", 4).as_hex()
        p2 = [pal2[1], pal2[0], pal2[2], pal2[3]]
        light_pal = sns.color_palette(p2)
        facet_kws = dict(sharey=sharey, sharex=sharex)

        sns.violinplot(data=iso_data_type, height=height, y=iso_type, x='valence', hue="Type",
                       inner=None, palette=light_pal, ax=ax, bw=.20, cut=0, scale='count',
                       saturation=.8,
                       legend=False,
                       legend_out=False,
                       dodge=True)
        sns.stripplot(data=iso_data_type, y=iso_type, x='valence', hue="Type", palette=pal,
                      dodge=True,
                      legend=False, )

        ax.axhline(y=0.0, color='k', ls='--', alpha=.5)
        if calculation_type == 'response':
            ax.axhline(y=-self.mra_ref, color='k', ls='--', alpha=.5)
            ax.axhline(y=self.mra_ref, color='k', ls='--', alpha=.5)
        ax.set_xticklabels("")
        ax.set(xlabel=None)
        ax.set(ylabel=None)
        ax.get_legend().remove()

    def plot_violin_strip(self, calculation_type, iso_type, valence_level, frequency=0, sharey=True,
                          sharex=False, border=0.0):
        if calculation_type == 'response':
            iso_data_type = self.data_collection.detailed_iso_diff.query('omega==@frequency')
        elif calculation_type == 'energy':
            iso_data_type = self.data_collection.detailed_energy_diff

        sns.set(rc={"xtick.bottom": True, "ytick.left": True}, font_scale=1.5)
        pal = sns.color_palette("seismic_r", 4).as_hex()
        p1 = [pal[1], pal[0], pal[2], pal[3]]
        pal = sns.color_palette(p1)

        pal2 = sns.color_palette("coolwarm_r", 4).as_hex()
        p2 = [pal2[1], pal2[0], pal2[2], pal2[3]]
        light_pal = sns.color_palette(p2)
        facet_kws = dict(sharey=sharey, sharex=sharex)
        aspect = 0.3

        g = sns.catplot(
            aspect=aspect,
            data=iso_data_type,
            col="valence",
            col_order=valence_level,
            y=iso_type,
            hue="Type",
            kind='strip',
            sharey=sharey,
            dodge=True,
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
                        cut=0,
                        inner=None,
                        alpha=0.5,
                        saturation=.8,
                        legend_out=False,
                        dodge=True)
        j = 0;
        g.add_legend(title=None, ncol=4, loc='lower center', frameon=False, borderaxespad=0.0, )
        # sns.move_legend(g, title=None, ncol=1, loc='upper right', frameon=False,
        #                bbox_to_anchor=(1.0, 0.95), borderaxespad=0.0, )
        g.tight_layout()
        if calculation_type == 'response':
            g.set_axis_labels("", "Percent Error")
        else:
            g.set_axis_labels("", "Energy Error (a.u.)")
        g.map(plt.axhline, y=0, color='k', linewidth=3.0, alpha=.8, dashes=(2, 1)).set_titles(
            "{col_name}Z")
        plt.subplots_adjust(right=0.97, bottom=0.20)
        g.figure.subplots_adjust(wspace=0.00, hspace=0.0)
        for ax in g.axes.flat:
            for spine in ax.spines.values():
                spine.set_linewidth(1)
                spine.set_color('black')
        return g

    def iso_quartile_plot(self, iso_type, frequency, basis, ax, Q, outer_boundary=0):
        if iso_type == 'alpha' or iso_type == 'gamma':
            iso_data = self.data_collection.detailed_iso_diff.query('omega.isin(@frequency)')
            if iso_type == 'alpha':
                aQ = self.get_basis_iso_data(self.all_basis)[0][basis].abs().quantile(Q)
                outer = set(iso_data.query('basis==@basis & alpha.abs() > @aQ').molecule.unique())
            else:
                aQ = self.get_basis_iso_data(self.all_basis)[1][basis].abs().quantile(Q)
                outer = set(iso_data.query('basis==@basis & gamma.abs() > @aQ').molecule.unique())
        else:
            iso_data = self.data_collection.detailed_energy_diff.copy()
            aQ = self.get_basis_energy_data(self.all_basis)[basis].quantile(Q)
            iso_type = 'error'
            outer = set(iso_data.query('basis==@basis & error > @aQ').molecule.unique())

        data = iso_data.query('basis==@basis').copy()
        data["dist_region"] = np.nan
        data.loc[data["molecule"].isin(outer), "dist_region"] = 'Outlier'
        data['molecule'] = data['molecule'].apply(format_formula)
        data['molecule'] = data['molecule'].astype('object')
        data = data.sort_values('dist_region').sort_values(iso_type)
        dm = data.query('dist_region=="Outlier"')
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
        if iso_type == 'error':
            ax.set(xlabel='Energy Error (a.u.)')
        else:
            ax.set(xlabel='Percent Error')

        ax.axvline(x=0.0, color='k', ls='--', alpha=.5)
        ax.axvline(x=self.mra_ref, color='k', ls='--', alpha=.3, )
        if iso_type != 'error':
            ax.axvline(x=-self.mra_ref, color='k', ls='--', alpha=.3, )

        if outer_boundary != 0:
            ax.axvline(x=outer_boundary, color='black', ls='--', alpha=.3, )
            ax.axvline(x=-outer_boundary, color='black', ls='--', alpha=.3, )

        sns.histplot(data=dm, y='molecule', hue='mol_system', weights=iso_type, ax=ax,
                     discrete=True,
                     binwidth=1, fill=True)
        ax.set(ylabel=None)

    def plot_valence_outliers_large(self, iso_type, frequency, nl, quartile, share_x,
                                    size: tuple, outer_boundary: float = 0.0):
        width = size[0]
        height = size[1]

        fig, axes = plt.subplots(nrows=3, ncols=4, figsize=(width, height),
                                 squeeze=True, sharey=False, sharex=share_x)
        z = 0
        for i, n in enumerate(nl):
            DZ = ['aug-cc-pV{}Z'.format(n), 'aug-cc-pCV{}Z'.format(n),
                  'd-aug-cc-pV{}Z'.format(n), 'd-aug-cc-pCV{}Z'.format(n), ]
            for j, basis in enumerate(DZ):
                axij = axes[i, j]
                self.iso_quartile_plot(iso_type, frequency, basis, axij, quartile, outer_boundary)

                if i == 2:
                    axij.set(xlabel='Percent Error')
                else:
                    axij.set(xlabel=None)

                if z != 0:
                    axij.get_legend().remove()
                else:
                    sns.move_legend(axij, "lower center", bbox_to_anchor=(.7, 1.05), ncol=3,
                                    title=None,
                                    frameon=False)
                z += 1
        if iso_type == 'energy':
            sup_title = r'Outliers for $E$'
        elif iso_type == 'alpha':
            sup_title = r'{} percentile for $\alpha(\omega_{})$ '.format(
                quartile, frequency[0])
        else:
            sup_title = r'{} percentile for $\gamma(\omega_{})$'.format(
                quartile, frequency[0])
        fig.suptitle(sup_title, x=.8, y=.995, )

        plt.subplots_adjust(left=0.05, right=0.97, bottom=0.05, top=.95)
        # fig.subplots_adjust(wspace=0.05, hspace=0.05)
        return fig

    def iso_plot_ax(self, ax, v_level, b_type, iso_type,
                    omegas=[0, 1, 2, 3, 4, 5, 6, 7, 8], border=0.0,
                    pal='colorblind'):

        iso_diff_detailed = self.data_collection.detailed_iso_diff.copy()
        data = iso_diff_detailed.query('omega.isin(@omegas) & valence==@v_level & Type==@b_type')
        sns.stripplot(x="omega", y=iso_type, data=data, dodge=True, ax=ax, palette=pal,
                      hue='mol_system', size=2.5, jitter=True, legend=False)

        ax.axhline(y=0.0, color='k', ls='--', alpha=.5)
        ax.axhline(y=-self.mra_ref, color='k', ls='--', alpha=.5)
        ax.axhline(y=self.mra_ref, color='k', ls='--', alpha=.5)

        ax.set(xlabel=None)
        ax.set(ylabel=None)

    def cluster_iso_plot_ax(self, iso_diff_detailed, ax, v_level, b_type, iso_type,
                            omegas=[0, 1, 2, 3, 4, 5, 6, 7, 8], border=0.0,
                            pal='colorblind'):

        data = iso_diff_detailed.query('omega.isin(@omegas) & valence==@v_level & Type==@b_type')
        sns.stripplot(x="omega", y=iso_type, data=data, dodge=True, ax=ax, palette=pal,
                      hue='cluster', size=2.0, jitter=True, legend=False)

        ax.axhline(y=0.0, color='k', ls='--', alpha=.5)
        ax.axhline(y=-self.mra_ref, color='k', ls='--', alpha=.5)
        ax.axhline(y=self.mra_ref, color='k', ls='--', alpha=.5)

        ax.set(xlabel=None)
        ax.set(ylabel=None)

    def freq_iso_plot_cluster(self, iso_diff_detailed, v_level, iso_type, mol_set="all",
                              sharey=False,
                              omegas=[0, 1, 2, 3, 4, 5, 6, 7, 8], border=0.0,
                              pal='colorblind'):
        if mol_set == "all":
            data = iso_diff_detailed.query('omega.isin(@omegas) & valence.isin(@v_level)')
        mdata = data.copy()
        mdata['valence'] = data.valence.cat.remove_unused_categories()
        mdata['Type'] = data.Type.cat.remove_unused_categories()
        sns.set(rc={"xtick.bottom": True, "ytick.left": True}, font_scale=1.5)
        g = sns.FacetGrid(data=mdata, col='Type', row='valence', margin_titles=True,
                          sharex=True, sharey=sharey, despine=False, legend_out=True,
                          palette=pal)

        g.map_dataframe(sns.stripplot, x="omega", y=iso_type, hue="cluster", dodge=True,
                        size=4.25, jitter=True, legend="full",
                        palette=pal, )

        g.add_legend(title=None, ncol=3, loc='lower center', frameon=False, borderaxespad=0.0, )
        # sns.move_legend(g, title=None, ncol=1, loc='upper right', frameon=False,
        #                bbox_to_anchor=(1.0, 0.95), borderaxespad=0.0, )
        g.map(plt.axhline, y=0, color='k', dashes=(2, 1), zorder=0).tight_layout()
        g.map(plt.axhline, y=self.mra_ref, color='k', dashes=(2, 1), zorder=0).tight_layout()
        g.map(plt.axhline, y=-self.mra_ref, color='k', dashes=(2, 1), zorder=0).tight_layout()
        g.map(plt.axhline, y=1, color='k', dashes=(2, 1), zorder=0).tight_layout()
        g.map(plt.axhline, y=-1, color='k', dashes=(2, 1), zorder=0).tight_layout()
        g.set_axis_labels(r"$\omega_i$", "Percent Error")
        g.figure.subplots_adjust(wspace=0.00, hspace=0.0)
        g.set_titles(row_template="{row_name}", col_template="{col_name}")
        if border != 0.0:
            g.map(plt.axhline, y=border, color='black', dashes=(2, 1), zorder=0)
            g.map(plt.axhline, y=-border, color='black', dashes=(2, 1), zorder=0)
        if iso_type == 'alpha':
            title = r'Error in $\alpha(\omega)$'
        else:
            title = r'Error in $\gamma(\omega)$'

        for ax in g.axes.flat:
            for spine in ax.spines.values():
                spine.set_linewidth(1)
                spine.set_color('black')
        plt.subplots_adjust(right=0.97, bottom=0.12)
        g.figure.subplots_adjust(wspace=0.00, hspace=0.0)
        return g

    def freq_iso_plot_v2(self, v_level, iso_type, mol_set="all", sharey=False,
                         omegas=[0, 1, 2, 3, 4, 5, 6, 7, 8], border=0.0,
                         pal='colorblind'):
        iso_diff_detailed = self.data_collection.detailed_iso_diff.copy()
        if mol_set == "all":
            data = iso_diff_detailed.query('omega.isin(@omegas) & valence.isin(@v_level)')
        mdata = data.copy()
        mdata['valence'] = data.valence.cat.remove_unused_categories()
        mdata['Type'] = data.Type.cat.remove_unused_categories()
        sns.set(rc={"xtick.bottom": True, "ytick.left": True}, font_scale=1.5)
        g = sns.FacetGrid(data=mdata, col='Type', row='valence', margin_titles=True,
                          sharex=True, sharey=sharey, despine=False, legend_out=True,
                          palette=pal)

        g.map_dataframe(sns.stripplot, x="omega", y=iso_type, hue="mol_system", dodge=True,
                        size=4.25, jitter=True, legend="full",
                        palette=pal, )

        g.add_legend(title=None, ncol=3, loc='lower center', frameon=False, borderaxespad=0.0, )
        # sns.move_legend(g, title=None, ncol=1, loc='upper right', frameon=False,
        #                bbox_to_anchor=(1.0, 0.95), borderaxespad=0.0, )
        g.map(plt.axhline, y=0, color='k', dashes=(2, 1), zorder=0).tight_layout()
        g.set_axis_labels(r"$\omega_i$", "Percent Error")
        g.figure.subplots_adjust(wspace=0.00, hspace=0.0)
        g.set_titles(row_template="{row_name}", col_template="{col_name}")
        if border != 0.0:
            g.map(plt.axhline, y=border, color='black', dashes=(2, 1), zorder=0)
            g.map(plt.axhline, y=-border, color='black', dashes=(2, 1), zorder=0)
        if iso_type == 'alpha':
            title = r'Error in $\alpha(\omega)$'
        else:
            title = r'Error in $\gamma(\omega)$'

        for ax in g.axes.flat:
            for spine in ax.spines.values():
                spine.set_linewidth(1)
                spine.set_color('black')
        plt.subplots_adjust(right=0.97, bottom=0.12)
        g.figure.subplots_adjust(wspace=0.00, hspace=0.0)
        return g

    def freq_iso_plot(self, v_level, iso_type, mol_set="all", sharey=False,
                      thresh=.1, omegas=[0, 1, 2, 3, 4, 5, 6, 7, 8], border=0.0, pC=True,
                      aspect=1.3, pal='colorblind'):
        iso_diff_detailed = self.data_collection.detailed_iso_diff.copy()
        if mol_set == "all":
            data = iso_diff_detailed.query('omega.isin(@omegas) & valence==@v_level')
            row_order = ['V', 'CV']
        elif mol_set[0] == "First-row":
            data = iso_diff_detailed.query(
                'polarization=="V" & omega.isin(@omegas) & valence==@v_level & mol_system.isin('
                '@mol_set)').copy()
            data['polarization'] = data.polarization.cat.remove_unused_categories()
            row_order = ['V']
        else:
            row_order = ['V', 'CV']
            data = iso_diff_detailed.query(
                'omega.isin(@omegas) & valence==@v_level & mol_system.isin(@mol_set)')
        mdata = data.query('valence==@v_level')
        sns.set(rc={"xtick.bottom": True, "ytick.left": True})
        sns.set(font_scale=2.0)

        if mol_set == 'all':
            if pC == True:

                g = sns.catplot(data=mdata,
                                x="omega",
                                row="polarization",
                                row_order=row_order,
                                y=iso_type,
                                hue="mol_system",
                                col="augmentation",
                                kind='strip',
                                palette=pal,
                                legend="full",
                                sharey=sharey,
                                dodge=True,
                                height=5,
                                aspect=aspect,
                                size=5,
                                )
            else:
                g = sns.catplot(data=mdata,
                                x="omega",
                                row_order=row_order,
                                y=iso_type,
                                hue="mol_system",
                                col="augmentation",
                                kind='strip',
                                legend="full",
                                sharey=sharey,
                                dodge=True,
                                aspect=aspect
                                )
        else:
            aspect = 1.0
            set1 = self.__select_basis_outliers(mdata, 'aug-cc-pV{}Z'.format(v_level), thresh)
            set2 = self.__select_basis_outliers(mdata, 'd-aug-cc-pCV{}Z'.format(v_level),
                                                thresh * .5)
            set1 = set1.union(set2)
            mdata = mdata.query('molecule.isin(@set1)')
            mdata.dropna()
            fwk = {"sharey": sharey, 'despine': True, }
            if iso_type == 'alpha':
                sizes = mdata.alpha
            else:
                sizes = mdata.gamma

            g = sns.relplot(data=mdata,
                            x="omega",
                            row="polarization",
                            row_order=row_order,
                            y=iso_type,
                            hue="molecule",
                            style="molecule",
                            col="augmentation",
                            kind='line',
                            markers=True,
                            facet_kws=fwk,
                            palette='Paired',
                            aspect=aspect,
                            # sizes="abs"
                            )
            g.map(plt.axhline, y=0, color='k', dashes=(2, 1), zorder=0). \
                set_axis_labels(r"$\omega_i$", "Percent Error"). \
                set_titles("{col_name}-cc-p{row_name}" + "{}Z".format(v_level)).tight_layout(
                w_pad=0)
            g.fig.suptitle(iso_type)
        g.map(plt.axhline, y=0, color='k', dashes=(2, 1), zorder=1). \
            set_axis_labels(r"$\omega_i$", "Percent Error").tight_layout(w_pad=0)
        if border != 0.0:
            g.map(plt.axhline, y=border, color='black', dashes=(2, 1), zorder=0)
            g.map(plt.axhline, y=-border, color='black', dashes=(2, 1), zorder=0)
        if pC == True:
            g.set_titles("{col_name}-cc-p{row_name}" + "{}Z".format(v_level))
        else:
            g.set_titles("{col_name}-cc-pV" + "{}Z".format(v_level))

        if iso_type == 'alpha':
            title = r'Error in $\alpha(\omega)$ at {}Z'.format(v_level)
        else:
            title = r'Error in $\gamma(\omega)$ at {}Z'.format(v_level)

        g.fig.suptitle(t=title, x=0.5, y=0.93,
                       ha='left', va='bottom', fontweight='bold')

        return g

    def plot_iso_valence_convergence_inset(self, ax, mol, iso_type, valence, omega, sharey=False):
        data = self.data_collection.detailed_iso_diff.query(
            'molecule==@mol & omega.isin(@omega) & valence==@valence').copy()

        data[iso_type] = data[iso_type].abs()
        sns.lineplot(data=data,
                     x='valence',
                     hue='Type',
                     style='Type',
                     y=iso_type,
                     legend=False,
                     ax=ax)

        ax.axhline(y=0.0, color='k', ls='--', alpha=.5)
        ax.axhline(y=-self.mra_ref, color='k', ls='--', alpha=.5)
        ax.axhline(y=self.mra_ref, color='k', ls='--', alpha=.5)

        ax.set(xlabel=None)
        ax.set(ylabel=None)

    def plot_iso_valence_convergence_v2(self, mol, iso_type, valence, omega, sharey=False):
        data = self.data_collection.detailed_iso_diff.query(
            'molecule==@mol & omega.isin(@omega) & valence==@valence')
        facet_kws = {"sharey": sharey, 'despine': True, }

        f = sns.relplot(data=data,
                        x='valence',
                        kind='line',
                        style='Type',
                        hue='Type',
                        facet_kws=facet_kws,
                        y=iso_type)
        f.map(plt.axhline, y=0, color='k', dashes=(2, 1), zorder=0)
        f.set_titles("{col_name}-cc-p(C)VnZ", x=.7).tight_layout(w_pad=0)
        f.fig.suptitle(r'$\Delta\{}(\omega)$ in {}'.format(iso_type, mol))
        f.set_xlabels('Valence [n]')
        f.set_ylabels('Percent Error')
        f.tight_layout()
        return f

    def plot_iso_valence_convergence(self, mol, iso_type, valence, omega, sharey=False):
        data = self.data_collection.detailed_iso_diff.query(
            'molecule==@mol & omega.isin(@omega) & valence==@valence')
        facet_kws = {"sharey": sharey, 'despine': True, }

        f = sns.relplot(data=data,
                        x='valence',
                        kind='scatter',
                        style='polarization',
                        col='augmentation',
                        facet_kws=facet_kws,
                        y=iso_type)
        f.map(plt.axhline, y=0, color='k', dashes=(2, 1), zorder=0)
        f.set_titles("{col_name}-cc-p(C)VnZ").tight_layout(w_pad=0)
        f.fig.suptitle(r'$\Delta\{}(\omega)$ in {}'.format(iso_type, mol))
        f.set_xlabels('Valence [n]')
        f.set_ylabels('Percent Error')
        return f

    def plot_alpha_component_convergence(self, mol, ij=['xx', 'yy', 'zz'], valence=['D', 'T', 'Q'],
                                         omega=[0, 1, 2, 3, 4, 5, 6, 7, 8], sharey=False):
        mol_data = self.data_collection.detailed_eigen_diff.query(
            ' molecule == @mol & ij.isin(@ij) & valence.isin(@valence) & omega.isin(@omega)')
        g = sns.relplot(data=mol_data,
                        x=mol_data.valence,
                        y='alpha',
                        hue='ij',
                        col='augmentation',
                        style='polarization',
                        kind='line',
                        markers=True,
                        facet_kws={'sharey': sharey},
                        dashes=True,
                        )
        for i, ax in enumerate(g.axes):
            for j, axi in enumerate(ax):
                axi.tick_params(axis='x', rotation=0)
                axi.grid(which="both")
                axi.tick_params(which="both", top="on", left="on", right="on", bottom="on", )
                axi.minorticks_on()
        g.map(plt.axhline, y=0, color='k', dashes=(2, 1), zorder=0).set_axis_labels("Valence",
                                                                                    "Error").set_titles(
            "{col_name}-cc-pV(C)nZ").tight_layout(w_pad=0)
        g.fig.suptitle(mol)
        return g

    def plot_alpha_component_convergence(self, mol, ij=['xx', 'yy', 'zz'], valence=['D', 'T', 'Q'],
                                         omega=[0, 1, 2, 3, 4, 5, 6, 7, 8], sharey=False):
        mol_data = self.data_collection.detailed_eigen_diff.query(
            ' molecule == @mol & ij.isin(@ij) & valence.isin(@valence) & omega.isin(@omega)')
        g = sns.relplot(data=mol_data,
                        x=mol_data.valence,
                        y='alpha',
                        hue='ij',
                        col='augmentation',
                        style='polarization',
                        kind='line',
                        markers=True,
                        facet_kws={'sharey': False},
                        dashes=True,
                        )
        for i, ax in enumerate(g.axes):
            for j, axi in enumerate(ax):
                axi.tick_params(axis='x', rotation=0)
                axi.grid(which="both")
                axi.tick_params(which="both", top="on", left="on", right="on", bottom="on", )
                axi.minorticks_on()
        g.map(plt.axhline, y=0, color='k', dashes=(2, 1), zorder=0).set_axis_labels("Valence",
                                                                                    "Error").set_titles(
            "{col_name}-cc-pV(C)nZ").tight_layout(w_pad=0)
        g.fig.suptitle(mol)
        return g

    def plot_alpha_eigen(self, mol, ij=['xx', 'yy', 'zz'], valence=['D', 'T', 'Q'],
                         omega=[0, 1, 2, 3, 4, 5, 6, 7, 8], sharey=False):
        mol_data = self.data_collection.detailed_alpha_eigen.query(
            ' molecule == @mol & ij.isin(@ij) & valence.isin(@valence) & omega.isin(@omega)')
        g = sns.relplot(data=mol_data,
                        x=mol_data.valence,
                        y='alpha',
                        hue='ij',
                        col='augmentation',
                        style='polarization',
                        kind='line',
                        markers=True,
                        facet_kws={'sharey': False},
                        dashes=True,
                        )
        for i, ax in enumerate(g.axes):
            for j, axi in enumerate(ax):
                axi.tick_params(axis='x', rotation=0)
                axi.grid(which="both")
                axi.tick_params(which="both", top="on", left="on", right="on", bottom="on", )
                axi.minorticks_on()
        g.map(plt.axhline, y=0, color='k', dashes=(2, 1), zorder=0).set_axis_labels("Valence",
                                                                                    "Error").set_titles(
            "{col_name}-cc-pV(C)nZ").tight_layout(w_pad=0)
        g.fig.suptitle(mol)
        return g

    # selects outliers by change in percent error from the static to the
    def __select_basis_outliers(self, data, basis, thresh):
        om = [0, 8]
        test = data.query('omega.isin(@om) & basis==@basis')

        ma = {}
        for mol in test.molecule.unique():
            dmol = test.query('molecule==@mol')
            a0 = dmol.query('omega==0').alpha.iloc[0]
            a8 = dmol.query('omega==8').alpha.iloc[0]
            ma[mol] = (a0 - a8)
        ma = pd.Series(ma)
        out_mols = ma[ma.abs() > thresh]
        out_mols = set(out_mols.index)
        return out_mols


def background_with_norm(s, vmax=1e-1):
    linthresh = 1e-2
    linscale = 1
    cmap = cm.coolwarm
    norm = mcolors.SymLogNorm(linthresh=linthresh, linscale=linscale, base=10, vmin=-vmax,
                              vmax=vmax * 2)
    return ['background-color: {:s}'.format(mcolors.to_hex(c.flatten())) for c in
            cmap(norm(s.values))]


DZ = ['aug-cc-pVDZ', 'aug-cc-pCVDZ', 'd-aug-cc-pVDZ', 'd-aug-cc-pCVDZ']
TZ = ['aug-cc-pVTZ', 'aug-cc-pCVTZ', 'd-aug-cc-pVTZ', 'd-aug-cc-pCVTZ']
QZ = ['aug-cc-pVQZ', 'aug-cc-pCVQZ', 'd-aug-cc-pVQZ', 'd-aug-cc-pCVQZ']


class Tabler:
    def __init__(self, data: BasisMRADataCollection, cmap='BrBG', thresh=5, basis_list=None):

        if basis_list is None:
            self.basis_list = [
                'aug-cc-pVDZ', 'aug-cc-pVTZ', 'aug-cc-pVQZ',
                'aug-cc-pCVDZ', 'aug-cc-pCVTZ', 'aug-cc-pCVQZ',
                'd-aug-cc-pVDZ', 'd-aug-cc-pVTZ', 'd-aug-cc-pVQZ',
                'd-aug-cc-pCVDZ', 'd-aug-cc-pCVTZ', 'd-aug-cc-pCVQZ', ]
            self.basis_list = DZ + TZ + QZ
        else:
            self.basis_list = basis_list

        self.data = data
        self.cmap = cmap
        self.thresh = thresh

    def get_basis_data(self):
        df = self.data.iso_diff_data.query('omega==8')
        alphab = {}
        gammab = {}
        for b in self.basis_list:
            bdata = df.query('basis==@b')
            alphab[b] = bdata.set_index('molecule').alpha
            gammab[b] = bdata.set_index('molecule').gamma
        alpha_df = pd.DataFrame(alphab)
        gamma_df = pd.DataFrame(gammab)
        return alpha_df, gamma_df

    def get_basis_edata(self):
        df = self.data.energy_diff
        e = {}
        for b in self.basis_list:
            bdata = df.query('basis==@b')
            e[b] = bdata.set_index('molecule').error
        e_df = pd.DataFrame(e)
        return e_df

    def __style_full(self, df, fmt="{:.2e}"):
        first_row, second, flourine = partition_molecule_list(df.index)

        first = first_row + flourine
        one = df.query("molecule.isin(@first)")
        two = df.query("molecule.isin(@second)")

        vmax = df.iloc[:, :-1].max().max()

        # define a function to format the index values

        def format_index(idx):
            return r'\ce{' + '{}'.format(idx) + '}'

        one.index = one.index.map(format_index)
        two.index = two.index.map(format_index)

        bnorm = lambda x: background_with_norm(x, vmax)

        col1 = one.columns[:-1]
        style_one = one.style.apply(bnorm, subset=col1)
        style_one.format(fmt, subset=col1)
        style_one.format("{:.4g}", subset=one.columns[-1])
        style_one.applymap_index(
            lambda v: "rotatebox:{45}--rwrap--latex;", level=0, axis=1)

        col2 = two.columns[:-1]
        style_two = two.style.apply(bnorm, subset=col2)
        style_two.format(fmt, subset=col2)
        style_two.format("{:.4g}", subset=two.columns[-1])
        style_two.applymap_index(
            lambda v: "rotatebox:{45}--rwrap--latex;", level=0, axis=1)

        return style_one, style_two

    def __style_summary(self, df, fmt="{:.2e}"):
        summary = self.__get_summary(df)
        vmax = df.max().max()
        bnorm = lambda x: background_with_norm(x, vmax)
        # Define the maximum data value (in absolute terms) for normalization
        styled_df = summary.style.apply(bnorm)
        styled_df.format(fmt)
        return styled_df

    def get_styled_summary(self, iso_type):
        if iso_type == 'energy':
            return self.__style_summary(self.get_basis_edata())
        elif iso_type == 'alpha':
            return self.__style_summary(self.get_basis_data()[0])
        elif iso_type == 'gamma':
            return self.__style_summary(self.get_basis_data()[1])

    def get_styled_molecules(self, iso_type):
        if iso_type == 'energy':
            basis_e_data = self.get_basis_edata()
            mra_e_data = self.data.energy_df.query('basis=="mra"').set_index('molecule').error
            e_data = pd.concat([basis_e_data, mra_e_data.rename('MRA-ref')], axis=1)
            return self.__style_full(e_data)
        elif iso_type == 'alpha':
            basis_a_data = self.get_basis_data()[0]

            mra_a_data = self.data.iso_data.query('basis=="MRA" & omega==0').set_index(
                'molecule').alpha
            a_data = pd.concat([basis_a_data, mra_a_data.rename('MRA-ref')], axis=1)
            return self.__style_full(a_data)
        elif iso_type == 'gamma':
            basis_g_data = self.get_basis_data()[1]
            mra_g_data = self.data.iso_data.query('basis=="MRA" & omega==0').set_index(
                'molecule').gamma
            g_data = pd.concat([basis_g_data, mra_g_data], axis=1)
            return self.__style_full(g_data)

    def __get_summary(self, bdata):
        summary = bdata.describe().iloc[1:]
        new_idx = {}
        for i in summary.index:
            if i[-1] == '%':
                new_idx[i] = i[:-1] + '\\' + i[-1]
            else:
                new_idx[i] = i
        summary.rename(index=new_idx, inplace=True)
        return summary.T

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
