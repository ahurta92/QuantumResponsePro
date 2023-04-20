from .BasisMRADataAssembler import *
from pathlib import Path
import glob
import json


def get_polar_df(molecules, xc, op, database, basis):
    mra_data = get_mra_polar_data(molecules, xc, op, database)
    b_data = get_basis_polar_data(molecules, basis, xc, op, database)
    a_data = pd.concat([mra_data, b_data])
    a_data = a_data.reset_index(drop=True)
    return a_data


def save_compressed_polar(df, database_dir, data_file):
    compress_data_dir = database_dir + '/compress_data'
    file_name = compress_data_dir + "/" + data_file
    df.to_feather(file_name)


basis_sets = ['aug-cc-pVDZ', 'aug-cc-pVTZ', 'aug-cc-pVQZ',
              'aug-cc-pCVDZ', 'aug-cc-pCVTZ', 'aug-cc-pCVQZ',
              'd-aug-cc-pVDZ', 'd-aug-cc-pVTZ', 'd-aug-cc-pVQZ',
              'd-aug-cc-pCVDZ', 'd-aug-cc-pCVTZ', 'd-aug-cc-pCVQZ']

xc = 'hf'
op = 'dipole'

development = "/mnt/data/madness_data/post_acs/development"
august = "/mnt/data/madness_data/post_watoc/august"


# create a Path object with the path to the file


class BasisMRADataCollection:
    """
    ResponseDataBundle: A class designed to manage and organize a collection of DataFrames related to response properties for molecules. This class simplifies the process of comparing and analyzing data from various sources, such as MADNESS and Dalton quantum chemistry packages, by consolidating the response properties information into a single, easy-to-use structure.
    """

    def __init__(self, data_dir):
        self.data_dir = data_dir
        self.xc = 'hf'
        self.op = 'dipole'
        self.molecules = []

        all_data_path = self.data_dir.joinpath('all_polar_data.feather')

        if all_data_path.is_file():
            self.all_polar_data = pd.read_feather(all_data_path)
            self.molecules = list(self.all_polar_data.molecule.unique())
        else:
            for g in glob.glob('/*.mol', root_dir=self.data_dir.joinpath('molecules')):
                m = g.split('/')
                mol = m[-1].split('.')[0]
                self.molecules.append(mol)
            self.molecules = self.__report_convergence()[0]
            self.all_polar_data = get_polar_df(self.molecules, xc, op, august, basis_sets)
            self.all_polar_data.to_feather(all_data_path)

        iso_data_path: Path = Path(august + '/iso_polar_data.feather')

        if iso_data_path.is_file():
            self.iso_data = pd.read_feather(iso_data_path)
        else:
            self.iso_data = get_iso_data(self.all_polar_data)
            self.iso_data.to_feather(iso_data_path)
            self.iso_data[self.iso_data['gamma'].abs() < 1e-2] = 0

        iso_diff_path: Path = Path(august + '/iso_diff_data.feather')

        if iso_diff_path.is_file():
            self.iso_diff_data = pd.read_feather(iso_diff_path)
        else:
            self.iso_diff_data = create_iso_diff_df(self.iso_data)
            self.iso_diff_data.reset_index().to_feather(iso_diff_path)

        ij_df_path: Path = Path(august + '/ij_diff_data.feather')
        if ij_df_path.is_file():
            self.ij_diff = pd.read_feather(ij_df_path)
            self.ij_diff.reset_index(inplace=True)
        else:
            self.ij_diff = create_component_diff_df(self.all_polar_data)
            self.ij_diff.reset_index(inplace=True)
            self.ij_diff.to_feather(ij_df_path)

        energy_path: Path = Path(august + '/compress_data/energy_E.feather')
        if ij_df_path.is_file():
            self.energy_df = pd.read_feather(energy_path)
        else:
            print("Get energy df")
            pass

        self.detailed_iso_diff = make_detailed_df(self.iso_diff_data)
        self.detailed_ij_diff = make_detailed_df(self.ij_diff)
        self.detailed_energy_diff = make_detailed_df(self.energy_df)

    def __report_convergence(self):
        converged = []
        not_converged = []
        not_found = []
        type_error = []
        json_error = []
        for mol in self.molecules:
            try:
                check_mol = FrequencyData(mol, self.xc, self.op, self.database_dir)
                if check_mol.converged.all() and check_mol.converged.sum() == self.num_freq:
                    converged.append(mol)
                else:
                    not_converged.append(mol)
            except FileNotFoundError as f:
                not_found.append(mol)
            except TypeError as f:
                type_error.append(mol)
            except json.decoder.JSONDecodeError as j:
                json_error.append(mol)

        num_c = len(converged)
        num_n = len(not_converged)
        num_nf = len(not_found)
        num_json_e = len(json_error)
        num_type_e = len(type_error)
        total = num_c + num_n + num_nf + num_json_e + num_type_e
        not_converged = []
        part_converged = []
        if True:
            for mol in not_converged:
                check = FrequencyData(mol, self.xc, self.op, self.database_dir)
                if check.converged.any():
                    # print(mol,'\n',check.converged)
                    part_converged.append(mol)
                else:
                    not_converged.append(mol)
        num_not_converged = len(not_converged)
        num_part_converged = len(part_converged)
        print("converged : ", num_c)
        return converged, part_converged, not_converged, not_found, type_error, json_error


# selects outliers by change in percent error from the static to the
def select_basis_outliers(data, basis, thresh):
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


class HFDatabasePlots:
    def mol_iso_convergence(self, iso_data, mol, basis):
        width = 2 * 5
        height = 2 * 4
        g = plt.subplots(nrows=2, ncols=2, figsize=(width, height), constrained_layout=True,
                         sharey=False)
        fig = g[0]
        ax = g[1]
        data = get_frequency_compare(iso_data, mol, basis)
        title = [r'$\alpha(\omega)$', r'$\gamma(\omega)$']

        for i, axi in enumerate(ax):
            for j, axij in enumerate(axi):
                dij = data[i][j]
                axij.set_title(title[j])
                dij.plot(ax=axij, ls='dashdot')
                axij.set_xlabel(r'$\omega_i$')
        return g

    def mol_component_convergence(self, mol, ij_diff_detailed):
        ij = ['xx', 'yy', 'zz']
        mol_data = ij_diff_detailed.query(' molecule == @mol & ij==@ij')
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
                                                                                    "Percent Error").set_titles(
            "{col_name}-cc-pV(C)nZ").tight_layout(w_pad=0)
        g.fig.suptitle(mol)
        return g

    def plot_violin_swarm(self, data, iso_type, xdata, hdata, pals, ax):
        dotsize = 1.5

        if iso_type == 'alpha':
            ydata = data.alpha
        else:
            ydata = data.gamma
        color1 = pals[0]
        color2 = pals[1]

        p1 = sns.violinplot(x=xdata, y=ydata, ax=ax, split=False, scale="count", hue=hdata,
                            inner='quartile', cut=0, palette=color1, )

        p11 = sns.stripplot(x=data.valence, y=ydata, ax=ax, size=dotsize, hue=hdata,
                            dodge=True, palette=color2)

    def make_plots(self, datas, xdata, hdata, pals, outliers, axes):
        iso_types = ['alpha', 'gamma']
        titles = [r'$\alpha(\omega)$', r'$\gamma(\omega)$']

        for i, axi in enumerate(axes):
            out_mols = outliers[i]
            data = datas.query("not molecule.isin(@out_mols)")
            self.plot_violin_swarm(data, iso_types[i], xdata, hdata, pals, axi)
            axi.tick_params(axis='x', rotation=0)
            axi.grid(which="both")
            axi.tick_params(which="both", top="on", left="on", right="on", bottom="on", )
            axi.minorticks_on()
            axi.set_title(titles[i])
            axi.set_ylabel('Percentage Error')
            axi.set_xlabel('[n]')
            handles, labels = axi.get_legend_handles_labels()
            nlabels = ['D', 'T', 'Q']
            labels = nlabels
            # axi.legend(handles[:2], hdata.unique(), title=hdata.name)
            if i == 0:
                axi.set_ylabel('Percentage Error')
            else:
                axi.legend('', frameon=False)
                axi.set_ylabel(None)
            for l in axi.lines:
                l.set_linestyle('--')
                l.set_linewidth(1.2)
                l.set_color('green')
                l.set_alpha(0.8)
            for l in axi.lines[1::3]:
                l.set_linestyle('-')
                l.set_linewidth(1.5)
                l.set_color('black')
                l.set_alpha(0.8)
            axi.axhline(y=0, linewidth=1.2, ls="--", color="r")

    def make_plot(self, datas, xdata, hdata, pals, out_mols, iso_type, axes):
        iso_types = ['alpha', 'gamma']
        titles = [r'$\alpha(\omega)$', r'$\gamma(\omega)$']
        if iso_type == 'alpha':
            title = titles[0]
        else:
            title = titles[1]

        data = datas.query("not molecule.isin(@out_mols)")
        axi = axes
        self.plot_violin_swarm(data, iso_type, xdata, hdata, pals, axi)
        axi.tick_params(axis='x', rotation=0)
        axi.grid(which="both")
        axi.tick_params(which="both", top="on", left="on", right="on", bottom="on", )
        axi.minorticks_on()
        axi.set_title(title)
        axi.set_ylabel('Percentage Error')
        axi.set_xlabel('[n]')
        handles, labels = axi.get_legend_handles_labels()
        nlabels = ['D', 'T', 'Q']
        labels = nlabels
        axi.legend(handles[:4], hdata.unique(), title=hdata.name)
        axi.set_ylabel('Percentage Error')
        for l in axi.lines:
            l.set_linestyle('--')
            l.set_linewidth(1.2)
            l.set_color('green')
            l.set_alpha(0.8)
        for l in axi.lines[1::3]:
            l.set_linestyle('-')
            l.set_linewidth(1.5)
            l.set_color('black')
            l.set_alpha(0.8)
        axi.axhline(y=0, linewidth=1.2, ls="--", color="r")

    def freq_iso_plot(data, iso_diff_detailed, v_level, iso_type, mol_set="all", sharey=False,
                      thresh=.1):
        all_o = [0, 1, 2, 3, 4, 5, 6, 7, 8]
        omegas = all_o
        if mol_set == "all":
            data = iso_diff_detailed.query('omega.isin(@omegas) & valence==@v_level')
        else:
            data = iso_diff_detailed.query(
                'omega.isin(@omegas) & valence==@v_level & mol_system.isin(@mol_set)')

        mdata = data.query('valence==@v_level')
        sns.set(rc={"xtick.bottom": True, "ytick.left": True})
        aspect = 1.3
        if mol_set == 'all':
            g = sns.catplot(data=mdata,
                            x="omega",
                            row="polarization",
                            y=iso_type,
                            hue="mol_system",
                            col="augmentation",
                            kind='strip',
                            sharey=sharey,
                            dodge=True,
                            aspect=aspect
                            )
        else:
            aspect = 1.0
            set1 = select_basis_outliers(mdata, 'aug-cc-pV{}Z'.format(v_level), thresh)
            set2 = select_basis_outliers(mdata, 'd-aug-cc-pCV{}Z'.format(v_level), thresh * .1)
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
            set_titles("{col_name}-cc-p{row_name}" + "{}Z".format(v_level)).tight_layout(w_pad=0)
        g.fig.suptitle(iso_type)
        return g
