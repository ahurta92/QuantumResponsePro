import math
from pathlib import Path
import pandas as pd
import matplotlib as mpl
import os
import numpy as np
from quantumresponsepro.BasisMRADataAssembler import *
import plotly.graph_objects as go
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from quantumresponsepro.BasisMRADataCollection import get_quad_df, get_polar_df


def round_to_n_significant_figures(x, n=3):
    if x == 0:
        return 0
    else:
        return round(x, -int(np.floor(np.log10(abs(x)))) + (n - 1))


class QuadraticDatabase:
    def __init__(self, mols, basis_sets, xc, op, freq, database, overwrite=False):
        self.q_df = None
        self.basis_set_error_df = None
        self.freq = freq
        self.molecules = mols
        self.basis_sets = basis_sets
        self.xc = xc
        self.op = op
        self.database_path = database
        self.overwrite = overwrite

        self.initialize_dfs()

    def initialize_dfs(self):
        data_frame_attributes = ['q_df', 'basis_set_error_df']

        for attr_name in data_frame_attributes:
            filename = self.database_path.joinpath(f"{attr_name}.csv")

            if not self.overwrite and os.path.exists(filename):
                df = pd.read_csv(filename)
            else:
                df = self.generate_df(attr_name)

            setattr(self, attr_name, df)
        self.__make_all_dfs_detailed()

    def generate_df(self, attr_name):
        # Replace this section with your specific data generating code
        if attr_name == 'q_df':
            return self.__generate_quad_df()
        elif attr_name == 'basis_set_error_df':
            # Implement your data generating logic here
            return self.__generate_basis_error_df()

    def save_dfs(self):
        data_frame_attributes = ['q_df', 'basis_set_error_df']

        for attr_name in data_frame_attributes:
            df = getattr(self, attr_name)
            filename = self.database_path.joinpath(f"{attr_name}.csv")
            df.to_csv(filename, index=False)

    def __make_all_dfs_detailed(self):
        data_frame_attributes = ['q_df', 'basis_set_error_df']

        for attr_name in data_frame_attributes:
            df = getattr(self, attr_name)  # Access the DataFrame attribute
            print(df)
            modified_df = make_detailed_df(df)  # Apply the function
            setattr(self, attr_name, modified_df)  #

    def __generate_quad_df(self):
        q_df = get_quad_df(self.molecules, self.xc, self.op, self.database_path, self.basis_sets)
        q_df_i = q_df.reset_index()
        # truncate the Afreq Bfreq and Cfreq to 3 decimal places

        q_df_i['Afreq'] = q_df_i['Afreq'].apply(lambda x: round(x, 3))
        q_df_i['Bfreq'] = q_df_i['Bfreq'].apply(lambda x: round(x, 3))
        q_df_i['Cfreq'] = q_df_i['Cfreq'].apply(lambda x: round(x, 3))

        # Function to round to n significant figures

        # Apply the function to the float column
        q_df_i['Beta'] = q_df_i['Beta'].apply(
            lambda x: round_to_n_significant_figures(x, 3))

        # q_df_i['Beta'] = q_df_i['Beta'].apply(lambda x: round(x, 2))

        return q_df_i.copy()

    def __generate_basis_error_df(self):

        basis_error_df = pd.DataFrame()
        for mol in self.molecules:
            mol_data, afreq = self.__get_mol_and_freq_data(mol)
            for basis in self.basis_sets:
                for b in range(len(self.freq)):
                    bf = afreq[b]
                    for c in range(b, len(self.freq)):
                        cf = afreq[c]
                        try:
                            df_i = self.__get_quad_basis_error_data_df(mol, basis, mol_data, bf, cf)
                            basis_error_df = pd.concat([basis_error_df, df_i])
                        except IndexError as e:
                            print(e)
                            print(mol, basis)
                            continue
        basis_error_df.dropna(inplace=True)
        return basis_error_df.copy()

    def __get_mol_and_freq_data(self, mol):
        mol_data = self.q_df.query('molecule == @mol').copy()
        freq = list(mol_data['Afreq'].unique())
        return mol_data, freq

    def __get_molecules_basis_error_df(self, mol, basis, Bfreq_i, Cfreq_j):
        mol_data, freq = self.__get_mol_and_freq_data(mol)
        Bfreq = freq[Bfreq_i]
        Cfreq = freq[Cfreq_j]
        b_data = self.__get_quad_basis_error_data_df(mol, basis, mol_data, Bfreq,
                                                     Cfreq)
        return b_data

    def __get_basis_df(self, mol, basis, Bfreq_i, Cfreq_j):
        mol_data, freq = self.__get_mol_and_freq_data(mol)
        Bfreq = freq[Bfreq_i]
        Cfreq = freq[Cfreq_j]
        mra_data = self.__get_basis_freq_df(basis, mol_data, Bfreq, Cfreq)
        return mra_data

    def __get_dataset_df(self, mol, basis_sets, Bfreq_i, Cfreq_j):
        df = pd.DataFrame()
        for basis in basis_sets:
            b_data = self.__get_molecules_basis_error_df(mol, basis, Bfreq_i, Cfreq_j)
            df = pd.concat([df, b_data])

        return df

    def __get_molecule_frequency_data(self, mol, basis_sets, Bfreqs, Cfreqs):
        '''
        returns beta dataframe at a given frequency for a set of basis sets
        '''
        df = pd.DataFrame()

        for basis in basis_sets:
            for b in Bfreqs:
                for c in Cfreqs:
                    df = pd.concat([df, self.__get_basis_df(mol, basis, b, c)])

        return df

    def __get_quad_basis_error_data_df(self, mol, basis, mol_data, Bfreq, Cfreq):
        """
        This function returns
        """
        b_data = mol_data.query('basis == @basis & Bfreq == @Bfreq & Cfreq== @Cfreq').copy()

        # remove the nan values
        b_data.dropna(inplace=True)
        b_data.drop_duplicates(subset=['ijk'], inplace=True)
        # remove Beta less than 1e-4
        m_data = mol_data.query('basis == "MRA" & Bfreq==@Bfreq & Cfreq == @Cfreq ').copy()
        # remove the Beta less than .001
        m_data.dropna(inplace=True)
        # remove the Beta less than 1e-4
        m_data = m_data.query('Beta.abs() > .005').copy()
        BC_index = ['Afreq', 'Bfreq', 'Cfreq', 'ijk']

        m_data.drop_duplicates(inplace=True, subset=['ijk'])
        m_data.set_index(BC_index, inplace=True)
        b_data.set_index(BC_index, inplace=True)

        # b_data.drop_duplicates(inplace=True, subset=['ijk'])
        b_data.drop_duplicates()

        b_data.drop_duplicates()
        # drop rows with Beta less than 1e-4
        b_data = b_data.query('Beta.abs() > .005').copy()
        b_data['Beta_MRA'] = m_data['Beta']
        b_data['Beta_sign'] = np.sign(b_data['Beta'] * b_data['Beta_MRA'])
        b_data['Beta_MRA_signed'] = b_data['Beta_MRA'] * b_data['Beta_sign']
        b_data['Beta_diff'] = (b_data['Beta'] - b_data['Beta_MRA']) / b_data['Beta_MRA'] * 100

        b_data['Beta'] = b_data['Beta_diff']
        b_data.drop(columns=['Beta_MRA', 'Beta_sign', 'Beta_MRA_signed', 'Beta_diff'], inplace=True)

        return b_data

    def __get_basis_freq_df(self, basis, mol_data, Bfreq, Cfreq):
        """
        This function returns
        """
        b_data = mol_data.query('basis == @basis & Bfreq == @Bfreq & Cfreq== @Cfreq').copy()

        # remove the nan values
        b_data.dropna(inplace=True)
        b_data.drop_duplicates(subset=['ijk'], inplace=True)
        # remove Beta less than 1e-4
        b_data.drop_duplicates()
        # drop rows with Beta less than 1e-4
        b_data = b_data.query('Beta.abs() > .005').copy()

        return b_data




class PolarizabilityData:

    def __init__(self, molecules, xc, op, basis_sets, database, overwrite=False):
        self.energy_diff_df = None
        self.energy_df = None
        self.iso_data = None
        self.eigen_diff = None
        self.alpha_eigen = None
        self.polar_data = None
        self.molecules = molecules
        self.xc = xc
        self.op = op
        self.basis_sets = basis_sets
        self.database = database
        self.overwrite = overwrite

        self.initialize_dfs()

    def initialize_dfs(self):
        data_frame_attributes = ['energy_df', 'energy_diff_df', 'polar_data', 'alpha_eigen',
                                 'eigen_diff', 'iso_data']

        for attr_name in data_frame_attributes:
            filename = self.database.joinpath(f"{attr_name}.csv")
            if not self.overwrite and os.path.exists(filename):
                df = pd.read_csv(filename)
            else:
                df = self.generate_df(attr_name)

            setattr(self, attr_name, df)
        self.__make_all_dfs_detailed()

    def generate_df(self, attr_name):
        if attr_name == 'energy_df':
            return get_energy_data(self.molecules, self.xc, self.op, self.basis_sets, self.database)
        elif attr_name == 'energy_diff_df':
            return get_energy_diff_data(mols, self.xc, self.op, basis_sets,
                                        polarizability_database)
        elif attr_name == 'polar_data':
            return get_polar_df(mols, self.xc, self.op, polarizability_database, basis_sets)
        elif attr_name == 'alpha_eigen':
            return get_ij_eigen(self.polar_data)
        elif attr_name == 'eigen_diff':
            return create_component_diff_df(self.alpha_eigen)
        elif attr_name == 'iso_data':
            return get_iso_data(self.polar_data)

        # Add similar cases for other DataFrame attributes

    def save_dfs(self):
        data_frame_attributes = ['energy_df', 'energy_diff_df', 'polar_data', 'alpha_eigen',
                                 'eigen_diff', 'iso_data']

        for attr_name in data_frame_attributes:
            df = getattr(self, attr_name)
            filename = self.database.joinpath(f"{attr_name}.csv")
            df.to_csv(filename, index=False)

    def __make_all_dfs_detailed(self):
        data_frame_attributes = ['energy_df', 'energy_diff_df', 'polar_data', 'alpha_eigen',
                                 'eigen_diff', 'iso_data']

        for attr_name in data_frame_attributes:
            df = getattr(self, attr_name)  # Access the DataFrame attribute
            print(df)
            modified_df = make_detailed_df(df)  # Apply the function
            setattr(self, attr_name, modified_df)  #




def plot_component_diff(data, valence, component, ax, color_pal=None):
    data['valence'] = pd.Categorical(data['valence'], valence)
    # make the polarization order ['V','CV']
    # make the polarization order ['V','CV']
    try:
        data['polarization'] = pd.Categorical(data['polarization'], ['V', 'CV'])
    except ValueError as v:
        print(v)
        print(data)
        pass
    # font scale *2
    g = sns.lineplot(data=data, x='valence', y='diff', hue='Type',
                     style=component,
                     markers=True, ax=ax, markersize=10, legend='full', dashes=False,
                     palette=color_pal)

    handles, labels = g.get_legend_handles_labels()
    N = math.ceil(len(handles) / 2)
    print(N)
    filtered_handles = handles[N:]
    filtered_labels = labels[N:]
    g.legend(filtered_handles, filtered_labels)

    ax.get_xaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
    ax.get_yaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
    ax.minorticks_on()
    ax.grid(which='major', color='w', linewidth=1.0)
    ax.grid(which='minor', color='w', linewidth=0.5)
    ax.set_xlabel('Valence')
    ax.set_ylabel('Percent Error')
    ax.axhline(y=0, color='k')

    # for the line plots, only keep the "style" legend not the hue legend
    handles, labels = ax.get_legend_handles_labels()

    # add the legends for the components only

    return g


def plot_energy_diff(data, valence, ax, color_pal=None):
    data['valence'] = pd.Categorical(data['valence'], valence)
    # make the polarization order ['V','CV']
    try:
        data['polarization'] = pd.Categorical(data['polarization'], ['V', 'CV'])
    except ValueError as v:
        print(v)
        print(data)
        pass

    # show minor ticks and grid lines  using with
    g = sns.lineplot(data=data, x='valence', y='error', hue='Type',
                     markers=True, ax=ax, palette=color_pal)

    ax.get_xaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
    ax.get_yaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
    ax.minorticks_on()
    ax.grid(which='major', color='w', linewidth=1.0)
    ax.grid(which='minor', color='w', linewidth=0.5)

    ax.set_xlabel('Valence')
    ax.set_ylabel('Error (a.u.)')
    ax.axhline(y=0, color='k')

    return g


pal = sns.color_palette("seismic_r", 4).as_hex()
p1 = [pal[1], pal[0], pal[2], pal[3]]
pal = sns.color_palette(p1)

pal2 = sns.color_palette("coolwarm_r", 4).as_hex()
p2 = [pal2[1], pal2[0], pal2[2], pal2[3]]
p3 = [pal2[1], pal2[2], ]
light_pal = sns.color_palette(p2)
simple_pal = sns.color_palette(p3)

sns.set_theme('paper', 'darkgrid', palette=light_pal, font='sans-serif', font_scale=1.0)


class basis_set_analaysis:

    def __init__(self, qda: QuadraticDatabase, pda: PolarizabilityData):
        self.qda = qda
        self.pda = pda

    def basis_set_analysis_plot(self, molecule, valence, type=None, ):

        r_plot, axes = plt.subplots(1, 3, figsize=(9, 4), frameon=True, layout='constrained')
        if type is None:
            pal = light_pal
        else:
            pal = simple_pal

        e_diff = self.pda.energy_diff_df.query('molecule == @molecule')
        a_diff = self.pda.eigen_diff.query('molecule == @molecule')
        # rename alpha to diff
        a_diff.rename(columns={'alpha': 'diff'}, inplace=True)
        q_diff = self.qda.basis_set_error_df.query('molecule == @molecule')
        q_diff.rename(columns={'Beta': 'diff'}, inplace=True)

        plot_energy_diff(e_diff, valence, axes[0], color_pal=pal)
        alpha_plot = plot_component_diff(a_diff, valence, 'ij', axes[1],
                                         color_pal=pal)
        beta_plot = plot_component_diff(q_diff, valence, 'ijk', axes[2],
                                        color_pal=pal)

        axes[0].set_title('Energy')

        axes[1].set_title(r'$\alpha_{ij}(0;0,0)$')
        axes[1].set_yscale('symlog', linthresh=1e-2, linscale=0.25)

        axes[2].set_title(r'$\beta_{ijk}(0;0,0)$')
        # axes[2].set_yscale('symlog', linthresh=1e-2, linscale=0.25)

        return r_plot, axes

database_path = Path('/mnt/data/madness_data/august_no_symmetry')
# get the unique frequencies
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)

basis_sets = ['aug-cc-pVDZ', 'aug-cc-pVTZ', 'aug-cc-pVQZ', 'd-aug-cc-pVDZ', 'd-aug-cc-pVTZ',
              'd-aug-cc-pVQZ'] + ['aug-cc-pCVDZ', 'aug-cc-pCVTZ', 'aug-cc-pCVQZ', 'd-aug-cc-pCVDZ',
                                  'd-aug-cc-pCVTZ', 'd-aug-cc-pCVQZ']

mols = ['CH3SH', 'H2O', 'CH3OH', 'C2H2', 'C2H4', 'CH3F', 'CH3OH','NaLi' ]
mol = 'H2O'
xc = 'hf'
op = 'dipole'

freq = [0, 8]
rdb = QuadraticDatabase(mols, basis_sets, xc, op, freq, database_path)

polarizability_database = Path('/mnt/data/madness_data/post_watoc/august')

polar_data = PolarizabilityData(mols, 'hf', 'dipole', basis_sets=basis_sets,
                                database=polarizability_database,
                                overwrite=True)
print(polar_data.eigen_diff)
mol = 'NaLi'

b_plotter = basis_set_analaysis(rdb, polar_data)
r_plot, axes = b_plotter.basis_set_analysis_plot(mol, ['D', 'T', 'Q', ], type=None)
r_plot.show()
# r_plot.savefig(paper_path.joinpath("nacl_plot.png"), dpi=1000)
