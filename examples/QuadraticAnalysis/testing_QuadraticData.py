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
        self.vector_q_df = None
        self.vector_basis_set_error_df = None
        self.freq = freq
        self.molecules = mols
        self.basis_sets = basis_sets
        self.xc = xc
        self.op = op
        self.database_path = database
        self.overwrite = overwrite

        self.initialize_dfs()

    def initialize_dfs(self):
        data_frame_attributes = ['q_df', 'basis_set_error_df', 'vector_q_df',
                                 'vector_basis_set_error_df']

        for attr_name in data_frame_attributes:
            filename = self.database_path.joinpath(f"{attr_name}.csv")

            if not self.overwrite and os.path.exists(filename):
                df = pd.read_csv(filename)
            else:
                df = self.generate_df(attr_name)

            setattr(self, attr_name, df)

    def generate_df(self, attr_name):
        # Replace this section with your specific data generating code
        if attr_name == 'q_df':
            return self.__generate_quad_df()
        elif attr_name == 'basis_set_error_df':
            # Implement your data generating logic here
            return self.__generate_basis_error_df()
        elif attr_name == 'vector_q_df':
            return self.__generate_vector_q_df()
        elif attr_name == 'vector_basis_set_error_df':
            return self.__generate_vector_basis_error_df()

    def save_dfs(self):
        data_frame_attributes = ['q_df', 'basis_set_error_df', 'vector_q_df',
                                 'vector_basis_set_error_df']

        for attr_name in data_frame_attributes:
            df = getattr(self, attr_name).copy()
            filename = self.database_path.joinpath(f"{attr_name}.csv")
            df.to_csv(filename, index=False)

    def make_all_dfs_detailed(self):
        data_frame_attributes = ['q_df', 'basis_set_error_df', 'vector_q_df']

        for attr_name in data_frame_attributes:
            df = getattr(self, attr_name)  # Access the DataFrame attribute
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
                            df_i = self.__get_quad_basis_error_data_df(basis, mol_data, bf, cf)
                            basis_error_df = pd.concat([basis_error_df, df_i])
                        except IndexError as e:
                            print(e)
                            print(mol, basis)
                            continue
        basis_error_df.dropna(inplace=True)
        basis_error_df.reset_index(inplace=True)
        return basis_error_df.copy()

    def __generate_vector_q_df(self):

        vector_df = []
        for mol in self.molecules:
            mol_data, afreq = self.__get_mol_and_freq_data(mol)
            for basis in self.basis_sets + ['MRA']:
                for b in range(len(self.freq)):
                    bf = afreq[b]
                    for c in range(b, len(self.freq)):
                        cf = afreq[c]
                        try:
                            df_i = self.__get_vector_representation(basis, mol_data, bf, cf)
                            print(df_i)
                            vector_df.append(df_i)
                        except IndexError as e:
                            print(e)
                            print(mol, basis)
                            continue
        vector_df = pd.concat(vector_df, axis=1).T
        print(vector_df)
        vector_df.dropna(inplace=True)
        vector_df.reset_index(inplace=True)
        return vector_df.copy()

    def __generate_vector_basis_error_df(self):

        basis_error_df = pd.DataFrame()
        for mol in self.molecules:
            mol_data, afreq = self.__get_vector_mol_and_freq_data(mol)
            for basis in self.basis_sets:
                for b in range(len(self.freq)):
                    bf = afreq[b]
                    for c in range(b, len(self.freq)):
                        cf = afreq[c]
                        try:
                            df_i = self.__get_quad_vector_basis_error_data_df(basis, mol_data, bf,
                                                                              cf)
                            basis_error_df = pd.concat([basis_error_df, df_i])
                        except IndexError as e:
                            print(e)
                            print(mol, basis)
                            continue
                        except TypeError as t:
                            print(t)
                            print(mol, basis, bf, cf)
                            continue

        basis_error_df.dropna(inplace=True)
        basis_error_df.reset_index(inplace=True)
        return basis_error_df.copy()

    def __get_mol_and_freq_data(self, mol):
        mol_data = self.q_df.query('molecule == @mol').copy()
        freq = list(mol_data['Afreq'].unique())
        return mol_data, freq

    def __get_vector_mol_and_freq_data(self, mol):
        mol_data = self.vector_q_df.query('molecule == @mol').copy()
        freq = list(mol_data['Afreq'].unique())
        return mol_data, freq

    def __get_molecules_basis_error_df(self, mol, basis, Bfreq_i, Cfreq_j):
        mol_data, freq = self.__get_mol_and_freq_data(mol)
        Bfreq = freq[Bfreq_i]
        Cfreq = freq[Cfreq_j]
        b_data = self.__get_quad_basis_error_data_df(basis, mol_data, Bfreq,
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

    def __get_quad_basis_error_data_df(self, basis, mol_data, Bfreq, Cfreq):
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
        b_data['Beta_diff'] = (b_data['Beta'] - b_data['Beta_MRA'])  # / b_data['Beta_MRA'] * 100
        print(b_data)

        b_data['Beta'] = b_data['Beta_diff']
        b_data.drop(columns=['Beta_MRA', 'Beta_sign', 'Beta_MRA_signed', 'Beta_diff'], inplace=True)

        return b_data

    def __get_quad_vector_basis_error_data_df(self, basis, mol_data, Bfreq, Cfreq):
        """
        This function returns
        """
        b_data = mol_data.query('basis == @basis & Bfreq == @Bfreq & Cfreq== @Cfreq').copy()

        # remove the nan values
        b_data.dropna(inplace=True)
        b_data.drop_duplicates(subset=['component'], inplace=True)
        # remove Beta less than 1e-4
        m_data = mol_data.query('basis == "MRA" & Bfreq==@Bfreq & Cfreq == @Cfreq ').copy()
        # remove the Beta less than .001
        m_data.dropna(inplace=True)
        # remove the Beta less than 1e-4
        m_data = m_data.query('Beta.abs() > .005').copy()
        BC_index = ['Bfreq', 'Cfreq', 'component']

        m_data.drop_duplicates(inplace=True, subset=['component'])
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
        b_data['Beta_diff'] = (b_data['Beta'] - b_data['Beta_MRA'])  # / b_data['Beta_MRA'] * 100

        print(b_data)
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

    def __get_vector_representation(self, basis, mol_data, Bfreq, Cfreq):
        """
        This function returns
        """

        b_data = mol_data.query('basis == @basis & Bfreq == @Bfreq & Cfreq== @Cfreq').copy()
        molecule = b_data['molecule'].unique()[0]
        b_data.drop_duplicates(inplace=True, subset=['ijk'])
        b_data.set_index('ijk', inplace=True)
        # the vector representation of the hyperpolarizability tensor accumulates components
        # whose output are in the same direction.
        # for example ax=beta(xxx)+beta(xyy)+beta(xzz) and ay=beta(yxx)+beta(yyy)+beta(yzz) and
        # az=beta(zxx)+beta(zyy)+beta(zzz)
        # the vector representation is then [ax,ay,az]

        ax = 0
        ay = 0
        az = 0
        ax_indices = ['XXX', 'XYY', 'XZZ']
        ay_indices = ['YXX', 'YYY', 'YZZ']
        az_indices = ['ZXX', 'ZYY', 'ZZZ']

        for i in range(3):
            print(i)
            try:
                f = b_data.loc[ax_indices[i]].Beta
                print(f)
                ax += f
            except KeyError as e:
                print(e)
                print('failed')
        for i in range(3):
            try:
                f = b_data.loc[ay_indices[i]].Beta
                ay += f
            except KeyError as e:
                print(e)
        for i in range(3):
            try:
                f = b_data.loc[az_indices[i]].Beta
                print(f)
                az += f
            except KeyError as e:
                print(e)

        index = ['molecule', 'basis', 'Afreq', 'Bfreq', 'Cfreq', 'component', 'Beta']
        Afreq = Bfreq + Cfreq
        # create the pd.Series
        vectorx = pd.Series([molecule, basis, Afreq, Bfreq, Cfreq, 'x', ax], index=index)
        vectory = pd.Series([molecule, basis, Afreq, Bfreq, Cfreq, 'y', ay], index=index)
        vectorz = pd.Series([molecule, basis, Afreq, Bfreq, Cfreq, 'z', az], index=index)

        vector = pd.concat([vectorx, vectory, vectorz], axis=1)

        return vector


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

    def generate_df(self, attr_name):
        if attr_name == 'energy_df':
            return get_energy_data(self.molecules, self.xc, self.op, self.basis_sets, self.database)
        elif attr_name == 'energy_diff_df':
            return get_energy_diff_data(mols, self.xc, self.op, basis_sets,
                                        self.database)
        elif attr_name == 'polar_data':
            return get_polar_df(mols, self.xc, self.op, self.database, basis_sets)
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
            df = getattr(self, attr_name).copy()
            filename = self.database.joinpath(f"{attr_name}.csv")
            df.to_csv(filename, index=False)

    def make_all_dfs_detailed(self):
        data_frame_attributes = ['energy_df', 'energy_diff_df', 'polar_data', 'alpha_eigen',
                                 'eigen_diff', 'iso_data']

        for attr_name in data_frame_attributes:
            df = getattr(self, attr_name)  # Access the DataFrame attribute
            print(df)
            modified_df = make_detailed_df(df)  # Apply the function
            setattr(self, attr_name, modified_df)  #


def plot_component_diff(data, valence, component, ax, color_pal=None):
    data = data.copy()
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
        self.qda.make_all_dfs_detailed()
        self.pda = pda
        self.pda.make_all_dfs_detailed()

    def basis_set_analysis_plot(self, molecule, valence, omega_1=0.0, omega_2=0.0, type=None, ):

        r_plot, axes = plt.subplots(1, 3, figsize=(9, 4), frameon=True, layout='constrained')

        if type is None:
            pal = light_pal
        else:
            pal = simple_pal

        e_diff = self.pda.energy_diff_df.query('molecule == @molecule').copy()
        a_diff = self.pda.eigen_diff.query('molecule == @molecule & omega==@omega_1').copy()
        # rename alpha to diff
        a_diff.rename(columns={'alpha': 'diff'}, inplace=True)
        q_diff = self.qda.vector_basis_set_error_df.query('molecule == @molecule & '
                                                          'Bfreq==@omega_1 & '
                                                          'Cfreq==@omega_2').copy()
        q_diff.rename(columns={'Beta': 'diff'}, inplace=True)

        plot_energy_diff(e_diff, valence, axes[0], color_pal=pal)
        alpha_plot = plot_component_diff(a_diff, valence, 'ij', axes[1],
                                         color_pal=pal)

        print(q_diff)
        beta_plot = plot_component_diff(q_diff, valence, 'component', axes[2],
                                        color_pal=pal)

        axes[0].set_title('Energy')

        axes[1].set_title(r'$\alpha_{ij}(0;0,0)$')
        axes[1].set_yscale('symlog', linthresh=1e-2, linscale=0.25)

        axes[2].set_title(r'$\beta_{ijk}(0;0,0)$')
        # axes[2].set_yscale('symlog', linthresh=1e-2, linscale=0.25)

        return r_plot, axes


class QuadVisualization:

    def __init__(self, data: QuadraticDatabase):
        self.quad_data = data

    def unit_sphere_representation(self, data, basis, omega_1, omega_2, radius=1, num_points=1000,
                                   basis_error=False, colormap='Blues', sizeref=0.5):

        data = self.quad_data.q_df.query('basis==@basis & Afreq==@omega_1 & Bfreq==@omega_2')
        # Initialize results
        results = []
        # Generate 1000 points on the unit sphere
        points = self.generate_points_on_sphere(radius, num_points)

        beta_tensor = self.beta_df_np(data)
        # Compute projection
        for p in points:
            bi = self.beta_proj(beta_tensor, p)
            results.append(bi)
        results = np.array(results)

        p1 = points.copy()
        p2 = points.copy() + results

        x = p1[:, 0]
        y = p1[:, 1]
        z = p1[:, 2]
        u = p2[:, 0]
        v = p2[:, 1]
        w = p2[:, 2]

        quiver = go.Cone(x=x, y=y, z=z, u=u, v=v, w=w, colorscale=colormap, sizemode="scaled",
                         sizeref=sizeref, anchor='tail')
        # return the quiver
        return quiver

    def unit_sphere_representation_basis(self, mol, basis_sets, omega_1, omega_2, radius=1,
                                         num_points=1000,
                                         colormap='Blues', sizeref=0.5, shift=1.0,
                                         sizemode='scaled', basis_error=False):

        if not basis_error:

            data = self.quad_data.q_df.query(
                'molecule==@mol & Cfreq==@omega_1 & Bfreq==@omega_2 ')
        else:
            data = self.quad_data.basis_set_error_df.query(
                'molecule==@mol & Cfreq==@omega_1 & Bfreq==@omega_2 ')

        x = np.zeros((len(basis_sets) * num_points))
        y = np.zeros((len(basis_sets) * num_points))
        z = np.zeros((len(basis_sets) * num_points))
        u = np.zeros((len(basis_sets) * num_points))
        v = np.zeros((len(basis_sets) * num_points))
        w = np.zeros((len(basis_sets) * num_points))

        text = []

        for i, basis in enumerate(basis_sets):
            for j in range(num_points):
                text.append(basis)

        # minimal_index_ijk = ['XXX', 'YYY', 'ZZZ', 'XYZ', 'XXY', 'YXY', 'ZXX', 'ZXZ', 'YYZ', 'YZZ']

        # for each basis set generate the xyzuvw data
        for i, basis in enumerate(basis_sets):

            bdata = data.query('basis==@basis & Afreq==@omega_1 & Bfreq==@omega_2')
            bdata.reset_index(inplace=True)

            bdata.set_index('ijk', inplace=True)
            bdata.sort_index(inplace=True)
            # bdata = bdata.loc[minimal_index_ijk]
            # Initialize results
            results = []
            # Generate 1000 points on the unit sphere
            points = self.generate_points_on_sphere(radius, num_points)

            beta_tensor = self.beta_df_np(bdata)
            # Compute projection
            for p in points:
                bi = self.beta_proj(beta_tensor, p)
                results.append(bi)
            results = np.array(results)
            p1 = points.copy()
            p2 = points.copy() + results
            x[i * num_points:(i + 1) * num_points] = p1[:, 0] + shift * (i - len(basis_sets) / 2)
            y[i * num_points:(i + 1) * num_points] = p1[:, 1]
            z[i * num_points:(i + 1) * num_points] = p1[:, 2]

            u[i * num_points:(i + 1) * num_points] = p2[:, 0] + shift * (i - len(basis_sets) / 2)
            v[i * num_points:(i + 1) * num_points] = p2[:, 1]
            w[i * num_points:(i + 1) * num_points] = p2[:, 2]

            # also generate text basis set name data to add to quiver
            # append the data to the quiver

        quiver = go.Cone(x=x, y=y, z=z, u=u, v=v, w=w, colorscale=colormap, sizemode=sizemode,
                         sizeref=sizeref, text=text, opacity=1.0)
        # get the vector representation of beta

        # return the quiver
        return quiver

    def vector_representation_basis(self, mol, basis_sets, omega_1, omega_2, colormap='Greens',
                                    sizeref=0.5, shift=1.0, sizemode='scaled'):

        data = self.quad_data.vector_q_df.query('molecule==@mol').copy()
        print(data)
        x = np.zeros(len(basis_sets))
        y = np.zeros(len(basis_sets))
        z = np.zeros(len(basis_sets))
        u = np.zeros(len(basis_sets))
        v = np.zeros(len(basis_sets))
        w = np.zeros(len(basis_sets))
        text = []

        for i, basis in enumerate(basis_sets):
            text.append(basis)

        for i, basis in enumerate(basis_sets):
            bdata = data.query('basis==@basis & Bfreq==@omega_1 & Cfreq==@omega_2')
            av = bdata.copy()
            av.set_index('component', inplace=True)
            av.drop_duplicates(inplace=True)

            x[i] = 0 + shift * (i - len(basis_sets) / 2)
            y[i] = 0 - shift
            z[i] = 0
            u[i] = av.loc['x'].Beta + shift * (i - len(basis_sets) / 2)
            v[i] = av.loc['y'].Beta - shift
            w[i] = av.loc['z'].Beta

        quiver2 = go.Cone(x=x, y=y, z=z, u=u, v=v, w=w, colorscale='Greens', sizemode=sizemode,
                          sizeref=sizeref, text=text, opacity=1.0)
        return quiver2

    def unit_sphere_representation_basis_error(self, mol, basis_sets, omega_1, omega_2,
                                               radius=1,
                                               num_points=1000,
                                               colormap='Blues', sizeref=0.5, shift=1.0,
                                               sizemode='scaled'):

        basis_error_df = self.quad_data.basis_set_error_df

        x = np.zeros((len(basis_sets) * num_points))
        y = np.zeros((len(basis_sets) * num_points))
        z = np.zeros((len(basis_sets) * num_points))
        u = np.zeros((len(basis_sets) * num_points))
        v = np.zeros((len(basis_sets) * num_points))
        w = np.zeros((len(basis_sets) * num_points))

        text = []

        for i, basis in enumerate(basis_sets):
            for j in range(num_points):
                text.append(basis)

        # minimal_index_ijk = ['XXX', 'YYY', 'ZZZ', 'XYZ', 'XXY', 'YXY', 'ZXX', 'ZXZ', 'YYZ', 'YZZ']

        # for each basis set generate the xyzuvw data
        for i, basis in enumerate(basis_sets):
            bdata = basis_error_df.query('basis==@basis & Afreq==@omega_1 & Bfreq==@omega_2')
            bdata.reset_index(inplace=True)

            bdata.set_index('ijk', inplace=True)
            bdata.sort_index(inplace=True)
            # bdata = bdata.loc[minimal_index_ijk]
            # Initialize results
            results = []
            # Generate 1000 points on the unit sphere
            points = self.generate_points_on_sphere(radius, num_points)

            beta_tensor = self.beta_df_np(bdata)
            # Compute projection
            for p in points:
                bi = self.beta_proj(beta_tensor, p)
                results.append(bi)
            results = np.array(results)
            p1 = points.copy()
            p2 = points.copy() + results
            x[i * num_points:(i + 1) * num_points] = p1[:, 0] + shift * (i - len(basis_sets) / 2)
            y[i * num_points:(i + 1) * num_points] = p1[:, 1]
            z[i * num_points:(i + 1) * num_points] = p1[:, 2]

            u[i * num_points:(i + 1) * num_points] = p2[:, 0] + shift * (i - len(basis_sets) / 2)
            v[i * num_points:(i + 1) * num_points] = p2[:, 1]
            w[i * num_points:(i + 1) * num_points] = p2[:, 2]

            # also generate text basis set name data to add to quiver
            # append the data to the quiver

        quiver = go.Cone(x=x, y=y, z=z, u=u, v=v, w=w, colorscale=colormap, sizemode=sizemode,
                         sizeref=sizeref, text=text, opacity=1.0, anchor='tail')
        # return the quiver
        return quiver

    def generate_points_on_sphere(self, radius, num_points):
        points = np.zeros((num_points, 3))
        phi = math.pi * (3. - math.sqrt(5.))  # Golden angle

        for i in range(num_points):
            y = radius - (i / float(num_points - 1)) * 2 * radius  # y goes from 1 to -1
            ri = math.sqrt(radius ** 2 - y * y)  # radius at y

            theta = phi * i  # Golden angle increment

            x = math.cos(theta) * ri
            z = math.sin(theta) * ri

            points[i, :] = [x, y, z]

        return points

    def beta_df_np(self, beta_df):
        beta_df = beta_df.copy()
        print(beta_df)
        xyz_to_012 = {'X': 0, 'Y': 1, 'Z': 2}
        beta_tensor = np.zeros((3, 3, 3))
        beta_df.reset_index(inplace=True)
        beta_df.set_index('ijk', inplace=True)
        beta_df.sort_index(inplace=True)
        # for each row in beta_zero
        for i, row in beta_df.iterrows():
            # get the indices of the row
            indices = [xyz_to_012[x] for x in i]
            # set the value of the tensor at the indices to the value of the row
            # by kleinman symmetry
            beta_tensor[indices[0], indices[1], indices[2]] = row['Beta']
            beta_tensor[indices[0], indices[2], indices[1]] = row['Beta']
            beta_tensor[indices[2], indices[0], indices[1]] = row['Beta']

        return beta_tensor

    def beta_proj(self, beta, E):
        return np.tensordot(beta, np.outer(E, E), axes=([0, 1], [0, 1])).T

        # write a function which takes in a molecule and returns the geometry and symbols
        # from the MadnessResponse Class

    def get_molecule_geometry(self, molecule, database_path):
        mra_mol = MadnessResponse(molecule, 'hf', 'dipole', database_path)
        molecule_dict = mra_mol.ground_info['molecule']
        geometry = molecule_dict['geometry']
        symbols = molecule_dict['symbols']
        return geometry, symbols

        # we will use plotly for this

    def ball_and_stick(self, molecule):
        mra_mol = MadnessResponse(molecule, self.quad_data.xc, self.quad_data.op,
                                  self.quad_data.database_path)
        molecule_dict = mra_mol.ground_info['molecule']
        geometry = molecule_dict['geometry']
        symbols = molecule_dict['symbols']

        x = []
        y = []
        z = []
        for atom in geometry:
            x.append(atom[0] * 1.0)
            y.append(atom[1] * 1.0)
            z.append(atom[2] * 1.0)

        # create a dictionary of colors for each atom
        colors = {'O': 'red', 'H': 'white', 'C': 'black', 'N': 'blue', 'Cl': 'green',
                  'Na': 'purple',
                  'F': 'orange', 'S': 'yellow', 'P': 'pink', 'I': 'brown', 'Br': 'cyan',
                  'B': 'grey',
                  'Ne': 'magenta', 'He': 'grey', 'Ar': 'grey', 'Kr': 'grey', 'Xe': 'grey',
                  'Li': 'grey', 'Mg': 'grey', 'Al': 'grey', 'Si': 'grey', 'K': 'grey', 'Ca': 'grey',
                  'Ti': 'grey', 'Be': 'grey', 'Fe': 'grey', 'Cu': 'grey', 'Zn': 'grey',
                  'Ag': 'grey', }
        colors = [colors[s] for s in symbols]
        sizes = [50 for s in symbols]
        # Create scatter plot for atoms
        scatter = go.Scatter3d(x=x, y=y, z=z, mode='markers', text=symbols,
                               marker=dict(color=colors, size=sizes, opacity=1.0, ))

        return scatter

    def plot_basis_sphere_and_vector(self, mol, basis_sets, omega_1=0, omega_2=0, radius=1,
                                     num_points=1000,
                                     sizeref=0.5, shift=1.0, sizemode='scaled', scene_length=5.0,
                                     basis_error=False):
        scatter = self.ball_and_stick(mol)
        quiver = self.unit_sphere_representation_basis(mol, basis_sets,
                                                       omega_1, omega_2, radius=radius,
                                                       num_points=num_points,
                                                       sizeref=sizeref, shift=shift,
                                                       sizemode=sizemode, basis_error=basis_error)
        quiver2 = self.vector_representation_basis(mol, basis_sets,
                                                   omega_1, omega_2, sizeref=0.5,
                                                   shift=shift,
                                                   sizemode='scaled')

        # Create the figure and plot
        fig = go.Figure(data=[scatter, quiver, quiver2])
        fig.update_layout(
            scene=dict(aspectmode='cube', xaxis=dict(range=[-scene_length, scene_length]),
                       yaxis=dict(range=[-scene_length, scene_length]),
                       zaxis=dict(range=[-scene_length, scene_length])))
        fig.show()


database_path = Path('/mnt/data/madness_data/august_no_symmetry')
# get the unique frequencies
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)

basis_sets = ['aug-cc-pVDZ', 'aug-cc-pVTZ', 'aug-cc-pVQZ', 'd-aug-cc-pVDZ', 'd-aug-cc-pVTZ',
              'd-aug-cc-pVQZ'] + ['aug-cc-pCVDZ', 'aug-cc-pCVTZ', 'aug-cc-pCVQZ', 'd-aug-cc-pCVDZ',
                                  'd-aug-cc-pCVTZ', 'd-aug-cc-pCVQZ']

mols = ['CH3SH', ]
xc = 'hf'
op = 'dipole'
mol = 'CH3SH'

freq = [0, ]
overwrite = False
rdb = QuadraticDatabase(mols, basis_sets, xc, op, freq, database_path, overwrite=overwrite)

print(rdb.q_df.query('molecule==@mol & Afreq==0.0'))
print(rdb.basis_set_error_df.query('molecule==@mol & Afreq==0.0'))
