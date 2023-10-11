import math

import matplotlib as mpl
import os
import plotly.graph_objects as go
import vtk
from quantumresponsepro.BasisMRADataAssembler import *
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
        self.beta_hrs_df = None
        self.bhrs_basis_set_error_df = None
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
                                 'vector_basis_set_error_df', 'beta_hrs_df',
                                 'bhrs_basis_set_error_df']

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
        elif attr_name == 'beta_hrs_df':
            return self.__generate_bhrs_df()
        elif attr_name == 'bhrs_basis_set_error_df':
            return self.__generate_bhrs_basis_error_df()

    def save_dfs(self):
        data_frame_attributes = ['q_df', 'basis_set_error_df', 'vector_q_df',
                                 'vector_basis_set_error_df', 'beta_hrs_df',
                                 'bhrs_basis_set_error_df']

        for attr_name in data_frame_attributes:
            df = getattr(self, attr_name).copy()
            filename = self.database_path.joinpath(f"{attr_name}.csv")
            df.to_csv(filename, index=False)

    def make_all_dfs_detailed(self):
        data_frame_attributes = ['q_df', 'basis_set_error_df', 'vector_q_df',
                                 'vector_basis_set_error_df', 'beta_hrs_df',
                                 'bhrs_basis_set_error_df']

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
        vector_df.dropna(inplace=True)
        vector_df.reset_index(inplace=True)
        return vector_df.copy()

    def bhrs_df(self, mol, basis_set, b, c):
        mol_data, afreq = self.__get_mol_and_freq_data(mol)
        Bfreq = afreq[b]
        Cfreq = afreq[c]
        beta_hrs = self.__get_beta_hrs(basis_set, mol_data, Bfreq, Cfreq)

        return beta_hrs

    def beta_tensor(self, mol, basis_set, Bfreq, Cfreq):
        mol_data, afreq = self.__get_mol_and_freq_data(mol)
        df_i = self.__get_beta_tensor(basis_set, mol_data, Bfreq, Cfreq)
        return df_i

    def __generate_bhrs_df(self):
        # Generate Hyperrayleigh Scattering from beta data

        vector_df = []
        for mol in self.molecules:
            mol_data, afreq = self.__get_mol_and_freq_data(mol)
            for basis in self.basis_sets + ['MRA']:
                for b in range(len(self.freq)):
                    bf = afreq[b]
                    for c in range(b, len(self.freq)):
                        cf = afreq[c]
                        try:
                            df_i = self.__get_beta_hrs(basis, mol_data, bf, cf)
                            vector_df.append(df_i)
                        except IndexError as e:
                            print(e)
                            print(mol, basis)
                            continue
        bshr_df = pd.concat(vector_df, axis=1).T
        return bshr_df

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

    def __generate_bhrs_basis_error_df(self):

        basis_error_df = pd.DataFrame()
        for mol in self.molecules:
            mol_data, afreq = self.__get_bhrs_mol_and_freq_data(mol)
            for basis in self.basis_sets:
                for b in range(len(self.freq)):
                    bf = afreq[b]
                    for c in range(b, len(self.freq)):
                        cf = afreq[c]
                        try:
                            df_i = self.__get_quad_bshrs_basis_error_data_df(basis, mol_data, bf,
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

    def __get_bhrs_mol_and_freq_data(self, mol):
        mol_data = self.beta_hrs_df.query('molecule == @mol').copy()
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

    def __get_quad_bshrs_basis_error_data_df(self, basis, mol_data, Bfreq, Cfreq):
        """
        This function returns
        """
        b_data = mol_data.query('basis == @basis & Bfreq == @Bfreq & Cfreq== @Cfreq').copy()

        # remove the nan values
        b_data.dropna(inplace=True)
        # remove Beta less than 1e-4
        m_data = mol_data.query('basis == "MRA" & Bfreq==@Bfreq & Cfreq == @Cfreq ').copy()
        # remove the Beta less than .001
        m_data.dropna(inplace=True)
        # remove the Beta less than 1e-4
        BC_index = ['Bfreq', 'Cfreq']

        m_data.set_index(BC_index, inplace=True)
        m_data.drop_duplicates(inplace=True, )
        b_data.drop_duplicates(inplace=True)

        b_data = b_data.query('Beta.abs() > .005').copy()
        m_data = m_data.query('Beta.abs() > .005').copy()
        # b_data.drop_duplicates(inplace=True, subset=['ijk'])

        basis_beta = b_data.Beta.values[0]
        mra_beta = m_data.Beta.values[0]
        print('basis_beta: ', basis_beta)
        print('mra_beta: ', mra_beta)
        rel_error = (basis_beta - mra_beta) / mra_beta * 100
        print('b_data: ', b_data)
        print('rel error: ', rel_error)
        b_data['Beta'] = rel_error
        # insert rel error into b_data series

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
        print(b_data)
        # the vector representation of the hyperpolarizability tensor accumulates components
        # whose output are in the same direction.
        # for example ax=beta(xxx)+beta(xyy)+beta(xzz) and ay=beta(yxx)+beta(yyy)+beta(yzz) and
        # az=beta(zxx)+beta(zyy)+beta(zzz)
        # the vector representation is then [ax,ay,az]
        print(b_data)

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

    def __get_beta_hrs(self, basis, mol_data, Bfreq, Cfreq):
        """
        This function returns
        """
        b_data = mol_data.query('basis == @basis & Bfreq == @Bfreq & Cfreq== @Cfreq').copy()
        molecule = b_data['molecule'].unique()[0]
        b_data.drop_duplicates(inplace=True, subset=['ijk'])
        b_data.set_index('ijk', inplace=True)
        beta = self.beta_df_np(b_data)
        dipolar, octupolar, beta_hrs = self.__get__beta_dipolar_and_octupolar(beta)
        # set index
        index = ['molecule', 'basis', 'Afreq', 'Bfreq', 'Cfreq', 'dipolar', 'octupolar', 'Beta']
        Afreq = Bfreq + Cfreq
        beta_hrs = pd.Series([molecule, basis, Afreq, Bfreq, Cfreq, dipolar, octupolar,
                              beta_hrs], index=index)

        return beta_hrs

    def __get__beta_dipolar_and_octupolar(self, beta):

        term1 = 0
        term2 = 0
        term3 = 0
        term4 = 0
        term5 = 0

        for i in range(3):
            term1 += beta[i, i, i] ** 2
            for j in range(3):
                if i != j:
                    term2 += beta[i, i, i] * beta[i, j, j]
                    term3 += beta[j, i, i] ** 2
                for k in range(3):
                    if i != j and i != k and j != k:
                        term4 += beta[j, i, i,] * beta[j, k, k]
                        term5 += beta[i, j, k] ** 2

        dipolar = 3 / 5 * term1 + 6 / 5 * term2 + 3 / 5 * term3 + 3 / 5 * term4
        octupolar = 2 / 5 * term1 - 6 / 5 * term2 + 12 / 5 * term3 - 3 / 5 * term4 + term5
        beta_hrs = np.sqrt(10 / 45 * dipolar + 10 / 105 * octupolar)
        return (dipolar, octupolar, beta_hrs)

    def __get__beta_octupolar(self, beta):

        term1 = 0
        term2 = 0
        term3 = 0
        term4 = 0

        for i in range(3):
            term1 += beta(i, i, i) ** 2
            for j in range(3):
                if i != j:
                    term2 += beta(i, i, i) * beta(i, j, j)
                    term3 += beta(j, i, i) ** 2
                for k in range(3):
                    if i != j and i != k and j != k:
                        term4 += beta(j, i, i, ) * beta(j, k, k)
        term1 = 3 / 5 * term1
        term2 = 6 / 5 * term2
        term3 = 3 / 5 * term2
        term4 = 3 / 5 * term2

    def __get_beta_tensor(self, basis, mol_data, Bfreq, Cfreq):
        """
        This function returns
        """
        b_data = mol_data.query('basis == @basis & Bfreq == @Bfreq & Cfreq== @Cfreq').copy()
        molecule = b_data['molecule'].unique()[0]
        b_data.drop_duplicates(inplace=True, subset=['ijk'])
        b_data.set_index('ijk', inplace=True)
        print(b_data)
        # the vector representation of the hyperpolarizability tensor accumulates components
        # whose output are in the same direction.
        # for example ax=beta(xxx)+beta(xyy)+beta(xzz) and ay=beta(yxx)+beta(yyy)+beta(yzz) and
        # az=beta(zxx)+beta(zyy)+beta(zzz)
        # the vector representation is then [ax,ay,az]
        print(b_data)
        beta = self.beta_df_np(b_data)
        print(beta)

        return beta

    def beta_df_np(self, beta_df):

        beta_df = beta_df.copy()
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

    def get_beta_points(self, mol, basis_sets, omega_1, omega_2, radius=1,
                        num_points=1000,
                        basis_error=False, xshift=0.0,
                        yshift=0.0, zshift=0.0):

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
            print(bdata)
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
            p2 = results
            x[i * num_points:(i + 1) * num_points] = p1[:, 0]
            y[i * num_points:(i + 1) * num_points] = p1[:, 1]
            z[i * num_points:(i + 1) * num_points] = p1[:, 2]

            u[i * num_points:(i + 1) * num_points] = p2[:, 0]
            v[i * num_points:(i + 1) * num_points] = p2[:, 1]
            w[i * num_points:(i + 1) * num_points] = p2[:, 2]

            coords = np.zeros((len(basis_sets) * num_points, 3))
            vector_vals = np.zeros((len(basis_sets) * num_points, 3))

            coords[:, 0] = x
            coords[:, 1] = y
            coords[:, 2] = z
            vector_vals[:, 0] = u
            vector_vals[:, 1] = v
            vector_vals[:, 2] = w
            # make an array of the points and vectors
            return coords, vector_vals

    def beta_to_vtk(self, mol, basis_sets, omega_1=0.0, omega_2=0.0, radius=1,
                    num_points=1000,
                    basis_error=False, xshift=0.0,
                    yshift=0.0, zshift=0.0):

        for basis in basis_sets:
            coords, vector_vals = self.get_beta_points(mol, [basis], omega_1, omega_2, radius,
                                                       num_points,
                                                       basis_error,
                                                       0,
                                                       0, 0)

            # Initialize VTK points and vector data structures
            points = vtk.vtkPoints()
            vectors = vtk.vtkFloatArray()
            vectors.SetNumberOfComponents(3)
            vectors.SetName("Hyperpolarizability")

            # Assuming `coords` is a list of (x, y, z) coordinates and `vector_vals` is a list of (vx, vy, vz) vectors
            for (x, y, z), (vx, vy, vz) in zip(coords, vector_vals):
                print(vx, vy, vz)
                points.InsertNextPoint(x, y, z)
                vectors.InsertNextTuple((vx, vy, vz))

            # Set up the polydata object
            polydata = vtk.vtkPolyData()
            polydata.SetPoints(points)
            polydata.GetPointData().SetVectors(vectors)

            # Write to VTK file
            writer = vtk.vtkXMLPolyDataWriter()

            file_name = self.q_df.database_path.joinpath(f'{mol}_{basis}_{omega_1}_'
                                                         f'{omega_2}.vtp')
            writer.SetFileName(file_name)
            writer.SetInputData(polydata)
            writer.Write()


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
            return get_energy_diff_data(self.molecules, self.xc, self.op, self.basis_sets,
                                        self.database)
        elif attr_name == 'polar_data':
            return get_polar_df(self.molecules, self.xc, self.op, self.database, self.basis_sets)
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


def plot_BHRS_diff(data, valence, ax, color_pal=None):
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


class basis_set_analysis:

    def __init__(self, qda: QuadraticDatabase, pda: PolarizabilityData):
        self.qda = qda
        self.qda.make_all_dfs_detailed()
        self.pda = pda
        self.pda.make_all_dfs_detailed()

    def basis_set_analysis_plot(self, molecule, valence, omega_1=0.0, omega_2=0.0, type=None,
                                q_measure='BHRS'):

        r_plot, axes = plt.subplots(1, 3, figsize=(9, 4), frameon=True, layout='constrained')

        if type is None:
            pal = light_pal
        else:
            pal = simple_pal

        e_diff = self.pda.energy_diff_df.query('molecule == @molecule').copy()
        a_diff = self.pda.eigen_diff.query('molecule == @molecule & omega==@omega_1').copy()
        # rename alpha to diff
        a_diff.rename(columns={'alpha': 'diff'}, inplace=True)

        plot_energy_diff(e_diff, valence, axes[0], color_pal=pal)
        alpha_plot = plot_component_diff(a_diff, valence, 'ij', axes[1],
                                         color_pal=pal)

        axes[0].set_title('Energy')

        axes[1].set_title(r'$\alpha_{ij}(0;0,0)$')
        axes[1].set_yscale('symlog', linthresh=1e-2, linscale=0.25)

        if q_measure == 'BHRS':

            q_diff = self.qda.bhrs_basis_set_error_df.query('molecule == @molecule & '
                                                            'Bfreq==@omega_1 & '
                                                            'Cfreq==@omega_2').copy()
            q_diff.rename(columns={'Beta': 'diff'}, inplace=True)
            print(q_diff)
            beta_plot = plot_BHRS_diff(q_diff, valence, axes[2],
                                       color_pal=pal)
            axes[2].set_title(
                r'$\beta_{HRS}$' + '({};{},{})'.format(omega_1 + omega_2, omega_1, omega_2))
        elif q_measure == 'vector':
            q_diff = self.qda.vector_basis_set_error_df.query('molecule == @molecule & '
                                                              'Bfreq==@omega_1 & '
                                                              'Cfreq==@omega_2').copy()

            q_diff.rename(columns={'Beta': 'diff'}, inplace=True)
            beta_plot = plot_component_diff(a_diff, valence, 'component', axes[2],
                                            color_pal=pal)
            print(q_diff)
            axes[2].set_title(
                r'$V$' + '({};{},{})'.format(omega_1 + omega_2, omega_1, omega_2))
        elif q_measure == 'ijk':
            q_diff = self.qda.basis_set_error_df.query('molecule == @molecule & '
                                                       'Bfreq==@omega_1 & '
                                                       'Cfreq==@omega_2').copy()
            q_diff.rename(columns={'Beta': 'diff'}, inplace=True)
            print(q_diff)
            beta_plot = plot_component_diff(a_diff, valence, 'ijk', axes[1],
                                            color_pal=pal)
            axes[2].set_title(
                r'$\beta_{ijk}$' + '({};{},{})'.format(omega_1 + omega_2, omega_1, omega_2))

        # axes[2].set_yscale('symlog', linthresh=1e-2, linscale=0.25)
        # axes[2].set_yscale('symlog', linthresh=1e-2, linscale=0.25)
        r_plot.suptitle(molecule)

        return r_plot, axes


class QuadVisualization:

    def __init__(self, data: QuadraticDatabase):
        self.quad_data = data

    def unit_sphere_representation(self, data, basis, omega_1, omega_2, radius=1, num_points=1000,
                                   basis_error=False, colormap='Jet', sizeref=0.5):

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
                                         colormap='Jet', sizeref=0.5,
                                         sizemode='scaled', basis_error=False, xshift=0.0,
                                         yshift=0.0, zshift=0.0):

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
            print(bdata)
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
            x[i * num_points:(i + 1) * num_points] = p1[:, 0] + xshift * (i % 4) - (2 * xshift)
            y[i * num_points:(i + 1) * num_points] = p1[:, 1] + yshift * (math.floor(i / 4)) - (
                    2 * yshift)
            z[i * num_points:(i + 1) * num_points] = p1[:, 2] + zshift

            u[i * num_points:(i + 1) * num_points] = p2[:, 0]
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
                                    sizeref=0.5, sizemode='scaled', xshift=0.0,
                                    yshift=0.0, zshift=0.0, colorbar_title='Beta'):

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

            x[i] = 0 + xshift * (i % 4)
            y[i] = 0 + yshift * (i * math.floor(i / 4))
            z[i] = 0 + zshift
            u[i] = av.loc['x'].Beta
            v[i] = av.loc['y'].Beta
            w[i] = av.loc['z'].Beta

        quiver2 = go.Cone(x=x, y=y, z=z, u=u, v=v, w=w, colorscale='Jet', sizemode=sizemode,
                          sizeref=sizeref, text=text, opacity=1.0,
                          colorbar=dict(title=colorbar_title))
        return quiver2

    def unit_sphere_representation_basis_error(self, mol, basis_sets, omega_1, omega_2,
                                               radius=1,
                                               num_points=1000,
                                               colormap='Jet', sizeref=0.5,
                                               sizemode='scaled', xshift=0.0, yshift=0.0,
                                               zshift=0.0):

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
            x[i * num_points:(i + 1) * num_points] = p1[:, 0] + xshift * (i % 4)
            y[i * num_points:(i + 1) * num_points] = p1[:, 1] + yshift * (math.floor(i / 4))
            z[i * num_points:(i + 1) * num_points] = p1[:, 2] + zshift

            u[i * num_points:(i + 1) * num_points] = p2[:, 0]
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

    def ball_and_stick(self, molecule, xshift=1.0, yshift=0.0, zshift=0.0):
        mra_mol = MadnessResponse(molecule, self.quad_data.xc, self.quad_data.op,
                                  self.quad_data.database_path)
        molecule_dict = mra_mol.ground_info['molecule']
        geometry = molecule_dict['geometry']
        symbols = molecule_dict['symbols']

        x = []
        y = []
        z = []
        for atom in geometry:
            x.append(atom[0] * 1.0 + xshift)
            y.append(atom[1] * 1.0 + yshift)
            z.append(atom[2] * 1.0 + zshift)

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
                                     sizeref=0.5, sizemode='scaled', scene_length=5.0,
                                     basis_error=False, xshift=0.0, yshift=0.0, zshift=0.0):

        scatter = self.ball_and_stick(mol, xshift, zshift=zshift)
        quiver = self.unit_sphere_representation_basis(mol, basis_sets,
                                                       omega_1, omega_2, radius=radius,
                                                       num_points=num_points,
                                                       sizeref=sizeref, xshift=xshift,
                                                       sizemode=sizemode, basis_error=basis_error)
        quiver2 = self.vector_representation_basis(mol, basis_sets,
                                                   omega_1, omega_2, sizeref=0.5,
                                                   sizemode='scaled', xshift=xshift, yshift=yshift,
                                                   zshift=0)

        # Create the figure and plot
        fig = go.Figure(data=[scatter, quiver, quiver2])
        fig.update_layout(
            scene=dict(aspectmode='cube', xaxis=dict(range=[-scene_length, scene_length]),
                       yaxis=dict(range=[-scene_length, scene_length]),
                       zaxis=dict(range=[-scene_length, scene_length])))
        fig.show()

    def HyperpolarizabilityScene(self, mol, basis_sets, omega_1=0, omega_2=0, radius=1,
                                 num_points=1000,
                                 sizeref=0.5, sizemode='scaled', scene_length=5.0,
                                 xshift=0.0, yshift=0.0, zshift=0.0,
                                 ):

        molecule_viz = self.ball_and_stick(mol, xshift=0.0, yshift=0.0, zshift=3 * zshift)
        MRA_quiver = self.unit_sphere_representation_basis(mol, ['MRA'],
                                                           omega_1, omega_2, radius=radius,
                                                           num_points=num_points,
                                                           sizeref=2 * sizeref, xshift=0.0,
                                                           yshift=0.0, zshift=5.0,
                                                           sizemode=sizemode, basis_error=False)
        MRA_vector = self.vector_representation_basis(mol, basis_sets,
                                                      omega_1, omega_2, sizeref=3.5,
                                                      sizemode='absolute', xshift=0,
                                                      yshift=0, zshift=1.5 * zshift)

        basis_Error = self.unit_sphere_representation_basis(mol, basis_sets,
                                                            omega_1, omega_2, radius=radius,
                                                            num_points=num_points,
                                                            sizeref=5 * sizeref,
                                                            sizemode=sizemode, basis_error=True,
                                                            xshift=xshift, yshift=yshift,
                                                            zshift=0.0)

        # Create the figure and plot
        fig = go.Figure(data=[molecule_viz, MRA_quiver, MRA_vector, basis_Error])
        fig.update_layout(
            scene=dict(aspectmode='cube', xaxis=dict(range=[-scene_length, scene_length]),
                       yaxis=dict(range=[-scene_length, scene_length]),
                       zaxis=dict(range=[-scene_length, scene_length])))
        fig.show()
        return fig

    def get_beta_points(self, mol, basis_sets, omega_1=0.0, omega_2=0.0, radius=1,
                        num_points=1000,
                        basis_error=False, xshift=0.0,
                        yshift=0.0, zshift=0.0):

        # get the frequencies from the database

        if not basis_error:
            print('omega_1', omega_1)
            print('omega_2', omega_2)

            data = self.quad_data.q_df.query(
                'molecule==@mol & Bfreq==@omega_1 & Cfreq==@omega_2 ')

        else:
            data = self.quad_data.basis_set_error_df.query(
                'molecule==@mol & Cfreq==@omega_1 & Bfreq==@omega_2 ')

        print(data)
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

            bdata = data.query('basis==@basis')
            print(bdata)
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
            p2 = results
            x[i * num_points:(i + 1) * num_points] = p1[:, 0] + xshift * (i % 4) - (2 * xshift)
            y[i * num_points:(i + 1) * num_points] = p1[:, 1] + yshift * (math.floor(i / 4)) - (
                    2 * yshift)
            z[i * num_points:(i + 1) * num_points] = p1[:, 2] + zshift

            u[i * num_points:(i + 1) * num_points] = p2[:, 0]
            v[i * num_points:(i + 1) * num_points] = p2[:, 1]
            w[i * num_points:(i + 1) * num_points] = p2[:, 2]

            coords = np.zeros((len(basis_sets) * num_points, 3))
            vector_vals = np.zeros((len(basis_sets) * num_points, 3))

            coords[:, 0] = x
            coords[:, 1] = y
            coords[:, 2] = z
            vector_vals[:, 0] = u
            vector_vals[:, 1] = v
            vector_vals[:, 2] = w
            # make an array of the points and vectors
            return coords, vector_vals

    def beta_to_vtk(self, mol, basis_sets, omega_1=0, omega_2=0, radius=1,
                    num_points=1000,
                    basis_error=False, xshift=0.0,
                    yshift=0.0, zshift=0.0):

        for basis in basis_sets:

            coords, vector_vals = self.get_beta_points(mol, [basis], omega_1=omega_1,
                                                       omega_2=omega_2, radius=radius,
                                                       num_points=num_points,
                                                       basis_error=basis_error,
                                                       xshift=xshift,
                                                       yshift=yshift, zshift=zshift)

            # Initialize VTK points and vector data structures
            points = vtk.vtkPoints()
            vectors = vtk.vtkFloatArray()
            vectors.SetNumberOfComponents(3)
            vectors.SetName("vectors")

            # Assuming `coords` is a list of (x, y, z) coordinates and `vector_vals` is a list of (vx, vy, vz) vectors
            for (x, y, z), (vx, vy, vz) in zip(coords, vector_vals):
                points.InsertNextPoint(x, y, z)
                vectors.InsertNextTuple((vx, vy, vz))

            # Set up the polydata object
            polydata = vtk.vtkPolyData()
            polydata.SetPoints(points)
            polydata.GetPointData().SetVectors(vectors)

            # Write to VTK file
            writer = vtk.vtkXMLPolyDataWriter()

            # set up vtk folder if it doesn't exist
            vtk_folder = self.quad_data.database_path.joinpath('vtk')
            if not vtk_folder.exists():
                vtk_folder.mkdir()
            # set up vtk/mol folder name
            mol_folder = self.quad_data.database_path.joinpath('vtk', mol)
            if not mol_folder.exists():
                mol_folder.mkdir()
            # i need to enforce omega_1 and omega_2 to have 3 decimal points

            file_name = mol_folder.joinpath(f'{mol}_{basis}_{omega_1:.3f}_{omega_2:.3f}.vtp')
            writer.SetFileName(file_name)
            writer.SetInputData(polydata)
            writer.Write()

    def beta_to_vtk_basis_error(self, mol, basis_sets, omega_1=0, omega_2=0, radius=1,
                                num_points=1000,
                                xshift=0.0,
                                yshift=0.0, zshift=0.0):

        for basis in basis_sets:
            coords, vector_vals = self.get_beta_points(mol, [basis], omega_1=omega_1,
                                                       omega_2=omega_2, radius=radius,
                                                       num_points=num_points,
                                                       basis_error=True,
                                                       xshift=xshift,
                                                       yshift=yshift, zshift=zshift)

            # Initialize VTK points and vector data structures
            points = vtk.vtkPoints()
            vectors = vtk.vtkFloatArray()
            vectors.SetNumberOfComponents(3)
            vectors.SetName("vectors")

            # Assuming `coords` is a list of (x, y, z) coordinates and `vector_vals` is a list of (vx, vy, vz) vectors
            for (x, y, z), (vx, vy, vz) in zip(coords, vector_vals):
                points.InsertNextPoint(x, y, z)
                vectors.InsertNextTuple((vx, vy, vz))

            # Set up the polydata object
            polydata = vtk.vtkPolyData()
            polydata.SetPoints(points)
            polydata.GetPointData().SetVectors(vectors)

            # Write to VTK file
            writer = vtk.vtkXMLPolyDataWriter()

            # set up vtk folder if it doesn't exist
            vtk_folder = self.quad_data.database_path.joinpath('vtk')
            if not vtk_folder.exists():
                vtk_folder.mkdir()
            # set up vtk/mol folder name
            mol_folder = self.quad_data.database_path.joinpath('vtk', mol)
            if not mol_folder.exists():
                mol_folder.mkdir()
            # i need to enforce omega_1 and omega_2 to have 3 decimal points

            file_name = mol_folder.joinpath(f'{mol}_{basis}_{omega_1:.3f}_{omega_2:.3f}.vtp')
            writer.SetFileName(file_name)
            writer.SetInputData(polydata)
            writer.Write()
