import os

from quantumresponsepro.BasisMRADataAssembler import *
from quantumresponsepro.BasisMRADataCollection import get_quad_df, get_polar_df


def round_to_n_significant_figures(x, n=3):
    if x == 0:
        return 0
    else:
        return round(x, -int(np.floor(np.log10(abs(x)))) + (n - 1))


def beta_df_np(beta_df):
    """
    Return the numpy array of the beta tensor
    """

    beta_df = beta_df.copy()
    xyz_to_012 = {'X': 0, 'Y': 1, 'Z': 2}
    beta_tensor = np.zeros((3, 3, 3))
    beta_df.reset_index(inplace=True)
    beta_df.set_index('ijk', inplace=True)
    beta_df.sort_index(inplace=True)
    # for each row in beta_zero
    for i, row in beta_df.iterrows():
        print(row)
        # get the indices of the row
        indices = [xyz_to_012[x] for x in i]
        # set the value of the tensor at the indices to the value of the row
        # by kleinman symmetry
        beta_tensor[indices[0], indices[1], indices[2]] = row['Beta']
        beta_tensor[indices[0], indices[2], indices[1]] = row['Beta']
        beta_tensor[indices[2], indices[0], indices[1]] = row['Beta']

    return beta_tensor


def get_beta_tensor(basis, mol_data, Bfreq, Cfreq):
    """
    Return the beta tensor for a given molecule, basis set, Bfreq and Cfreq
    """
    b_data = mol_data.query('basis == @basis & Bfreq == @Bfreq & Cfreq== @Cfreq').copy()
    b_data.drop_duplicates(inplace=True, subset=['ijk'])
    b_data.set_index('ijk', inplace=True)
    beta = beta_df_np(b_data)

    return beta


class QuadraticDatabase:
    def __init__(self, mols, basis_sets, xc, op, freq, database, overwrite=False):
        self.q_df = None
        self.basis_set_error_df = None
        self.vector_q_df = None
        self.vector_basis_set_error_df = None
        self.beta_hrs_df = None
        self.bhrs_basis_set_error_df = None
        self.pq_df = None
        self.freq = freq
        self.molecules = mols
        self.basis_sets = basis_sets
        self.xc = xc
        self.op = op
        self.database_path = database
        self.overwrite = overwrite

        self.initialize_dfs()

    def initialize_dfs(self):
        data_frame_attributes = ['q_df', 'vector_q_df', 'pq_df', ]

        for attr_name in data_frame_attributes:
            try:
                filename = self.database_path.joinpath(f"{attr_name}.csv")

                if not self.overwrite and os.path.exists(filename):
                    df = pd.read_csv(filename)
                else:
                    df = self.generate_df(attr_name)

                setattr(self, attr_name, df)
            except AttributeError as a:
                print(a)
                print(f"Could not initialize {attr_name}")
                continue
            except ValueError as v:
                print(v)
                print(f"Could not initialize {attr_name}")
                continue

    def generate_df(self, attr_name):
        # Replace this section with your specific data generating code
        if attr_name == 'q_df':
            return self.__generate_quad_df()
        elif attr_name == 'vector_q_df':
            return self.__generate_vector_q_df()
        elif attr_name == 'pq_df':
            return self.__generate_processed_quad_df()

    def save_dfs(self):
        data_frame_attributes = ['q_df', 'vector_q_df', 'pq_df'
                                 ]

        for attr_name in data_frame_attributes:
            df = getattr(self, attr_name).copy()
            filename = self.database_path.joinpath(f"{attr_name}.csv")
            df.to_csv(filename, index=False)

    def __generate_quad_df(self):
        q_df = get_quad_df(self.molecules, self.xc, self.op, self.database_path, self.basis_sets)
        print('q_df: ', q_df)
        print('q_df index: ', q_df.keys())
        q_df_i = q_df.copy()
        # truncate the Afreq Bfreq and Cfreq to 3 decimal places

        q_df_i['Afreq'] = q_df_i['Afreq'].apply(lambda x: round(x, 4))
        q_df_i['Bfreq'] = q_df_i['Bfreq'].apply(lambda x: round(x, 4))
        q_df_i['Cfreq'] = q_df_i['Cfreq'].apply(lambda x: round(x, 4))

        # Function to round to n significant figures

        # Apply the function to the float column
        q_df_i['Beta'] = q_df_i['Beta'].apply(
            lambda x: round_to_n_significant_figures(x, 6))

        # post process the molecule data
        f = []
        for mol in self.molecules:
            mol_data = q_df_i.query('molecule==@mol').copy()
            freq_map = mol_data.query('basis=="MRA"').Afreq.unique()
            mol_data['a'] = mol_data.Afreq.map(lambda x: mol_data.Afreq.unique().tolist().index(x))
            mol_data['b'] = mol_data.Bfreq.map(lambda x: mol_data.Afreq.unique().tolist().index(x))
            mol_data['c'] = mol_data.Cfreq.map(lambda x: mol_data.Afreq.unique().tolist().index(x))
            f.append(mol_data)
        q_df_i = pd.concat(f)

        return q_df_i.copy()

    def __generate_processed_quad_df(self):

        sample = self.q_df
        index = ['a', 'b', 'c', 'basis', 'molecule']
        sample = sample.set_index(index)
        # if beta is less than 1e-4 then set to zero
        sample['Beta'] = sample['Beta'].where(sample['Beta'].abs() > 5e-3, 0)
        xyz_to_012 = {'X': 0, 'Y': 1, 'Z': 2}
        beta_tensor = np.zeros((3, 3, 3))
        new_sample = pd.DataFrame()
        # for each individual index in sample
        for i in sample.index.unique():

            for j, row in sample.loc[i].iterrows():
                ijk = row.ijk
                beta = row.Beta
                # get the index of the tensor
                index_1 = [xyz_to_012[ijk[0]], xyz_to_012[ijk[1]], xyz_to_012[ijk[2]]]
                index_2 = [xyz_to_012[ijk[0]], xyz_to_012[ijk[2]], xyz_to_012[ijk[1]]]
                index_3 = [xyz_to_012[ijk[2]], xyz_to_012[ijk[0]], xyz_to_012[ijk[1]]]

                beta_tensor[tuple(index_1)] = beta
                beta_tensor[tuple(index_2)] = beta
                beta_tensor[tuple(index_3)] = beta
            # now create the reverse map
            d = {}
            for a in range(3):
                for b in range(3):
                    for c in range(3):
                        ijk = '{}{}{}'.format('XYZ'[a], 'XYZ'[b], 'XYZ'[c])
                        d[ijk] = beta_tensor[a, b, c]
            # now make dataframe with index i and d
            # where i is mulitindex of a,b,c,basis,molecule
            df = pd.Series(d)
            # set the index name to ijk
            df.index.name = 'ijk'
            # set the name of the series to Beta
            df.name = 'Beta'
            # set the index to be a,b,c,basis,molecule
            df = pd.DataFrame(df)
            # set multiindex to be a,b,c,basis,molecule
            df['a'] = i[0]
            df['b'] = i[1]
            df['c'] = i[2]
            df['basis'] = i[3]
            df['molecule'] = i[4]
            new_sample = pd.concat([new_sample, df])
        new_sample.reset_index(inplace=True)
        return new_sample

    def beta_tensor(self, mol, basis_set, Bfreq, Cfreq):
        mol_data, afreq = self.__get_mol_and_freq_data(mol)
        df_i = get_beta_tensor(basis_set, mol_data, Bfreq, Cfreq)
        return df_i

    def __get_mol_and_freq_data(self, mol):
        mol_data = self.q_df.query('molecule == @mol').copy()
        freq = list(mol_data['Afreq'].unique())
        return mol_data, freq

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

            beta_tensor = beta_df_np(bdata)
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
            try:
                filename = self.database.joinpath(f"{attr_name}.csv")
                if not self.overwrite and os.path.exists(filename):
                    df = pd.read_csv(filename)
                else:
                    df = self.generate_df(attr_name)

                setattr(self, attr_name, df)
            except AttributeError as a:
                print(a)
                print(f"Could not initialize {attr_name}")
                pass
            except ValueError as v:
                print(v)
                print(f"Could not initialize {attr_name}")
                pass

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
            try:
                df = getattr(self, attr_name).copy()
                filename = self.database.joinpath(f"{attr_name}.csv")
                df.to_csv(filename, index=False)
            except AttributeError as a:
                print(a)
                print(f"Could not save {attr_name}")
                pass

    def make_all_dfs_detailed(self):
        data_frame_attributes = ['energy_df', 'energy_diff_df', 'polar_data', 'alpha_eigen',
                                 'eigen_diff', 'iso_data']

        for attr_name in data_frame_attributes:
            df = getattr(self, attr_name)  # Access the DataFrame attribute
            print(df)
            modified_df = make_detailed_df(df)  # Apply the function
            setattr(self, attr_name, modified_df)  #
