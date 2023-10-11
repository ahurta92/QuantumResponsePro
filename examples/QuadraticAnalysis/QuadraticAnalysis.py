import math

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from pathlib import Path
from quantumresponsepro import MadnessResponse
from quantumresponsepro.BasisMRADataCollection import get_quad_df, get_polar_df


def get_basis_error_data_df(mol, basis, database_path, mol_data, Bfreq, Cfreq):
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


def get_basis_df(basis, mol_data, Bfreq, Cfreq):
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


def get_MRA_data_df(mol_data, Bfreq, Cfreq):
    # remove Beta less than 1e-4
    m_data = mol_data.query('basis == "MRA" & Bfreq==@Bfreq & Cfreq == @Cfreq ').copy()
    # remove the Beta less than .001
    m_data.dropna(inplace=True)
    m_data = m_data.query('Beta.abs() > .005').copy()

    return m_data


def round_to_n_significant_figures(x, n=3):
    if x == 0:
        return 0
    else:
        return round(x, -int(np.floor(np.log10(abs(x)))) + (n - 1))


class QuadraticAnalysisData:
    database_path = None
    mol = None
    basis_sets = None
    xc = None
    op = None
    q_df = None
    beta_df = None
    basis_set_error_df = None

    def __init__(self, mols, basis_sets, xc, op, database):
        self.molecules = mols
        self.basis_sets = basis_sets
        self.xc = xc
        self.op = op
        self.database_path = database

        alpha_df = get_polar_df(mols, xc, op, database_path, basis_sets)
        print(alpha_df)

        q_df = get_quad_df(mols, xc, op, database_path, basis_sets)
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

        self.q_df = q_df_i.copy()

    def get_mol_and_freq_data(self, mol):
        mol_data = self.q_df.query('molecule == @mol').copy()
        freq = list(mol_data['Afreq'].unique())
        return mol_data, freq

    def get_molecules_basis_error_df(self, mol, basis, Bfreq_i, Cfreq_j):
        mol_data, freq = self.get_mol_and_freq_data(mol)
        Bfreq = freq[Bfreq_i]
        Cfreq = freq[Cfreq_j]
        b_data = get_basis_error_data_df(mol, basis, self.database_path, mol_data, Bfreq,
                                         Cfreq)
        return b_data

    def get_mra_df(self, mol, Bfreq_i, Cfreq_j):
        mol_data, freq = self.get_mol_and_freq_data(mol)
        Bfreq = freq[Bfreq_i]
        Cfreq = freq[Cfreq_j]
        mra_data = get_MRA_data_df(mol_data, Bfreq, Cfreq)
        return mra_data

    def get_basis_df(self, mol, basis, Bfreq_i, Cfreq_j):
        mol_data, freq = self.get_mol_and_freq_data(mol)
        Bfreq = freq[Bfreq_i]
        Cfreq = freq[Cfreq_j]
        mra_data = get_basis_df(basis, mol_data, Bfreq, Cfreq)
        return mra_data

    def get_dataset_df(self, mol, basis_sets, Bfreq_i, Cfreq_j):
        df = pd.DataFrame()
        for basis in basis_sets:
            b_data = self.get_molecules_basis_error_df(mol, basis, Bfreq_i, Cfreq_j)
            df = pd.concat([df, b_data])

        return df

    def get_molecule_frequency_data(self, mol, basis_sets, Bfreqs, Cfreqs):
        '''
        returns beta dataframe at a given frequency for a set of basis sets
        '''
        df = pd.DataFrame()

        for basis in basis_sets:
            for b in Bfreqs:
                for c in Cfreqs:
                    df = pd.concat([df, self.get_basis_df(mol, basis, b, c)])

        return df

    def get_molecule_frequency_basis_error_data(self, mol, basis_sets, Bfreqs, Cfreqs):
        '''
        returns the basis set error data frame for a set of frequencies and basis sets
        '''
        df = pd.DataFrame()

        mol_data, freq = self.get_mol_and_freq_data(mol)

        for basis in basis_sets:
            for b in Bfreqs:
                bfreq = freq[b]
                for c in Cfreqs:
                    cfreq = freq[c]
                    b_data = get_basis_error_data_df(mol, basis, self.database_path, mol_data,
                                                     bfreq, cfreq)
                    df = pd.concat([df, b_data])

        return df


database_path = Path('/mnt/data/madness_data/august_no_symmetry')

# get the unique frequencies
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)


class QuadVisualization:

    def __init__(self, data: QuadraticAnalysisData):
        self.quad_data = data

    def unit_sphere_representation(self, data, basis, omega_1, omega_2, radius=1, num_points=1000,
                                   basis_error=False, colormap='Blues', sizeref=0.5):

        data = data.query('basis==@basis & Afreq==@omega_1 & Bfreq==@omega_2')
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
                                         sizemode='scaled'):

        data = self.quad_data.get_molecule_frequency_data(mol, basis_sets, [omega_1], [omega_2])

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

    def vector_representation_basis(self, basis_sets, omega_1, omega_2, colormap='Greens',
                                    sizeref=0.5, shift=1.0, sizemode='scaled'):

        data = self.quad_data.get_molecule_frequency_data(mol, basis_sets, [omega_1], [omega_2])
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
            bdata = data.query('basis==@basis & Afreq==@omega_1 & Bfreq==@omega_2')
            av = self.vector_representation(bdata)
            x[i] = 0 + shift * (i - len(basis_sets) / 2)
            y[i] = 0
            z[i] = 0
            u[i] = av[0] + shift * (i - len(basis_sets) / 2)
            v[i] = av[1]
            w[i] = av[2]

        quiver2 = go.Cone(x=x, y=y, z=z, u=u, v=v, w=w, colorscale='Blues', sizemode=sizemode,
                          sizeref=sizeref, text=text, opacity=1.0)
        return quiver2

    def unit_sphere_representation_basis_error(self, mol, basis_sets, omega_1, omega_2,
                                               radius=1,
                                               num_points=1000,
                                               colormap='Blues', sizeref=0.5, shift=1.0,
                                               sizemode='scaled'):

        basis_error_df = self.quad_data.get_molecule_frequency_basis_error_data(mol,
                                                                                basis_sets,
                                                                                [omega_1],
                                                                                [omega_2])

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
            print(bdata)
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

    def vector_representation(self, beta_df):
        beta_vector = np.zeros(3, float)
        beta_df = beta_df.copy()
        beta_df.reset_index(inplace=True)
        beta_df.set_index('ijk', inplace=True)

        ax_indices = ['XXX', 'XYY', 'XZZ']
        ay_indices = ['YXX', 'YYY', 'YZZ']
        az_indices = ['ZXX', 'ZYY', 'ZZZ']

        for i in range(3):
            print(i)
            try:
                f = beta_df.loc[ax_indices[i]].Beta
                print(f)
                beta_vector[0] += f
            except KeyError as e:
                print(e)
                print('failed')
        for i in range(3):
            try:
                f = beta_df.loc[ay_indices[i]].Beta
                print(f)
                beta_vector[1] += f
            except KeyError as e:
                print(e)
        for i in range(3):
            try:
                beta_vector[2] += beta_df.loc[az_indices[i]].Beta
            except KeyError as e:
                print(e)
        print(beta_vector)
        return beta_vector

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


basis_sets = ['aug-cc-pVDZ', 'aug-cc-pVTZ', 'aug-cc-pVQZ', 'd-aug-cc-pVDZ', 'd-aug-cc-pVTZ',
              'd-aug-cc-pVQZ']

mols = ['CH3SH', 'H2O', 'CH3OH', 'C2H2', 'C2H4', 'CH3F', 'CH3OH', 'NaCl']
mol = 'H2O'

qa = QuadraticAnalysisData(mols, basis_sets, 'hf', 'dipole', database_path)
viz = QuadVisualization(qa)
scatter = viz.ball_and_stick(mol)
quiver = viz.unit_sphere_representation_basis(mol, ['MRA', ],
                                              0, 0, num_points=1000,
                                              sizeref=5.0, shift=0.0, sizemode='scaled')
quiver2 = viz.vector_representation_basis(['MRA'], 0, 0, sizeref=.5, shift=5.0,
                                          sizemode='absolute')

# Create the figure and plot
fig = go.Figure(data=[scatter, quiver, quiver2])
fig.update_layout(scene=dict(aspectmode='cube', xaxis=dict(range=[-5, 5]),
                             yaxis=dict(range=[-5, 5]),
                             zaxis=dict(range=[-5, 5])))
fig.show()

quiver = viz.unit_sphere_representation_basis_error(mol,
                                                    [
                                                        'aug-cc-pVQZ', 'd-aug-cc-pVQZ'],
                                                    num_points=250, omega_1=0, omega_2=0,
                                                    sizeref=.5, shift=2.5, sizemode='absolute')

# Create the figure and plot
fig = go.Figure(data=[scatter, quiver, ])
fig.update_layout(scene=dict(aspectmode='cube', xaxis=dict(range=[-5, 5]),
                             yaxis=dict(range=[-5, 5]),
                             zaxis=dict(range=[-5, 5])))
fig.show()
basis_error_df = pd.DataFrame()
for mol in mols:
    for basis in basis_sets:
        try:
            df_i = qa.get_molecules_basis_error_df(mol, basis, 0, 0)
            basis_error_df = pd.concat([basis_error_df, df_i])
        except IndexError as e:
            print(e)
            print(mol, basis)
            continue
basis_error_df.dropna(inplace=True)
print(basis_error_df)

polar_keys = ['xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy', 'zz']

alpha_basis_error_df = pd.DataFrame()
for mol in mols:
    for basis in basis_sets:
        try:
            df_i = qa.get_molecules_basis_error_df(mol, basis, 0, 0)
            basis_error_df = pd.concat([basis_error_df, df_i])
        except IndexError as e:
            print(e)
            print(mol, basis)
            continue
basis_error_df.dropna(inplace=True)
print(basis_error_df)
