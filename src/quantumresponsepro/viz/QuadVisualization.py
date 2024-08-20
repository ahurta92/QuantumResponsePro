import math

import numpy as np
import plotly.graph_objects as go
import vtk

from quantumresponsepro import MadnessResponse


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
            print('av', av)
            av.set_index('component', inplace=True)
            av.drop_duplicates(inplace=True)
            print('av', av)

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
