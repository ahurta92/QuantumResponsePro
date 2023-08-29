import pandas as pd
import matplotlib as mpl
from quantumresponsepro import BasisMRAData, DaltonRunner

from quantumresponsepro import BasisMRADataAnalyzer
from quantumresponsepro import MadnessResponse

import seaborn as sns

from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# Create some example data
import numpy as np
from quantumresponsepro.BasisMRADataAssembler import make_detailed_df


def tensor_to_numpy(j):
    array = np.empty(j["size"])
    array[:] = j["vals"]
    return np.reshape(array, tuple(j["dims"]))


august_data = Path('/mnt/data/madness_data/post_watoc/august_high_prec')

database = BasisMRAData(august_data, new=False)

quad_data = database.all_quad_data.copy()

molecule = 'H2O'
quad_data['ijk'] = quad_data['A'] + quad_data['B'] + quad_data['C']
# sort by the component column
quad_data.sort_values(by='ijk', inplace=True)
# drop the ABC columns
quad_data.drop(columns=['A', 'B', 'C'], inplace=True)

h2o_data = quad_data.query('molecule==@molecule & basis=="MRA"').query(
    'Afreq==0 & Bfreq==0 & Cfreq==0')
print('mad_data', h2o_data)

h2o_basis_data = quad_data.query('molecule==@molecule & basis=="aug-cc-pVTZ"').query(
    'Afreq==0 & Bfreq==0 & Cfreq==0')
print('basis_data', h2o_basis_data.drop_duplicates())

h2o_mad = MadnessResponse(molecule, 'hf', 'dipole', august_data)

mad_dipole = tensor_to_numpy(h2o_mad.ground_info['scf_dipole_moment'])

runner = DaltonRunner(august_data, False)
b_ground, b_response = runner.get_frequency_result(molecule, 'hf', 'dipole', 'd-aug-cc-pVQZ')

basis_dipole = np.array(b_ground['dipole'], dtype=float)
print(
    'basis dipole:', basis_dipole
)

# normalize the dipole moment
basis_dipole /= np.linalg.norm(basis_dipole)
print('basis dipole:', basis_dipole)

mad_dipole /= np.linalg.norm(mad_dipole)
print('mad dipole:', mad_dipole)

# now compute the cross product between them
cross = np.cross(basis_dipole, mad_dipole)
# if the norm of the cross product is zero, then the vectors are parallel
# i need to see the vector to any vector that is perpendicular to either of them
if np.linalg.norm(cross) < 1e-5:
    print('parallel')
    random_vector = np.random.uniform(0, 1, 3)
    random_vector /= np.linalg.norm(random_vector)
    print('random vector:', random_vector)
    cross = np.cross(basis_dipole, random_vector)

# now compute the angle between them
angle = np.arccos(np.dot(basis_dipole, mad_dipole))

print('cross:', cross)
print('angle:', angle)

from scipy.spatial.transform import Rotation as R

r = R.from_rotvec(cross * angle)
rotation_matrix = r.as_matrix()
print(rotation_matrix)

# multiply the dipole moments together.  If the sum of the elements is positive then they are
# aligned, if negative then they are anti-aligned then they we multiple baisis set beta by -1
measure = basis_dipole * mad_dipole
print(measure)
print(np.sum(measure))
