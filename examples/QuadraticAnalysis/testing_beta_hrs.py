from pathlib import Path
from quantumresponsepro.BasisMRADataAssembler import *

from QuadraticData_class import PolarizabilityData
from QuadraticData_class import QuadVisualization
from QuadraticData_class import QuadraticDatabase
from QuadraticData_class import basis_set_analysis

database_path = Path('/mnt/data/madness_data/august_no_symmetry')
# get the unique frequencies
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)

basis_sets = ['aug-cc-pVDZ', 'aug-cc-pCVDZ', 'd-aug-cc-pVDZ', 'd-aug-cc-pCVDZ'] \
             + ['aug-cc-pVTZ', 'aug-cc-pCVTZ', 'd-aug-cc-pVTZ', 'd-aug-cc-pCVTZ'] + \
             ['aug-cc-pVQZ', 'aug-cc-pCVQZ', 'd-aug-cc-pVQZ', 'd-aug-cc-pCVQZ']
mols = ['H2O', 'CH3OH', 'C2H2', 'C2H4', 'CH3F', 'CH3OH', 'NaLi', 'LiCN']
xc = 'hf'
op = 'dipole'

freq = [0, 1, 2, 3, 4, 5, 6, 7, 8]
overwrite = False
rdb = QuadraticDatabase(mols, basis_sets, xc, op, freq, database_path, overwrite=overwrite)
mol = 'LiCN'
beta_hrs = rdb.bhrs_df(mol, 'MRA', 0, 0)
basis = 'aug-cc-pVDZ'
rdb.save_dfs()

polarizability_database = Path('/mnt/data/madness_data/post_watoc/august')
polar_data = PolarizabilityData(mols, 'hf', 'dipole', basis_sets=basis_sets,
                                database=polarizability_database,
                                overwrite=overwrite)
b_plotter = basis_set_analysis(rdb, polar_data)
for mol in mols:
    r_plot, axes = b_plotter.basis_set_analysis_plot(mol, ['D', 'T', 'Q', ], type=None)
    r_plot.show()

mol = 'LiCN'
print(rdb.beta_hrs_df.query('molecule ==@mol & basis == "MRA"'))
print(rdb.q_df.query('molecule ==@mol & basis == "MRA"'))

viz = QuadVisualization(rdb)

om1 = 0
freq = list(rdb.q_df.query('molecule ==@mol & basis == "MRA"').Bfreq.unique())
print(rdb.q_df.query('molecule ==@mol & basis == "MRA"'))
print(freq)
for om1 in range(0, 4):
    om_1 = freq[om1]
    viz.beta_to_vtk(mol,
                    ['MRA', ],
                    omega_1=om_1, omega_2=om_1, radius=5.0, basis_error=False)

basis_sets = ['aug-cc-pVDZ', 'aug-cc-pCVDZ', 'd-aug-cc-pVDZ', 'd-aug-cc-pCVDZ'] \
             + ['aug-cc-pVTZ', 'aug-cc-pCVTZ', 'd-aug-cc-pVTZ', 'd-aug-cc-pCVTZ'] + \
             ['aug-cc-pVQZ', 'aug-cc-pCVQZ', 'd-aug-cc-pVQZ', 'd-aug-cc-pCVQZ']
for om1 in range(0, 4):
    om_1 = freq[om1]
    viz.beta_to_vtk(mol,
                    basis_sets,
                    omega_1=om_1, omega_2=om_1, radius=1.0, basis_error=True)
