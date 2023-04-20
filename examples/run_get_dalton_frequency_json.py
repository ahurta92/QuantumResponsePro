from quantumresponsepro import DatabaseGenerator
from pathlib import Path
import os

dir_path = Path(os.path.dirname(os.path.realpath(__file__)))
database_path = dir_path.joinpath('../example_database')

dg = DatabaseGenerator(database_path)

single_VZ = ['aug-cc-pVDZ', 'aug-cc-pVTZ', 'aug-cc-pVQZ']
double_VZ = ['d-aug-cc-pVDZ', 'd-aug-cc-pVTZ', 'd-aug-cc-pVQZ']

single_CZ = ['aug-cc-pCVDZ', 'aug-cc-pCVTZ', 'aug-cc-pCVQZ']
double_CZ = ['d-aug-cc-pCVDZ', 'd-aug-cc-pCVTZ', 'd-aug-cc-pCVQZ']

basis_list = single_VZ + double_VZ + single_CZ + double_CZ

dalton_dipole_json = dg.get_dalton_frequency_json('hf', 'dipole', basis_list, run=True)
