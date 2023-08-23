from pathlib import Path
from quantumresponsepro import DatabaseGenerator

database_path = Path('/mnt/data/madness_data/fd_compare2')
basis = 'aug-cc-pVTZ'
xc = 'hf'
op = 'dipole'

dg = DatabaseGenerator(database_path)
freq_json = dg.get_frequency_json(9, xc, op, basis, .5)
