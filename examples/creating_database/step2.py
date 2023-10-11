from pathlib import Path
from quantumresponsepro import DatabaseGenerator

database_path = Path('/mnt/data/madness_data/nacl_ioniazation')
basis = 'aug-cc-pVTZ'
xc = 'hf'
op = 'dipole'

dg = DatabaseGenerator(database_path)
freq_json = dg.get_frequency_json(9, xc, op, basis, .5)
