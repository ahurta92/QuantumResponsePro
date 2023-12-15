from pathlib import Path
from quantumresponsepro import DatabaseGenerator, DaltonRunner

database_path = Path('/mnt/data/madness_data/database')
basis = 'aug-cc-pVTZ'
xc = 'hf'
op = 'dipole'

# the basis set determines frequency of the excited state
dg = DatabaseGenerator(database_path)
freq_json = dg.get_frequency_json(9, xc, op, basis, .5)
# creating frequency.json file

dr = DaltonRunner(database_path, True)

single=['aug-cc-p{zeta}Z'.format(zeta=zeta) for zeta in ['D', 'T', 'Q', '5', '6']]
double=['d-aug-cc-p{zeta}Z'.format(zeta=zeta) for zeta in ['D', 'T', 'Q', '5', '6']]
triple=['t-aug-cc-p{zeta}Z'.format(zeta=zeta) for zeta in ['D', 'T', 'Q', '5', '6']]
basis = single + double + triple

mol=['H2', 'Ne', 'H2O', 'SO2']
results = []
for b in basis:
    try:
        results[b] = dr.get_polar_json(mol, xc, op, b)
    except (FileNotFoundError, KeyError):
        print(f"File not found for {b}")
        continue

print(results)

