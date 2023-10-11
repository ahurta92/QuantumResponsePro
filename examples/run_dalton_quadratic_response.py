from pathlib import Path
from quantumresponsepro import DaltonRunner

mol = "H2O"

xc = "hf"
op = 'dipole'
basis = 'aug-cc-pVDZ'

database_path = Path("/mnt/data/madness_data/post_watoc/august")

# create a DaltonRunner object and set run_dalton to True
dr = DaltonRunner(database_path, True)

dr = dr.get_quad_json(mol, xc, op, basis)
print(dr['Quad'])
print(dr)
