from quantumresponsepro import DaltonRunner
from pathlib import Path

mol = "Ne"
basis = ["aug-cc-pCVDZ", "d-aug-cc-pCVDZ"]
xc = "hf"
op = 'dipole'

database = Path("/mnt/data/madness_data/post_watoc/august")

# create a DaltonRunner object and set run_dalton to True
dr = DaltonRunner(database, True)

# run Dalton for the given molecule, basis, and xc
results = {}
for b in basis:
    results[b] = dr.get_polar_json(mol, xc, op, b)
