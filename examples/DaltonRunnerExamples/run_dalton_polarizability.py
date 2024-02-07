from pathlib import Path

from quantumresponsepro import DaltonRunner

# Runs a Linear Response Dalton calculation for a given molecule and set of basis sets
# give a good name to this script
# define molecule, exchange-correlation functional, and operator


mol = "H2"
xc = "hf"
op = 'dipole'

database_path = Path("/mnt/data/madness_data/post_watoc/alrich_test_set")
# create a DaltonRunner object and set run_dalton to True
dr = DaltonRunner(database_path, True)

all_basis = []
zetas = ['D', 'T', 'Q', '5', '6']

for zeta in zetas:
    basis = ["cc-pV{}Z".format(zeta), "aug-cc-pV{}Z".format(zeta), 'aug-cc-pCV{}Z'.format(zeta),
             'd-aug-cc-pV{}Z'.format(zeta),
             'd-aug-cc-pCV{}Z'.format(zeta), 't-aug-cc-pV{}Z'.format(zeta), ]
    all_basis = all_basis + basis
triple = ['t-aug-cc-pVDZ', 't-aug-cc-pVTZ', 't-aug-cc-pVQZ', 't-aug-cc-pCVDZ',
          't-aug-cc-pCVTZ',
          't-aug-cc-pCVQZ', 't-aug-cc-pV5Z', 't-aug-cc-pV6Z']
results = {}
for b in ['def2-QZVPPD']:
    try:
        results[b] = dr.get_polar_json(mol, xc, op, b)
    except (FileNotFoundError, KeyError):
        print(f"File not found for {b}")
        continue
