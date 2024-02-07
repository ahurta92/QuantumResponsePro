import sys
import os
from pathlib import Path
from quantumresponsepro import DaltonRunner

# Runs a Quadratic Response Dalton calculation for a given molecule and set of basis sets
# Results are stored in Dalton directory from the current working directory

mol = sys.argv[1]
num_proc = sys.argv[2]

# Define basis sets
taug_basis = ['t-aug-cc-pVDZ', 't-aug-cc-pVTZ', 't-aug-cc-pVQZ']
q_aug_basis = ['q-aug-cc-pVDZ', 'q-aug-cc-pVTZ', 'q-aug-cc-pVQZ']
basis_list = taug_basis + q_aug_basis

BASEDIR = Path(str(os.getcwd()))
# Create DaltonRunner object
runner = DaltonRunner(BASEDIR, True)
# Set the number of processors
runner.Np = num_proc

for basis in basis_list:
    try:
        result = runner.get_quad_json(mol, 'hf', 'dipole', basis)
        print(result)
        print(BASEDIR)
        print(os.curdir)

    except FileNotFoundError as f:
        print(f)
        os.chdir(BASEDIR)
        pass
