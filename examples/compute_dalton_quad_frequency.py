import os
import sys
from quantumresponsepro import DaltonRunner
from pathlib import Path

mol = sys.argv[1]
num_proc = sys.argv[2]

# basis_list = ['aug-cc-pV5Z-uc','aug-cc-pV6Z']  # , 'aug-cc-pV5Z', 'aug-cc-pV6Z']
# d_basis_list = ['d-aug-cc-pV5Z-uc','d-aug-cc-pV6Z']  # , 'd-aug-cc-pV5Z', 'd-aug-cc-pV6Z']
# basis_list = ['q-aug-cc-pVDZ','q-aug-cc-pVTZ','q-aug-cc-pVQZ','q-aug-cc-pVD5',]  # , 'd-aug-cc-pV5Z', 'd-aug-cc-pV6Z']
aug_basis = ['aug-cc-pVDZ', 'aug-cc-pVTZ', 'aug-cc-pVQZ']
daug_basis = ['d-aug-cc-pVDZ', 'd-aug-cc-pVTZ', 'd-aug-cc-pVQZ']
aug_pc_basis = ['aug-cc-pCVDZ', 'aug-cc-pCVTZ', 'aug-cc-pCVQZ', ]
daug_pc_basis = ['d-aug-cc-pCVDZ', 'd-aug-cc-pCVTZ', 'd-aug-cc-pCVQZ', ]

basis_list = aug_basis + daug_basis + aug_pc_basis + daug_pc_basis
#basis_list = ['aug-cc-pVDZ', 'aug-cc-pVTZ']
BASEDIR = Path(str(os.getcwd()))
runner = DaltonRunner(BASEDIR, True)
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
