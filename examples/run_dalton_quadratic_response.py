import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

from quantumresponsepro import DaltonRunner, BasisMRADataCollection, BasisMRADataAnalyzer
from quantumresponsepro.BasisMRADataAssembler import make_detailed_df
from pathlib import Path

mol = "H2O"

xc = "hf"
op = 'dipole'
basis = 'aug-cc-pVDZ'

database_path = Path("/gpfs/projects/rjh/adrian/post_watoc/august")

# create a DaltonRunner object and set run_dalton to True
dr = DaltonRunner(database_path, True)

dr = dr.get_quad_json(mol, xc, op, basis)
print(dr['Quad'])
print(dr)
