from quantumresponsepro import MadnessResponse
from quantumresponsepro import FrequencyData
from pathlib import Path

mol = "He"
xc = "hf"
operator = "dipole"

data_dir_1 = Path("/mnt/data/madness_data/development")
data_dir_2 = Path("/mnt/data/madness_data/example_database")

mad_r2 = MadnessResponse(mol, xc, operator, data_dir_1)

p = mad_r2.compare_mra_convergence(data_dir_2, mad_r2.frequencies[8], 'alpha', ['xx', 'yy', 'zz'])
print(p)
