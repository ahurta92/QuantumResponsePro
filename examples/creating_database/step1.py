from pathlib import Path
from quantumresponsepro import DatabaseGenerator

# Step 1 is to generate the excited state data for the molecules in the database

dir = "/mnt/data/madness_data/nacl_ioniazation"

db_gen = DatabaseGenerator(Path(dir))

xc = "hf"
basis = 'aug-cc-pVTZ'

db_gen.get_dalton_excited_json(xc, [basis], True)

# freq_json = dg.get_frequency_json(9, xc, op, basis, .5)
