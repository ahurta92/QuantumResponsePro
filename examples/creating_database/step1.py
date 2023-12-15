from pathlib import Path
from quantumresponsepro import DatabaseGenerator

# Step 1 is to generate the excited state data for the molecules in the database

database_path = Path("/mnt/data/madness_data/database_director")

db_gen = DatabaseGenerator(database_path)
xc = "hf"
basis = 'aug-cc-pVTZ'
db_gen.get_dalton_excited_json(xc, [basis], True)
# looks for molecules in the database directory and generates the excited-state data from which to
# generate frequency.json input file.


# freq_json = dg.get_frequency_json(9, xc, op, basis, .5)
