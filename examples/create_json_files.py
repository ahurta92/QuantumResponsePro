import argparse
from pathlib import Path

from quantumresponsepro.utils import DatabaseGenerator

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Example script to generate excited_state.json.")
    parser.add_argument("arg1", help="directory containing the database")
    parser.add_argument("-o", "--optional", help="Optional argument", default="default_value")

    args = parser.parse_args()

    database_path = Path(args.arg1)
    print('database path', database_path)

    dg = DatabaseGenerator(database_path)
    excited_json = dg.get_dalton_excited_json('hf', ['aug-cc-pVDZ'], run=True)
    try:
        freq_json = dg.get_frequency_json(8, 'hf', 'dipole', 'aug-cc-pVTZ', .5)
    except KeyError as k:
        excited_json = dg.get_dalton_excited_json('hf', ['aug-cc-pVTZ'], run=True)
        freq_json = dg.get_frequency_json(8, 'hf', 'dipole', 'aug-cc-pVTZ', .5)
    excited_json = dg.get_dalton_excited_json('hf', ['aug-cc-pVDZ'], run=True)
