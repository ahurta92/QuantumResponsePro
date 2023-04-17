import argparse
from pathlib import Path

from utils import generate_dalton_excited_state_json

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Example script to generate excited_state.json.")
    parser.add_argument("arg1", help="directory containing the database")
    parser.add_argument("-o", "--optional", help="Optional argument", default="default_value")

    args = parser.parse_args()

    database_path = Path(args.arg1)
    print('database path', database_path)

    excited_state_json = generate_dalton_excited_state_json('hf', ['aug-cc-pVDZ'], database_path, False)
