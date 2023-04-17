import argparse
from pathlib import Path
import json

from quantumresponsepro.utils import generate_dalton_excited_state_json

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Example script to generate excited_state.json.")
    parser.add_argument("arg1", help="directory containing the database")
    parser.add_argument("-o", "--optional", help="Optional argument", default="default_value")

    args = parser.parse_args()

    database_path = Path(args.arg1)
    print('database path', database_path)

    # generates and saves molecule_excited_state_keys.json

    excited_state_json = generate_dalton_excited_state_json('hf', ['aug-cc-pVDZ', 'aug-cc-pVTZ'], database_path, True)
    print(excited_state_json)

    json_path = database_path.joinpath('json_data')
    if (not json_path.exists()):
        json_path.mkdir()

    dalton_excited_json_path = json_path.joinpath('dalton_excited_state.json')
    if dalton_excited_json_path.exists():
        old_excited_state_json = json.load(open(dalton_excited_json_path, 'r'))
        old_excited_state_json.update(excited_state_json)
        with open(dalton_excited_json_path, 'w') as f:
            json.dump(old_excited_state_json, f, indent=4)

    else:
        with open(dalton_excited_json_path, 'w') as f:
            json.dump(excited_state_json, f, indent=4)
