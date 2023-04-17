import argparse
import glob
from pathlib import Path

from utils import generate_excited_state_json

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Example script to generate excited_state.json.")
    parser.add_argument("arg1", help="directory containing the database")
    parser.add_argument("-o", "--optional", help="Optional argument", default="default_value")

    args = parser.parse_args()

    database_path = Path(args.arg1)

    molecule_path = database_path.joinpath("molecules")

    print('database path', database_path)
    print('molecule path', molecule_path)
    print('molecule path', molecule_path)

    mols = glob.glob('*.mol')

    # mol_files = glob.glob(f"{args.arg1}/*.mol")
