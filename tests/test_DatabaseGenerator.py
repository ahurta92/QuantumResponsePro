from quantumresponsepro import DatabaseGenerator
from pathlib import Path
import os
import json


def test_generate_frequeny_json():
    dir_path = Path(os.path.dirname(os.path.realpath(__file__)))
    database_path = dir_path.joinpath('../example_database')
    dg = DatabaseGenerator(database_path)
    freq_json = dg.get_frequency_json(8, 'hf', 'dipole', 'aug-cc-pVTZ')

    freq_test_path = database_path.joinpath('molecules/frequency.json')
    with open(freq_test_path, 'r') as f:
        freq_test = json.load(f)

    assert freq_json['Be']['hf']['dipole'] == freq_test['Be']['hf']['dipole']
