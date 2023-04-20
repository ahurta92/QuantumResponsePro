from quantumresponsepro import DatabaseGenerator
from pathlib import Path
import os
import json


def test_generate_frequeny_json():
    dir_path = Path(os.path.dirname(os.path.realpath(__file__)))
    database_path = dir_path.joinpath('../example_database')
    dg = DatabaseGenerator(database_path)
    freq_json = dg.get_frequency_json(8, 'hf', 'dipole', 'aug-cc-pVTZ', .5)

    freq_test_path = database_path.joinpath('molecules/frequency.json')
    with open(freq_test_path, 'r') as f:
        freq_test = json.load(f)

    assert freq_json['Be']['hf']['dipole'] == freq_test['Be']['hf']['dipole']


def test_get_dalton_dipole_json():
    op = 'dipole'
    dir_path = Path(os.path.dirname(os.path.realpath(__file__)))
    database_path = dir_path.joinpath('../example_database')
    dg = DatabaseGenerator(database_path)
    dipole_json = dg.get_dalton_frequency_json('hf', op, ['aug-cc-pVDZ'], run=True)

    freq_test_path = database_path.joinpath('molecules/dalton-dipole.json')
    with open(freq_test_path, 'r') as f:
        dipole_test = json.load(f)

    new = dipole_json['Be']['hf']['dipole']['aug-cc-pVDZ']['response']['values']['xx']
    old = dipole_test['Be']['hf']['dipole']['aug-cc-pVDZ']['response']['xx']

    assert new == old
