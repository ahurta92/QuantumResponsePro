import pytest
from quantumresponsepro import MadnessResponse


def test_reading_madness_response():
    # Replace with appropriate test input and expected output for your package
    mol = 'F2'
    op = 'dipole'
    xc = 'hf'

    mol_response = MadnessResponse(mol, xc, op, "../example_database")

    assert mol_response.converged[0]
    "Madness input translation failed"
