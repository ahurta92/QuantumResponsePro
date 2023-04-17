# QuantumResponsePro

QuantumResponsePro is a Python package for computing and analyzing dynamic response properties of molecules in TDDFT and
TDHF formalisms with arbitrary accuracy. It enables comparison of results computed using the MADNESS and Dalton quantum
chemistry packages by translating generic inputs into specific MADNESS and Dalton input files, and processing Dalton
output files into easily readable JSON files.

## Features

- Input and output translators for MADNESS and Dalton quantum chemistry packages.
- Frequency-dependent polarizability analysis for various molecules at multiple frequencies.
- Extensive basis set analysis based on augmentation, valence size, core polarization, etc.
- Compatibility with other response properties beyond polarizability.
- Includes example scripts and Jupyter notebooks for demonstration purposes.

## Installation

Install QuantumResponsePro via pip:

```bash
pip install quantumresponsepro
```

Alternatively, you can clone the repository and install it from the source:

```bash
git clone https://github.com/yourusername/QuantumResponsePro.git
cd QuantumResponsePro
pip install -e .
```

## Usage

### Generate a response database

1. Create a response database directory with the following subdirectories:
   molecules: contains the molecule definition madness.mol in xyz format
   json_data: stores database property json files (generated later)
   To run response property calculations, generate a frequency.json file with a list of frequencies for computing
   response properties. First, calculate the first excited-state energy for a given functional and save the results in
   dalton_excited.json.

Use this code:

```python
from pathlib import Path
from quantumresponsepro import DatabaseGenerator

dir = "path/to/directory"
db_gen = DatabaseGenerator(Path(dir))
xc = "hf"
basis = 'aug-cc-pVDZ'

db_gen.DatabaseGenerator.get_dalton_excited_json(xc, [basis], True)
```

The code generates input files for the first excited-state calculation, runs the calculation, and saves the results for
each molecule in dalton/xc/excited_state/
Output files have the prefix `excited_mol_basis`, with .out for output files and .json for json files.
All molecule-basis pair json files are combined into a single `dalton_excited.json` file in the **json_data** directory.

Once the excited-state calculations are complete and a `dalton_excited.json` file has been generated you
can now generate a `frequency.json` file with a list of frequencies for computing response properties.
The `frequency.json` file is used by the MADNESS `mad_freq` executable and dalton to set the frequencies for
the response property calculations.

### Generate frequency.json file

```python
from pathlib import Path
from quantumresponsepro import DatabaseGenerator

database_path = Path('path/to/directory')
dg = DatabaseGenerator(database_path)
excited_json = dg.get_dalton_excited_json('hf', ['aug-cc-pVDZ'], run=True)
try:
    freq_json = dg.get_frequency_json('hf', 'dipole', 'aug-cc-pVTZ')
except KeyError as k:
    excited_json = dg.get_dalton_excited_json('hf', ['aug-cc-pVTZ'], run=True)
    freq_json = dg.get_frequency_json('hf', 'dipole', 'aug-cc-pVTZ')
excited_json = dg.get_dalton_excited_json('hf', ['aug-cc-pVDZ'], run=True)
```

In this example, we generate a `frequency.json` file for the **HF** functional, **dipole** response, and **aug-cc-pVTZ**
basis set.
The `frequency.json` file defines the frequencies for response property calculations for a specific functional and
operator.
The frequencies are generated from the excited-state energies of the basis set.

The `get_frequency_json` function first checks if the `dalton_excited.json` file exists and whether data for a specific
functional and basis set is present. If found, it overwrites the frequency.json file. Otherwise, it raises an exception.

Before running the `get_frequency_json` function, make sure to run the `get_dalton_excited_json` function for the
desired
functional and basis set, as shown in the example.

### Run MADNESS and Dalton calculations with `frequency.json` file

3. Create a json file for the database properties. The json file should contain the following properties:

After generating the MADNESS and Dalton files, you can use the following code to generate the output response dataframe:

```python
from quantumresponsepro import (
    ResponseDataBundle,
)

response_dataframes = ResponseDataBundle(dir)

```

For more detailed examples, please refer to the examples directory and the provided Jupyter notebooks.

Contributing
We welcome contributions to QuantumResponsePro! Please fork the repository, create a new branch with your changes, and
submit a pull request. Make sure to follow the established code style and add appropriate tests and documentation for
your changes.

## License

QuantumResponsePro is released under the [MIT License](LICENSE.txt).
.



