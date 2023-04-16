# QuantumResponsePro

QuantumResponsePro is a Python package for computing and analyzing dynamic response properties of molecules in TDDFT and TDHF formalisms with arbitrary accuracy. It enables comparison of results computed using the MADNESS and Dalton quantum chemistry packages by translating generic inputs into specific MADNESS and Dalton input files, and processing Dalton output files into easily readable JSON files.

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

```python
from quantumresponsepro import (
    Dalton,
    MadnessResponseReader,
    ResponseDataBundle,
)

response_dataframes=ResponseDataBundle(dir)

```
For more detailed examples, please refer to the examples directory and the provided Jupyter notebooks.

Contributing
We welcome contributions to QuantumResponsePro! Please fork the repository, create a new branch with your changes, and submit a pull request. Make sure to follow the established code style and add appropriate tests and documentation for your changes.

## License

QuantumResponsePro is released under the [MIT License](LICENSE.txt).
.



