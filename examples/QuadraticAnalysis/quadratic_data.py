from quantumresponsepro import BasisMRAData
from quantumresponsepro import BasisMRADataAnalyzer
from quantumresponsepro import MadnessResponse

import seaborn as sns

from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# Create some example data
import numpy as np

fd_path = Path('/mnt/data/madness_data/fd_compare3')
thesis_path = Path('/home/adrianhurtado/projects/writing/thesis2023/Figures_v2')
paper_path = Path('/home/adrianhurtado/projects/writing/mra-tdhf-polarizability/Figures_v2')
tromso_path = Path('/home/adrianhurtado/projects/writing/tromso_poster/figures')

paper_path = paper_path

database = BasisMRAData(fd_path.joinpath('low-low'), new=True)
print(database.all_polar_data)

print(database.molecules)
print(database.basis_sets)

print(database.available_molecules)

nacl = MadnessResponse('NaCl', 'hf', 'dipole', fd_path.joinpath('low-low'))
print(nacl.quad_data)

print(database.all_quad_data)


class DataViewGenerator:
    def __init__(self, database):
        self.database = database
