from quantumresponsepro import BasisMRADataCollection
from quantumresponsepro import BasisMRADataAnalyzer
import seaborn as sns

from pathlib import Path

august = Path('/mnt/data/madness_data/post_watoc/august')
paper_path = Path('/home/adrianhurtado/projects/writing/mra-tdhf-polarizability/Figures_v2')
database = BasisMRADataCollection(august)
analyzer = BasisMRADataAnalyzer(database, .02)

vls = ['D', 'T', 'Q']

mol = "HOOH"
omega = [0, 8]
g = analyzer.plot_iso_valence_convergence(mol, 'alpha', ['D', 'T', 'Q'], omega)
g.fig.show()

g1 = analyzer.plot_alpha_component_convergence(mol, valence=['D', 'T', 'Q', '5'])
g1.fig.show()

mol = "NaCl"

g1 = analyzer.plot_alpha_component_convergence(mol, valence=['D', 'T', 'Q', '5'])
g1.fig.show()

mol = "F2"

g1 = analyzer.plot_alpha_component_convergence(mol, valence=['D', 'T', 'Q', '5'])
g1.fig.show()

g1 = analyzer.plot_alpha_component_convergence('NaH', valence=['D', 'T', 'Q', '5'])
g1.fig.show()

g1 = analyzer.plot_alpha_component_convergence('H2', valence=['D', 'T', 'Q', '5'])
g1.fig.show()

g1 = analyzer.plot_alpha_component_convergence('H2O', valence=['D', 'T', 'Q', '5'])
g1.fig.show()



g1 = analyzer.plot_alpha_component_convergence('NH3', valence=['D', 'T', 'Q', '5'])
g1.fig.show()

