from quantumresponsepro import BasisMRADataCollection
from quantumresponsepro import BasisMRADataAnalyzer
from quantumresponsepro import Tabler
import seaborn as sns

from pathlib import Path

august = Path('/mnt/data/madness_data/post_watoc/august')
paper_path = Path('/home/adrianhurtado/projects/writing/mra-tdhf-polarizability/Figures_v2')
database = BasisMRADataCollection(august)
tabler = Tabler(database)

col_format = 'l' + 'S[\
table-auto-round = true,\
retain-zero-exponent = false,\
table-number-alignment=center,\
table-format = 1.2e-1,]' * 12
e1 = tabler.get_styled_summary('energy')
e1.to_latex(
    paper_path.joinpath(Path("energy.tex")),
    hrules=True,
    convert_css=True,
    multicol_align='|c|',
    column_format=col_format,
    siunitx=True
)

a1 = tabler.get_styled_summary('alpha')
a1.to_latex(paper_path.joinpath(Path("alpha.tex")),
            convert_css=True,
            column_format=col_format,
            multicol_align='|c|',
            hrules=True,
            siunitx=True)
g1 = tabler.get_styled_summary('gamma')
g1.to_latex(paper_path.joinpath(Path("gamma.tex")),
            hrules=True,
            column_format=col_format,
            multicol_align='|c|',
            convert_css=True,
            siunitx=True)

alpha1, alpha2 = tabler.get_styled_molecules('alpha')
print(alpha1)
alpha1.to_latex(paper_path.joinpath(Path("alpha_molecules_1.tex")),
                hrules=True,
                multicol_align='|c|',
                convert_css=True,
                siunitx=True)

alpha2.to_latex(paper_path.joinpath(Path("alpha_molecules_2.tex")),
                hrules=True,
                multicol_align='|c|',
                convert_css=True,
                siunitx=True)
