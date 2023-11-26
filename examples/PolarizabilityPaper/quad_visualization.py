from pathlib import Path
from DataAnalysisClass import *


october_absolute_kain_path = Path('/mnt/data/madness_data/october_beta_absolute_kain/')
paper_path = Path('/home/adrianhurtado/projects/writing/mra-tdhf-polarizability/Figures_v2')
paper_path = Path('/home/adrianhurtado/projects/writing/thesis2023/PresentationMaterials')

single = ['aug-cc-pVDZ', 'aug-cc-pVTZ', 'aug-cc-pVQZ', 'aug-cc-pV5Z', 'aug-cc-pV6Z']
single_polarized = ['aug-cc-pCVDZ', 'aug-cc-pCVTZ', 'aug-cc-pCVQZ']
double = ['d-aug-cc-pVDZ', 'd-aug-cc-pVTZ', 'd-aug-cc-pVQZ', 'd-aug-cc-pV5Z', 'd-aug-cc-pV6Z']
double_polarized = ['d-aug-cc-pCVDZ', 'd-aug-cc-pCVTZ', 'd-aug-cc-pCVQZ']
all_basis_sets = single + single_polarized + double + double_polarized

beta_path = Path('/mnt/data/madness_data/august_no_symmetry')
beta_mol = '/mnt/data/madness_data/august_no_symmetry/molecules/*.mol'
beta_mos = ['FNO',
            'CH3SH',
            'N2H2',
            'CH3NH2',
            'HOCl',
            'SiO',
            'HBS',
            'HBO',
            'CH2BH',
            'PH3O',
            'ClF',
            'HF',
            'BH2Cl',
            'Li2',
            'SiH4',
            'SF2',
            'P2H4',
            'NaCN',
            'OCl2',
            'CH3Cl',
            'SCl2',
            'HCCF',
            'LiH',
            'BF',
            'SiH3F',
            'HCONH2',
            'CS',
            'S2H2',
            'FCN',
            'NaCl',
            'SiH3Cl',
            'N2H4',
            'NH2Cl',
            'LiCl',
            'NH3O',
            'NH2OH',
            'BH3',
            'SO2',
            'HCN',
            'CH3BH2',
            'NaLi',
            'HOF',
            'CH3F',
            'HNS',
            'CH3OH',
            'Na2',
            'O3',
            'CH2NH',
            'ClCN',
            'Mg2',
            'NH3',
            'H2O',
            'PH3',
            'BHF2',
            'LiCN',
            'NH2F',
            'CH4',
            'HCl',
            'HCCCl',
            'HCOOH',
            'HCHS',
            'LiH_s',
            'HCHO',
            'SH2',
            'CSO',
            'HOOH',
            'CO',
            'HCP',
            'NOCl',
            'HNO',
            'NaH',
            'OF2',
            'HNC',
            'BH2F',
            'LiBH4']

beta_path = Path('/mnt/data/madness_data/august_no_symmetry')
quad_data = QuadraticDatabase(beta_mos, all_basis_sets, 'hf', 'dipole', [0], beta_path, False)
quad_data.save_dfs()
# save data

aug=['aug-cc-pVDZ', 'aug-cc-pVTZ', 'aug-cc-pVQZ', ]
daug=['d-aug-cc-pVDZ', 'd-aug-cc-pVTZ', 'd-aug-cc-pVQZ', ]

from QuadVisualization import QuadVisualization
qv=QuadVisualization(quad_data)

qv.HyperpolarizabilityScene('H2O',aug+daug,xshift=1.0,yshift=1.0,zshift=1.0,scene_length=10)


