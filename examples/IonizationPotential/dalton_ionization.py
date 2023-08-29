import json

import seaborn as sns
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from quantumresponsepro import DaltonRunner
from pathlib import Path
from quantumresponsepro.BasisMRADataAssembler import make_detailed_df

mol = 'Cl'
charge = 1
xc = 'hf'
basis = 'cc-pVDZ'

database = Path("/mnt/data/madness_data/nacl_ioniazation")

# Create DaltonRunner object

dalton = DaltonRunner(database, run_new=True)

data = dalton.get_energy_json(mol, xc, basis, charge)


def compute_basis_ionization_potential(mol, xc, basis):
    neutral = dalton.get_energy_json(mol, xc, basis, 0)
    cation = dalton.get_energy_json(mol, xc, basis, 1)

    print(neutral)

    neutral_energy = neutral[basis]['ground']['calculationResults']['totalEnergy']['value']
    cation_energy = cation[basis]['ground']['calculationResults']['totalEnergy']['value']
    print(neutral_energy)

    return cation_energy - neutral_energy


def compute_basis_electron_affinity(mol, xc, basis):
    neutral = dalton.get_energy_json(mol, xc, basis, 0)
    anion = dalton.get_energy_json(mol, xc, basis, -1)

    neutral_energy = neutral[basis]['ground']['calculationResults']['totalEnergy']['value']
    anion_energy = anion[basis]['ground']['calculationResults']['totalEnergy']['value']

    return anion_energy - neutral_energy


na_ionization = compute_basis_ionization_potential('Na', 'hf', 'cc-pVDZ')
cl_electron_affinity = compute_basis_electron_affinity('Cl', 'hf', 'cc-pVDZ')

print("Na ionization potential: ", na_ionization)
print("Cl electron affinity: ", cl_electron_affinity)


# Lets grab the madness results for Na

#

def get_madness_energy(database_path, mol_name, xc, ):
    json_path = database_path.joinpath('{}/{}/moldft.calc_info.json'.format(xc, mol_name))
    with open(json_path) as json_file:
        json_data = json.load(json_file)

    return json_data['return_energy']


na_energy = get_madness_energy(database, 'Na', 'hf')
na_cation_energy = get_madness_energy(database, 'Na_m1', 'hf')

na_ionization = na_cation_energy - na_energy
print("Na ionization potential: ", na_ionization)

cl_energy = get_madness_energy(database, 'Cl', 'hf')
cl_anion_energy = get_madness_energy(database, 'Cl_p1', 'hf')

cl_electron_affinity = cl_anion_energy - cl_energy
print("Cl electron affinity: ", cl_electron_affinity)

#
no_aug = ['cc-pVDZ', 'cc-pVTZ', 'cc-pVQZ', 'cc-pV5Z', 'cc-pV6Z']
no_aug_core_polarized = ['cc-pCVDZ', 'cc-pCVTZ', 'cc-pCVQZ', 'cc-pC5Z', 'cc-pC6Z']
# Now let's test a number of basis sets
aug = ['aug-cc-pVDZ', 'aug-cc-pVTZ', 'aug-cc-pVQZ', 'aug-cc-pV5Z', 'aug-cc-pV6Z']
aug_core_polarized = ['aug-cc-pCVDZ', 'aug-cc-pCVTZ', 'aug-cc-pCVQZ', ]
daug = ['d-aug-cc-pVDZ', 'd-aug-cc-pVTZ', 'd-aug-cc-pVQZ', 'd-aug-cc-pV5Z', 'd-aug-cc-pV6Z']
daug_core_polarized = ['d-aug-cc-pCVDZ', 'd-aug-cc-pCVTZ', 'd-aug-cc-pCVQZ', ]

all_basis = aug + aug_core_polarized + daug + daug_core_polarized + no_aug + no_aug_core_polarized

na_data = {}
cl_data = {}
for basis in all_basis:
    try:
        na_ionization = compute_basis_ionization_potential('Na', 'hf', basis)
        na_data[basis] = na_ionization
    except KeyError as k:
        print("Key error: ", k)
        na_data[basis] = np.NaN
    except TypeError as t:
        print("Type error: ", t)
        na_data[basis] = np.NaN
    try:
        cl_electron_affinity = compute_basis_electron_affinity('Cl', 'hf', basis)
        cl_data[basis] = cl_electron_affinity
    except KeyError as k:
        print("Key error: ", k)
        cl_data[basis] = np.NaN
    except TypeError as t:
        print("Type error: ", t)
        cl_data[basis] = np.NaN

# now create a pandas dataframe

na_df = pd.DataFrame.from_dict(na_data, orient='index', columns=['ionization_potential'])
# add molecule column
na_df['molecule'] = 'Na'
# name the index 'basis'
na_df.index.name = 'basis'

cl_df = pd.DataFrame.from_dict(cl_data, orient='index', columns=['electron_affinity'])
cl_df['molecule'] = 'Cl'
cl_df.index.name = 'basis'
# reset the index
na_df.reset_index(inplace=True)
cl_df.reset_index(inplace=True)

# make a detailed dataframe


# now append madness data to the dataframe
na_df['madness'] = na_ionization
cl_df['madness'] = cl_electron_affinity

# compute the error
na_df['error'] = na_df['ionization_potential'] - na_df['madness']
cl_df['error'] = cl_df['electron_affinity'] - cl_df['madness']

print(na_df.query("basis == 'd-aug-cc-pVDZ'"))

na_df = make_detailed_df(na_df)
cl_df = make_detailed_df(cl_df)

# reorder the CV and V categories
print(na_df['Type'].unique())

# na_df['Type'] = na_df['Type'].cat.reorder_categories(
#    ['aug-cc-pVnZ', 'aug-cc-pCVnZ', 'd-aug-cc-pVnZ', 'd-aug-cc-pCVnZ'])
# cl_df['Type'] = cl_df['Type'].cat.reorder_categories(
#    ['aug-cc-pVnZ', 'aug-cc-pCVnZ', 'd-aug-cc-pVnZ', 'd-aug-cc-pCVnZ'])
# now seaborn relplot

# set the palette
pal = sns.color_palette("seismic_r", 4).as_hex()
p1 = [pal[1], pal[0], pal[2], pal[3]]
pal = sns.color_palette(p1)

pal2 = sns.color_palette("coolwarm_r", 4).as_hex()
p2 = [pal2[1], pal2[0], pal2[2], pal2[3]]
p3 = [pal2[1], pal2[2], ]
light_pal = sns.color_palette(p2)
simple_pal = sns.color_palette(p3)

sns.set_theme('paper', 'darkgrid', palette=light_pal, font='sans-serif', font_scale=1.0)

# 1*2 subfigure

r_plot, axes = plt.subplots(1, 2, figsize=(9, 4), frameon=True, layout='constrained',
                            sharex=False, sharey=False)
if type is None:
    pal = light_pal
else:
    pal = simple_pal

sns.lineplot(data=na_df, x='valence', y='error', hue='Type', style='Type', ax=axes[0], markers=True,
             markersize=10,
             palette=light_pal)
sns.lineplot(data=cl_df, x='valence', y='error', hue='Type', style='Type', ax=axes[1], markers=True,
             markersize=10, palette=light_pal)

pd.set_option('display.max_columns', None)

# now lets fix the titles and labels
axes[0].set_title('Na Ionization Potential')
axes[1].set_title('Cl Electron Affinity')

axes[0].set_ylabel('Error (Hartree)')
axes[0].set_xlabel('Valence')
axes[1].set_xlabel('Valence')
axes[1].set_ylabel('Error (Hartree)')

# now put this in my thesis


thesis_path = Path('/home/adrianhurtado/projects/writing/thesis2023/Figures_v2')
r_plot.savefig(thesis_path.joinpath("nacl_ionization_plot.png"), )
