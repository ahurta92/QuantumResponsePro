import numpy as np
import pandas as pd

from quantumresponsepro.BasisMRADataAssembler import get_mra_quad_data
from quantumresponsepro import BasisMRAData, DaltonRunner
from quantumresponsepro import MadnessResponse
from pathlib import Path
import json
from quantumresponsepro.BasisMRADataCollection import get_quad_df

database_path = Path('/mnt/data/madness_data/august_no_symmetry')

mol = "H2O"
dr = DaltonRunner(database_path, False)
beta_json = dr.get_quad_json(mol, 'hf', 'dipole', 'd-aug-cc-pVDZ')

mols = ['CH3SH', 'H2O', 'NaCl', 'CH3OH', ]
basis_sets = ['aug-cc-pVDZ', 'aug-cc-pVTZ', 'aug-cc-pVQZ', 'd-aug-cc-pVDZ', 'd-aug-cc-pVTZ',
              'd-aug-cc-pVQZ', ]


def get_basis_mra_data_df(mol, basis, database_path, mol_data, Bfreq, Cfreq):
    """
    This function returns
    """
    b_data = mol_data.query('basis == @basis & Bfreq == @Bfreq & Cfreq== @Cfreq').copy()

    # remove the nan values
    b_data.dropna(inplace=True)
    b_data.drop_duplicates(subset=['ijk'], inplace=True)
    # remove Beta less than 1e-4
    m_data = mol_data.query('basis == "MRA" & Bfreq==@Bfreq & Cfreq == @Cfreq ').copy()
    # remove the Beta less than .001
    m_data.dropna(inplace=True)
    # remove the Beta less than 1e-4
    m_data = m_data.query('Beta.abs() > .005').copy()
    BC_index = ['Afreq', 'Bfreq', 'Cfreq', 'ijk']

    # here we deal with the extra frequencies in MRA
    # we need to remove the extra frequencies

    m_data.set_index(BC_index, inplace=True)
    b_data.set_index(BC_index, inplace=True)

    # b_data.drop_duplicates(inplace=True, subset=['ijk'])
    b_data.drop_duplicates()

    b_data['abs'] = b_data['Beta'].abs()
    b_data.sort_values(by=['abs'], inplace=True)

    m_data['abs'] = m_data['Beta'].abs()
    m_data.sort_values(by=['abs'], inplace=True)

    print(m_data)
    b_data.drop_duplicates()
    # drop rows with Beta less than 1e-4
    b_data = b_data.query('Beta.abs() > .005').copy()

    print(b_data)
    m_vals = m_data['Beta'].rename('Beta_MRA').reset_index(drop=True)
    m_vals.index = b_data.index
    print(m_vals)
    # b_data = pd.concat([b_data, m_data['Beta'].rename('Beta_MRA')], axis=1, ignore_index=False)
    b_data['Beta_MRA'] = m_data['Beta']
    b_data['Beta_sign'] = np.sign(b_data['Beta'] * b_data['Beta_MRA'])
    b_data['Beta_MRA_signed'] = b_data['Beta_MRA'] * b_data['Beta_sign']
    b_data['Beta_diff'] = (b_data['Beta'] - b_data['Beta_MRA_signed'])

    return b_data


def compare_dipole_moments(mol, basis, xc, op, database):
    mr = MadnessResponse(mol, xc, op, database)
    mad_dipole = mr.get_dipole_moment()
    print(mad_dipole)
    dr = DaltonRunner(database, False)
    g_data = dr.get_quad_json(mol, xc, op, basis)
    print(g_data)


compare_dipole_moments('H2O', 'd-aug-cc-pVTZ', 'hf', 'dipole', database_path)

mol = 'CH3OH'
q_df = get_quad_df(mols, 'hf', 'dipole', database_path, basis_sets)
q_df_i = q_df.reset_index()
# truncate the Afreq Bfreq and Cfreq to 3 decimal places
q_df_i['Afreq'] = q_df_i['Afreq'].apply(lambda x: round(x, 3))
q_df_i['Bfreq'] = q_df_i['Bfreq'].apply(lambda x: round(x, 3))
q_df_i['Cfreq'] = q_df_i['Cfreq'].apply(lambda x: round(x, 3))
# get the unique frequencies
q_df = q_df_i.copy()
mol_data = q_df.query('molecule == @mol').copy()
freq = list(mol_data['Afreq'].unique())
print(freq)
freq_i = freq[1]

# Comparing the dipole moments to check the geometries of calculations

b_data = get_basis_mra_data_df('CH3SH', 'aug-cc-pVTZ', database_path, mol_data, freq[0],
                               freq[0])

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
print(b_data)
