import json
from quantumresponsepro import DaltonRunner
import os

from matplotlib import pyplot as plt
from quantumresponsepro import MadnessResponse, BasisMRADataAnalyzer
from quantumresponsepro.BasisMRADataCollection import get_polar_df, BasisMRADataCollection
from pathlib import Path
import pandas as pd
import re
from quantumresponsepro import DatabaseGenerator


def get_mad_series(data_dir, molecules):
    mad_zz_data = {}
    for mol in molecules:
        mad_mol = MadnessResponse(mol, 'hf', 'dipole', data_dir)
        mad_zz_data[mol] = mad_mol.polar_data.iloc[0]['zz']
    return pd.Series(mad_zz_data, name=data_dir.name)


fd_base = Path("/mnt/data/madness_data/fd_compare3")
db_gen = DatabaseGenerator(fd_base)
molecules = ['LiH', 'BH', 'FH', 'CO', 'BF', 'NaCl', ]

fd_data = pd.read_csv(fd_base.joinpath('FD_compare_2.csv'), dtype={'alpha': float,
                                                                   'beta': float})
print(fd_data)

low_path = fd_base.joinpath('low-low')
high_path = fd_base.joinpath('high-high')

low_low = get_mad_series(fd_base.joinpath('low-low'), molecules)
high_high = get_mad_series(fd_base.joinpath('high-high'), molecules)

print(low_low)
print(high_high)

# grab
# only grab data for these molecules
fd_data = fd_data[fd_data['System'].isin(molecules)]
fd_data_copy = fd_data.copy()
print(fd_data)


def get_mad_beta_zz(path):
    beta_mad = {}
    for mol in molecules:
        beta_json_path = path.joinpath('hf/{}/beta.json'.format(mol))
        loaded_json = json.loads(beta_json_path.read_text())
        beta_mad[mol] = beta_json = \
            pd.DataFrame(loaded_json).query('A =="Z" & B=="Z" & C=="Z" ').Beta.values[0]
    # set the index to the molecules
    mzz = pd.Series(beta_mad, name=path.name)
    # name the index to 'System'
    mzz.index.name = 'System'

    return mzz


mad_beta_low = get_mad_beta_zz(fd_base.joinpath('low-low'))
mad_beta_high = get_mad_beta_zz(fd_base.joinpath('high-high'))

print(mad_beta_low)
print(mad_beta_high)

# append low and high beta values to fd_data
fd_data = pd.concat([fd_data.set_index('System'), low_low.rename('alpha-low'),
                     high_high.rename('alpha-high'),
                     mad_beta_low.rename('beta-low'),
                     mad_beta_high.rename(
                         'beta-high')],
                    axis=1)

# drop dipole column
fd_data.drop('dipole', axis=1, inplace=True)

# take the absolute value of all values
fd_data = fd_data.abs()
print(fd_data.keys())
# split the alpha and beta values into separate dataframes
fd_data_alpha = fd_data[['R ', 'alpha', 'alpha-low', 'alpha-high']].copy()
fd_data_beta = fd_data[['beta', 'beta-low', 'beta-high']].copy()
# compute the relative error of the alpha and beta values compared to the finite difference values


fd_data_alpha['alpha-low'] = (fd_data_alpha['alpha-low'] - fd_data_alpha['alpha']) / fd_data_alpha[
    'alpha'] * 100
fd_data_alpha['alpha-high'] = ((fd_data_alpha['alpha-high'] - fd_data_alpha['alpha']) /
                               fd_data_alpha['alpha']) * 100
fd_data_beta['beta-low'] = (fd_data_beta['beta-low'] - fd_data_beta['beta']) / fd_data_beta[
    'beta'] * 100
fd_data_beta['beta-high'] = (fd_data_beta['beta-high'] - fd_data_beta['beta']) / fd_data_beta[
    'beta'] * 100
print(fd_data)
print(fd_data_alpha)
print(fd_data_beta)

# now to stylize the tables for latex
thesis_path = Path('/home/adrianhurtado/projects/writing/thesis2023')
paper_path = thesis_path
# first we put \ce{} around the molecules
fd_data_alpha.index = fd_data_alpha.index.map(lambda x: "\\ce{" + x + "}")
fd_data_beta.index = fd_data_beta.index.map(lambda x: "\\ce{" + x + "}")
# now we style the tables
# place the beta data next to alpha data
fd_data_alpha = pd.concat([fd_data_alpha, fd_data_beta], axis=1)

# change R index int R (a.u.)
fd_data_alpha.rename(columns={'R ': 'R (a.u.)'}, inplace=True)
# change the column name of alpha and beta to \alpha and \beta
fd_data_alpha.rename(columns={'alpha': '$\\alpha_{zz}(FD)$'}, inplace=True)
fd_data_alpha.rename(columns={'beta': '$\\beta_{zzz}(FD)$'}, inplace=True)
# change the column name of alpha-low to 'MRA low' and alpha-high to 'MRA high'
fd_data_alpha.rename(columns={'alpha-low': '$\\alpha_{zz}(low)$'}, inplace=True)
fd_data_alpha.rename(columns={'alpha-high': '$\\alpha_{zz}(high)$'}, inplace=True)
# change the column name of beta-low to 'MRA low' and beta-high to 'MRA high'
fd_data_alpha.rename(columns={'beta-low': '$\\beta_{zzz}(low)$'}, inplace=True)
fd_data_alpha.rename(columns={'beta-high': '$\\beta_{zzz}(high)$'}, inplace=True)
# in the latex i'd like to place a verticle line between alpha and beta columns

# reorder the columns to place the R column first then alpha and beta column for FD and the low and high columns for MRA
fd_data_alpha = fd_data_alpha[
    ['R (a.u.)', '$\\alpha_{zz}(FD)$', '$\\beta_{zzz}(FD)$', '$\\alpha_{zz}(low)$',
     '$\\beta_{zzz}(low)$', '$\\alpha_{zz}(high)$', '$\\beta_{zzz}(high)$']]

# now we can write the latex tables

column_format = "|c|c|cc|cc|cc|"

# format R column to 4 decimal places
# fd_data_alpha['R (a.u.)'] = fd_data_alpha['R (a.u.)'].map(lambda x: "{:.4f}".format(x))

fd_data_alpha.style.format(
    {'$\\alpha_{zz}(high)$': "{:.2e}", '$\\beta_{zzz}(high)$': "{:.2e}",
     '$\\beta_{zzz}(low)$': "{:.2e}",
     "$\\alpha_{zz}(low)$": "{:.2e}", "R (a.u.)": "{:.4f}",
     }).to_latex(
    paper_path.joinpath(
        Path(
            "Tables/fd_compare_2.tex")),
    convert_css=True,
    multicol_align='|c|',
    hrules=True,
    siunitx=True,
    column_format=column_format, )
# insert the R values before the alpha values
print(fd_data_alpha)

# now read the basis set data for the molecules
basis_sets = ['aug-cc-pVDZ', 'aug-cc-pVTZ', 'aug-cc-pVQZ'] + ['d-aug-cc-pVDZ', 'd-aug-cc-pVTZ',
                                                              'd-aug-cc-pVQZ'
                                                              ] + ['d-aug-cc-pCVTZ',
                                                                   'd-aug-cc-pCVQZ']

runner = DaltonRunner(fd_base, False)

alpha_dict = {}
beta_dict = {}
for mol in molecules:
    alpha_dict[mol] = {}
    beta_dict[mol] = {}

    for basis in basis_sets:

        try:
            result = runner.get_quad_json(mol, 'hf', 'dipole', basis)

            result.rename(columns={'A-freq': 'Afreq', 'B-freq': 'Bfreq', 'C-freq': 'Cfreq'},
                            inplace=True)

            result.rename(columns={'Beta Value': 'Beta'}, inplace=True)
            beta_zzz = result.query('A =="Z" & B=="Z" & C=="Z" ').query(
                'Afreq==0.0 & '
                'Bfreq==0.0 & Cfreq ==0.0 ').Beta.values[0]
            polar_data = runner.get_polar_json(mol, 'hf', 'dipole', basis)
            alpha_zz = polar_data[basis]['response']['values']['zz'][0]
            alpha_dict[mol][basis] = alpha_zz
            beta_dict[mol][basis] = beta_zzz


        except (FileNotFoundError, TypeError) as f:
            print(f)
            pass
basis_alpha = pd.DataFrame(alpha_dict).T
basis_beta = pd.DataFrame(beta_dict).T

# now we compare the basis set data to the FD data by computing the relative error
# for each column i'd like to compute the percent error with respect to the FD data column



fd_data_copy.set_index('System', inplace=True)

fd_alpha = fd_data_copy['alpha'].copy()
fd_beta = fd_data_copy['beta'].copy()


# concat the basis set data with the fd data
basis_alpha = pd.concat([fd_alpha, basis_alpha], axis=1)
basis_beta = pd.concat([fd_beta, basis_beta], axis=1)




def compute_percent_error(df, target):
    df = df.copy().abs()
    for col in df.columns[1:]:  # Skip the first column
        percent_error = (df[col] - df[target]) / df[target] * 100
        df[f'{col}'] = percent_error
    return df


pe_alpha = compute_percent_error(basis_alpha, 'alpha')
pe_beta = compute_percent_error(basis_beta, 'beta')

# now we can write the latex tables
# first we put \ce{} around the molecules
pe_alpha.index = pe_alpha.index.map(lambda x: "\\ce{" + x + "}")
pe_beta.index = pe_beta.index.map(lambda x: "\\ce{" + x + "}")
# now we style the tables

# in this case we drop the alpha and beta columns
pe_alpha.drop(['alpha'], axis=1, inplace=True)
pe_beta.drop(['beta'], axis=1, inplace=True)

# now we format all the colums to 2 decimal places
pe_alpha = pe_alpha.applymap(lambda x: "{:.2f}".format(x))
pe_beta = pe_beta.applymap(lambda x: "{:.2f}".format(x))

# now we can write the latex tables
column_format = "|c|ccc|ccc|ccc|ccc|"

# rotate the column names by 45 degrees


# now we can write the latex tables
pe_alpha.style.applymap_index(lambda v: "rotatebox:{45}--rwrap--latex;", axis=1).to_latex(
    paper_path.joinpath(
        Path(
            "Tables/alpha_basis_compare_fd.tex")),
    convert_css=True,
    multicol_align='|c|',
    hrules=True,
    siunitx=True,
    column_format=column_format, )

pe_beta.style.applymap_index(lambda v: "rotatebox:{45}--rwrap--latex;", axis=1).to_latex(
    paper_path.joinpath("Tables/beta_basis_compare_fd.tex"),
    convert_css=True,
    multicol_align='|c|',
    hrules=True,
    siunitx=True,
    column_format=column_format, )

# compute the percent error of the basis set data compared to the fd data
