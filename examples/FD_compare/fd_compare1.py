import pandas as pd
from matplotlib import pyplot as plt
from pathlib import Path
from quantumresponsepro import MadnessResponse


def get_mad_series(data_dir, molecules):
    mad_zz_data = {}
    for mol in molecules:
        mad_mol = MadnessResponse(mol, 'hf', 'dipole', data_dir)
        mad_zz_data[mol] = mad_mol.polar_data.iloc[0]['zz']
    return pd.Series(mad_zz_data, name=data_dir.name)


fd_base = Path("/mnt/data/madness_data/fd_compare")


def get_fd_data(fd_base):
    fd_low_path = fd_base.joinpath("low-low")
    fd_highlow_path = fd_base.joinpath("high-low")
    fd_high_path = fd_base.joinpath("high-high")

    fd_data_path = Path("/home/adrianhurtado/projects/writing/mra-tdhf-polarizability/FD_data.csv")

    # read Finite Difference data from csv

    with open(fd_data_path, 'r') as f:
        fd_data = pd.read_csv(f)

    print(fd_data)

    molecules = fd_data['Molecule']
    print(molecules)

    lowlow = get_mad_series(fd_low_path, molecules)
    highlow = get_mad_series(fd_highlow_path, molecules)
    highhigh = get_mad_series(fd_high_path, molecules)

    mad_data = pd.concat([lowlow, highlow, highhigh], axis=1)
    # rename index to Molecules and reset_index
    mad_data.index.name = 'Molecule'
    mad_data.reset_index(inplace=True)
    mad_data = pd.concat([mad_data, fd_data.set_index('Molecule').reset_index(drop=True)], axis=1)

    # compute the percent error of madness data compared to finite difference data
    # FD data is in column FD
    # madness data is in other columns
    mad_data['Percent Error (low-low)'] = (mad_data['low-low'] - mad_data['FD']) / mad_data[
        'FD'] * 100
    mad_data['Percent Error (high-low)'] = (mad_data['high-low'] - mad_data['FD']) / mad_data[
        'FD'] * 100
    mad_data['Percent Error (high-high)'] = (mad_data['high-high'] - mad_data['FD']) / mad_data[
        'FD'] * \
                                            100

    # reformat the molecules to have \ce{} around them
    mad_data_out = mad_data.copy()
    mad_data['Molecule'] = mad_data['Molecule'].apply(lambda x: "\\ce{" + x + "}")

    print(mad_data)
    # create a copy with just the percent error columns and FD column
    mad_data_percent_error = mad_data[
        ['Molecule', 'Percent Error (low-low)', 'Percent Error (high-low)',
         'Percent Error (high-high)', 'FD']].copy()
    # rename the columns to be just low-low, high-low, high-high
    mad_data_percent_error.rename(columns={'Percent Error (low-low)': 'low-low',
                                           'Percent Error (high-low)': 'high-low',
                                           'Percent Error (high-high)': 'high-high'}, inplace=True)
    # set the index to be the molecule
    mad_data_percent_error.set_index('Molecule', inplace=True, drop=True)
    return mad_data_percent_error, mad_data_out


# print to latex table printing 6 significant figures for madness data and scientific notation
# for percent errors by creating pandas styler objects and then using the to_latex method

mad_data_percent_error, mad_data = get_fd_data(fd_base)
mad_data_styled = mad_data_percent_error.style.format({'low-low': '{:.2e}',
                                                       'high-low': '{:.2e}',
                                                       'high-high': '{:.2e}', })

mad_data_styled.to_latex(
    Path("/home/adrianhurtado/projects/writing/mra-tdhf-polarizability/FD_data.tex"),
    siunitx=True)
print(mad_data_styled)
fd_base = Path("/mnt/data/madness_data/fd_compare2")
mad_data_percent_error, mad_data_2 = get_fd_data(fd_base)

# reorder to columns to place the FD data first and then the madness data
mad_data_percent_error = mad_data_percent_error[['FD', 'low-low',
                                                 'high-high']]

# rename FD columns to \alpha_{zz}(FD)
mad_data_percent_error.rename(columns={'FD': '$\\alpha_{zz}(FD)$'}, inplace=True)
# rename low-low to MRA low
mad_data_percent_error.rename(columns={'low-low': 'MRA low'}, inplace=True)
# rename high-high to MRA high
mad_data_percent_error.rename(columns={'high-high': 'MRA high'}, inplace=True)

column_format = "|c|ccc|c|"
mad_data_styled = mad_data_percent_error.style.format({'MRA low': '{:.2e}',
                                                       'MRA high': '{:.2e}', })

mad_data_styled.to_latex(
    Path("/home/adrianhurtado/projects/writing/thesis2023/Tables/FD_data2.tex"),
    siunitx=True, column_format=column_format, multicol_align='|c|', hrules=True, )
print(mad_data_styled)

# for each row plot a line at the value of the polarizability and then the values of low-low,
# high-low, and high-high + values of finite difference
# i'd like to add fd data to each value in high-low, high-high, and low-low

fig = plt.figure()
ax = fig.add_subplot(111)
for i in range(len(mad_data)):
    row = mad_data_2.iloc[i]
    ax.plot([0, 1, 2, ], [row['Percent Error (low-low)'], row['Percent Error (high-low)'],
                          row['Percent Error (high-high)']],
            label=row['Molecule'])

    ax.set_xticks([0, 1, 2, ])
    ax.set_xticklabels(['low-low', 'high-low', 'high-high'])
    ax.set_ylabel('Polarizability (a.u.)')

    ax.legend()
ax.axhline(0, color='black')
plt.show()
