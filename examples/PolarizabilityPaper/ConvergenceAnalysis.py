import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from pathlib import Path
from quantumresponsepro import MadnessResponse

xc='hf'
fd_data_path = Path("/home/adrianhurtado/projects/writing/mra-tdhf-polarizability/FD_data.csv")

# read Finite Difference data from csv

with open(fd_data_path, 'r') as f:
    fd_data = pd.read_csv(f)

fd_beta_data_path = Path("/mnt/data/madness_data/fd_compare3/FD_compare_2.csv")
with open(fd_beta_data_path, 'r') as f:
    fd_beta_data = pd.read_csv(f)



def get_z_residuals(mol, path):
    m = MadnessResponse(mol, xc, 'dipole', path)
    freq_key = m.frequencies[0]
    rd = m.data['convergence'][freq_key]['density_residuals'].rename(
        columns={'x': 'dx', 'y': 'dy', 'z': 'dz'})
    try:
        rx = m.data['convergence'][freq_key]['x_relative_residuals'].loc[:, ['x', 'y', 'z']].rename(
            columns={'x': 'rx', 'y': 'ry', 'z': 'rz'})
    except KeyError as e:
        rx = m.data['convergence'][freq_key]['x_residuals'].loc[:, ['x', 'y', 'z']].rename(
            columns={'x': 'rx', 'y': 'ry', 'z': 'rz'})

    x = (m.data['convergence'][freq_key]['x'].loc[:, ['xx', 'yy', 'zz']].rename(columns={'xx':
                                                                                             'xnorm',
                                                                                         'yy': 'ynorm',
                                                                                         'zz':
                                                                                             'znorm'}))
    x = np.sqrt(x)
    alpha = m.data['convergence'][freq_key]['alpha'].loc[:, ['zz']]
    iters = m.num_iter_proto[m.frequencies[0]]
    iters = [i - 1 for i in iters]

    residuals = pd.concat(
        [rx.iloc[iters, 2], rd.iloc[iters, 2], x.iloc[iters, 2], alpha.iloc[iters]],
        axis=1)
    # reindex as low, medium, high
    # if index is size 2 then its just low and medium
    if len(residuals.index) == 2:
        residuals.index = ['Low', 'Medium']
    else:
        residuals.index = ['Low', 'Medium', 'High']
    # rename columns as max bsh and max density
    residuals.rename(columns={0: 'z bsh', 1: 'z density', 'zz': 'MRA'}, inplace=True)
    # add a column for the FD value of the molecule
    residuals['FD'] = fd_data.loc[fd_data['Molecule'] == mol, 'FD'].values[0]
    # get the percent error
    residuals['Percent Error'] = (residuals['MRA'] - residuals['FD']) / residuals['FD'] * 100
    residuals['Absolute Percent Error'] = residuals['Percent Error'].abs()
    residuals['Absolute Residual Z'] = residuals['rz'] * residuals['znorm']
    # add a column for molecule
    residuals['Molecule'] = mol
    return residuals


def get_database_residuals_z(mols, path):
    residuals = pd.DataFrame()
    for mol in mols:
        print(mol)
        try:
            residuals = pd.concat([residuals, get_z_residuals(mol, path)], axis=0)
        except (KeyError, FileNotFoundError) as e:
            print(e)
            continue
    # reset the index
    residuals.reset_index(inplace=True, drop=False)
    # rename the index to be protocol
    residuals.rename(columns={'index': 'Protocol'}, inplace=True)
    return residuals

    # plot the linear regression fit on the same axis
