import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from pathlib import Path
from quantumresponsepro import MadnessResponse
from DataAnalysisClass import make_detailed_df

xc='hf'
fd_data_path = Path("/home/ahurta92/Projects/Latex/mra-tdhf-polarizability/csv/FD_data.csv")

# read Finite Difference data from csv

with open(fd_data_path, 'r') as f:
    fd_data = pd.read_csv(f)

#fd_beta_data_path = Path("/mnt/data/madness_data/fd_compare3/FD_compare_2.csv")
#with open(fd_beta_data_path, 'r') as f:
#    fd_beta_data = pd.read_csv(f)



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
        residuals.index = ['Guess', 'Low']
    else:
        residuals.index = ['Guess', 'Low', 'High']
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
single = ['aug-cc-pVDZ', 'aug-cc-pVTZ', 'aug-cc-pVQZ', 'aug-cc-pV5Z', 'aug-cc-pV6Z']
single_polarized = ['aug-cc-pCVDZ', 'aug-cc-pCVTZ', 'aug-cc-pCVQZ']
double = ['d-aug-cc-pVDZ', 'd-aug-cc-pVTZ', 'd-aug-cc-pVQZ', 'd-aug-cc-pV5Z', 'd-aug-cc-pV6Z']
double_polarized = ['d-aug-cc-pCVDZ', 'd-aug-cc-pCVTZ', 'd-aug-cc-pCVQZ']
all_basis_sets = single + single_polarized+ double  + double_polarized


def create_z_data(data,fd_version='one'):
    z_data = data.polar_data.query('ij=="zz" & omega==0').copy()

    if fd_version=='one':
        z_data['FD'] = z_data.apply(lambda x: fd_data.query('Molecule=="{}"'.format(x.molecule))['FD']
                                    .values[0], axis=1)
    else:
        z_data['FD'] = z_data.apply(lambda x: fd_beta_data.query('Molecule=="{}"'.format(x
                                                                                         .molecule))['alpha']
                                    .values[0], axis=1)
    # compute the percent error relative to the FD value
    z_data['Percent Error'] = z_data.apply(lambda x: (x['alpha'] - x['FD']) / x['FD'] * 100, axis=1)
    # create a column which includes the MRA Percent error for each molecule
    z_data['MRA Percent Error'] = z_data.apply(
        lambda x: z_data.query('molecule=="{}"'.format(x.molecule))
        ['Percent Error'].values[0], axis=1)
    z_data['Absolute Percent Error'] = z_data.apply(lambda x: abs(x['Percent Error']), axis=1)
    # make a column to indicate basis is MRA or orther
    z_data['MRA'] = z_data.apply(lambda x: 'Other' if x['basis'] in all_basis_sets else 'MRA',
                                 axis=1)

    z_data = z_data.query('omega==0')
    return z_data


def compare_z_to_basis_set(z_data, y='Percent Error'):
    q_data = z_data.query('basis in @all_basis_sets & omega==0').copy()
    basis_data = make_detailed_df(q_data).query('omega==0').copy()
    # make the color palette go from light red to dark blue and a yellow for the MRA
    colors = sns.color_palette("muted", 4)
    # sns.set(rc={"xtick.bottom": True, "ytick.left": True}, font_scale=1.5)
    pal = sns.color_palette("seismic_r", 4).as_hex()
    p1 = [pal[1], pal[0], pal[2], pal[3]]
    pal = sns.color_palette(p1)

    # increase the jitter to make the points more legible
    g = sns.relplot(x='molecule', y=y, data=basis_data.query('omega==0'),
                    style='Type',
                    hue='valence',
                    palette=colors, s=150,
                    height=5, aspect=1.5, alpha=.5)

    # set y axis to be log scale
    for ax in g.axes.flat:
        if y == 'Percent Error':
            ax.set_yscale('symlog', linthresh=1e-3)
        else:
            ax.set_yscale('log')
        # add the MRA reference line for this data
        ax.scatter(z_data.query('MRA=="MRA" & Protocol=="Low"')['molecule'],
                   z_data.query('MRA=="MRA" & Protocol=="Low"')[y],
                   color='black', s=100, label='MRA', marker='^')

        ax.scatter(z_data.query('MRA=="MRA" & Protocol=="High"')['molecule'],
                   z_data.query('MRA=="MRA" & Protocol == "High"')[y],
                   color='black', s=50, label='MRA', marker='o')

        ax.set_xlabel('Molecule')
    return g

def create_z_beta_data(beta_data):

    z_data=beta_data.query('ijk=="ZZZ"').copy()
    # Add a column for the FD data
    # add a column for the FD data bases on the molecule
    z_data['FD'] = z_data['molecule'].apply(lambda x: fd_beta_data.query('Molecule == @x')['beta'].values[0])
    # remove HeNe
    z_data=z_data.query('molecule != "HeNe"').copy()
    # remove LiBr
    z_data=z_data.query('molecule != "LiBr"').copy()
    z_data['Percent Error'] = (z_data['Beta'] - z_data['FD']) / z_data['FD'] * 100
    z_data['Absolute Percent Error'] = z_data['Percent Error'].abs()

    z_data['MRA'] = z_data.apply(lambda x: 'Other' if x['basis'] in all_basis_sets else 'MRA',
                                 axis=1)
    return z_data

def compare_z_beta_to_basis_set(z_data,y='Percent Error'):

    q_data=z_data.query('basis in @all_basis_sets').copy()
    basis_data=make_detailed_df(q_data).query('Afreq==0').copy()
    # make the color palette go from light red to dark blue and a yellow for the MRA
    colors= sns.color_palette("muted", 4)
    # sns.set(rc={"xtick.bottom": True, "ytick.left": True}, font_scale=1.5)
    pal = sns.color_palette("seismic_r", 4).as_hex()
    p1 = [pal[1], pal[0], pal[2], pal[3]]
    pal = sns.color_palette(p1)

    # increase the jitter to make the points more legible
    g=sns.relplot(x='molecule', y=y, data=basis_data.query('Afreq==0'),
                  style='Type',
                  hue='valence',
                  palette=colors, s=250,
                  height=5,aspect=1.5,alpha=.5)

    # set y axis to be log scale
    for ax in g.axes.flat:
        if y=='Percent Error':
            ax.set_yscale('symlog', linthresh=1e-3)
        else:
            ax.set_yscale('linear')
            ax.set_yscale('log')
        # add the MRA reference line for this data
        # add the mra data to the plot
        ax.scatter(z_data.query('MRA=="MRA" & Protocol=="Low"')['molecule'],
                   z_data.query('MRA=="MRA" & Protocol=="Low"')[y],
                   color='black', s=200, label='MRA',marker='^')

        ax.scatter(z_data.query('MRA=="MRA" & Protocol=="High"')['molecule'],
                   z_data.query('MRA=="MRA" & Protocol == "High"')[y],
                   color='black', s=50, label='MRA',marker='o')
        # add a label for the MRA data

        ax.set_xlabel('Molecule')
    return g
