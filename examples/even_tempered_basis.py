import os

import pandas as pd
from matplotlib import pyplot as plt

from quantumresponsepro import EvenTemperedBasis
from pathlib import Path
from quantumresponsepro import MadnessResponse
from quantumresponsepro import DaltonRunner


# Path: examples/even_tempered_basis.py


def get_alpha_beta(n, a1, an):
    alpha = a1
    beta = (an / a1) ** (1 / (n - 1))
    alphas = [alpha * beta ** (i - 1) for i in range(1, n + 1)]
    return alphas


etb = EvenTemperedBasis('/home/adrianhurtado/projects/dalton/build/basis')
database_path = Path("/mnt/data/madness_data/post_watoc/august/")

dalton = DaltonRunner(database_path, True)

element_betas = {'H': (1776, .005), 'He': (4785, .04), 'O': (570800, .05), 'Be': (12000, .05)}
hydrogen_betas = {'H': {'s': (1775, .005), 'p': (8.6, .005), 'd': (4.5, .005), 'f': (2.0,
                                                                                     .005),
                        'g': (1, .005), 'h': (.8, .005)}}

# make a dictionary of number of elements for each orbital s,p,d,f,g,h,.. starting at 8 for s and
# 7 for p and so on

# The idea for the loop is to start with a small basis and then add elements to the basis until
# we reach convergence.  Convergence is measured by the difference in values between consecutive
# iterations.  The values we are going to track is the ground state energy and the polarizability
# in xx and zz directions.  We will converge each orbital type separately.  We will start with s
# orbitals and then move to p orbitals and so on.  We start with the value in the dictionary and
# then we add one element at a time until we reach convergence.  Once we converge an orbital type
# we move to the next one, keeping the previous ones within the basis.  We will stop when we
# reach the end of the dictionary.  During each iteration we will also track the value of the
# energy and the polarizability in xx and zz directions.  Also, it is important to tighten the
# convergence threshold as we move to higher orbitals since the convergence is expected to be
# smaller.

orb_dict = {'s': 4, 'p': 7, 'd': 6, 'f': 5, 'g': 4, 'h': 3, }
orb_dict = {'s': 10, 'p': 10, 'd': 10, 'f': 3, 'g': 3, 'h': 3}
mol = 'H2'
en = []

active_orbs = []
max_iter_per_orb = 10
e_res = 1e-4
a_res = .001

data = {'energy': [], 'alpha_xx': [], 'alpha_yy': [], 'alpha_zz': [], 'energy_residual': [],
        'alpha_xx_residual': [], 'alpha_yy_residual': [], 'alpha_zz_residual': []}
en = 0
gi = 10000
xx_i = 10000
yy_i = 10000
zz_i = 10000

e_res_i = 1.0
xx_res_i = 1.0
zz_res_i = 1.0
orb_j = 0

all_orbs = ['s', 'p', 'd', 'f', 'g', 'h']

for orb_j, orb in enumerate(all_orbs, 1):
    active_orbs = all_orbs[:orb_j + 1]
    i = 0
    print(i, orb, active_orbs)

    elem_dict = {}
    orb_exponents_dict = {}

    while (i < max_iter_per_orb) and e_res_i > e_res and xx_res_i > a_res and zz_res_i > a_res:
        print(e_res, e_res_i)

        # loop through active orbitals and generate exponents for each based on value in orb_dict
        for orb_a in active_orbs:
            orb_max_exponent = hydrogen_betas['H'][orb_a][0] * 1.5
            orb_min_exponent = hydrogen_betas['H'][orb_a][1] / 1.5
            orb_num = orb_dict[orb_a]
            exponents = get_alpha_beta(orb_num, orb_max_exponent, orb_min_exponent)
            orb_exponents_dict[orb_a] = exponents
        elem_dict['H'] = orb_exponents_dict
        etb.write_basis('test', elem_dict)

        filePath = database_path.joinpath(
            'dalton/hf/{}/dipole/freq_{}-EVT-test.out'.format(mol, mol))
        if os.path.exists(filePath):
            os.remove(filePath)
        ground, response = dalton.get_frequency_result(mol, 'hf', 'dipole', 'EVT-test')
        gn = ground['totalEnergy']
        xx = response.loc[0, 'xx']
        yy = response.loc[0, 'yy']
        zz = response.loc[0, 'zz']
        data['energy'].append(gn)
        data['alpha_xx'].append(xx)
        data['alpha_yy'].append(yy)
        data['alpha_zz'].append(zz)
        orb_dict[orb] += 1
        i += 1
        # print last iteration data
        print("iteration: ", i, "energy: ", gn, "alpha_xx: ", xx, "alpha_yy: ", yy, "alpha_zz: ",
              zz)
        e_res_i = abs(gn - gi)
        ax_res_i = abs(xx - xx_i)
        ay_res_i = abs(xx - yy_i)
        az_res_i = abs(zz - zz_i)

        # line printing the change in ground state energy
        print("energy residual: ", e_res_i, "alpha_xx residual: ", ax_res_i, "alpha_zz residual: "
              , az_res_i)
        data['energy_residual'].append(e_res_i)
        data['alpha_xx_residual'].append(ax_res_i)
        data['alpha_yy_residual'].append(ay_res_i)
        data['alpha_zz_residual'].append(az_res_i)
        gi = gn
        xx_i = xx
        yy_i = yy
        zz_i = zz
    e_res /= 5.0
    a_res /= 5.0
    # set the max iteration per orbital to be difference between orbital i and orbital i+1
    if orb_j < len(orb_dict) - 1:
        max_iter_per_orb = abs(orb_dict[orb] - orb_dict[list(orb_dict.keys())[orb_j + 1]])
    orb_j += 1

# plot energy and polarizability as a function of iterations


data = pd.DataFrame(data)
data.to_csv('h2_evp.csv')

fig = plt.figure()
plt.plot(data['energy'][len(data) - 10:])
plt.show()
fig = plt.figure()
# only plot the last 10
plt.plot(data['alpha_xx'][len(data) - 10:])
plt.plot(data['alpha_yy'][len(data) - 10:])
plt.plot(data['alpha_zz'][len(data) - 10:])
plt.show()

fig = plt.figure()
plt.plot(data['energy_residual'])
plt.title('Energy Residual')
plt.yscale('log')
plt.show()

fig = plt.figure()
plt.plot(data['alpha_xx_residual'])
plt.plot(data['alpha_yy_residual'])
plt.plot(data['alpha_zz_residual'])
plt.title('Polarizability Residual')
plt.yscale('log')
plt.show()

basis_data = pd.DataFrame(orb_dict)
